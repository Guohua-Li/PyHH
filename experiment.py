# modified on May 19, 2023

__doc__ = """
                      PyHH Neuron Simulator

PyHH is a small Python library for simulation of ion channel dynamics and the
electrical activity of cells. It provides two ways of simulating ion channels:
the Hodgkin-Huxley approach and the Markov approach.

Basic units used in PyHH are:
 (1). time:          ms
 (2). voltage:       mV
 (3). current:       pA
 (4). conductance:   nS
 (5). capacitance:   pF
 (6). length:        μm
 (7). area:          μm^2
 (8). concentration: mM
 (9). Specific membrane conductance (Gm):       nS/μm^2
 (10).Specific membrance capacitance (Cm):      pF/μm^2
 (11).Specific axial conductance (Ga):          nS/μ^m

Rate Constant = 1/ms

Symbols:
  transmembrane current density (J_ion): pA/μm^2
  total injected current density (J_injs): pA /μm^2
  maximum conductances per unit area (gMax): nS/μm^2
  axial current density (J_axial):
Some important relationship between units are
 pA = pF * mV / ms = mV / Gohm = nS * mV

 tau = 1/(alpha+beta)
 inf = alpha* tau

 alpha = inf/tau
 beta = (1-inf)/tau

Pre-defined objects:
 (1). NaC:   Sodium channel
 (2). KDR:   potassium channel
 (3). gL:    Leak Conductance


Finally, the Author assumes no responsibility for damage or loss of system
performance as a direct or indirect result of the use of this software.
***************************************************************************


Tree representation of a neuronal structure. Edges are labelled with electrical properties of the represented segment, e.g. functions of axial and membrane resistance/capacitance. Nodes are labelled by point-like current functions, e.g. synaptic input.  



cytoplasm conductivity: the reciprocal of electrical resistivity
or speciﬁc electrical resistance for the cytoplasm modelled as a resistor
    # The new value is taken from 
    # Wang, K., Zhao, Y., Chen, D. et al. Specific membrane capacitance, cytoplasm conductivity and
    # instantaneous Young’s modulus of single tumour cells. Sci Data 4, 170015
    # (2017). https://doi.org/10.1038/sdata.2017.15
    #self.Gx    = 958.0 # average of values from 5 types of cell, the unit is ns/μm

"""

import array
import time

from math import exp, log10
import matplotlib.pyplot as plt#
from matplotlib.gridspec import GridSpec

def isCompartment(x):
  return hasattr(x,'is_Compartment')


class Experiment:
  """
  Integration
  """
  def __init__(self, *preps):
    self.UNITS = []
    """for a in preps:
      if type (a) is list:
        self._add_list(a)
      else:
        self._add_cmp(a)"""
    for item in preps:
      if type (item) is list:
        for b in item: self._add_cmp(b)
          #if isCompartment(b):
          #  self.UNITS.append(b)
          #else:
          #  raise Exception("Input error")
      else:
        self._add_cmp(item)
        """if isCompartment(item):
          self.UNITS.append(item)
        else:
          raise Exception("Input error")"""

    self.ready = False
    self.Clock = 0.0
    self.dt    = 0.005
    for u in self.UNITS:
      #u.Vi = u.V0
      u._calc_surface()

      """
      if u.Length == None: # sphere
        u.Surface = 3.1415926 * (Surface.Diameter**2) # sphere
        for child in u.Children:
          u.Surface -= 3.1415926 * child.Diameter * child.Diameter / 4.0
      else:
        u.Surface = 3.1415926 * u.Diameter * u.Length # surface before handling the ends
        cross_area = 3.1415926 * u.Diameter * u.Diameter / 4.0 # this is the circle area
        if u.Parent:
          if u.Children == []:
            u.Surface = u.Surface + cross_area # This is a terminal segment.
        else:
          u.Surface = u.Surface + cross_area # one end is closed
          if u.Children == []:
            u.Surface = u.Surface + cross_area # the other end is closed"""

      u._calc_gx()

    N = len(self.UNITS)
    print ('_________________________________________\n')
    print ('%i compartment(s) supplied'%(N))

  def _add_list(self, a):
        for b in a:
          if isCompartment(b):
            self.UNITS.append(b)
          else:
            raise Exception("Input error")

  def _add_cmp(self, unit):
        if isCompartment(unit):
          self.UNITS.append(unit)
        else:
          raise Exception("Input error")

  def calc_J_axial(self): # get J_axial (axial current density)
    for u in self.UNITS:
      if u.Parent:
        u.J_axial = (u.Vi - u.Parent.Vi) * u.gx/u.Surface
      else:
        u.J_axial = 0.0

  def calc_J_total(self, i): # net current density to this compartment, excluding J_ion and J_cap
    for u in self.UNITS:
      inj = 0.0 if u.iClamper == None else u.iClamper.Command[i] # clamper current density
      u.J_total = inj - u.J_axial  # current leaving the compartment
      for chd in u.Children:       # current coming from the child compartments
        u.J_total += chd.J_axial

  def _prep(self, steps):
    self.C_Clampers = [] # collect the ligand clamper
    self.V_Clampers = [] # collect the ligand clamper
    self.I_Clampers = [] # collect the ligand clamper
    for u in self.UNITS:
      if u.iClamper: self.I_Clampers.append(u.iClamper)
      if u.vClamper: self.V_Clampers.append(u.vClamper)
      for ch in u.channel_list:
        if ch.Tag == 'LGIC':
          if ch.Ligand.Clamper:
            self.C_Clampers.append(ch.Ligand.Clamper)

    self.VarRates = [] # collect the variable rates (voltage-dependent, so far)
    for u in self.UNITS:
      for var_rate in u.VoltGatedRates: # RateList
        self.VarRates.append(var_rate)

    self.Units_vClamped = [u for u in self.UNITS if u.vClamper]
    self.Units_iClamped = [u for u in self.UNITS if u.iClamper]
    self.Units_Others = list(set(self.UNITS)-set(self.Units_vClamped))
    #Units_CurrentRecorder = [u for u in self.UNITS if u.iRecorder]

    for u in self.UNITS:
      u.dt = self.dt
      u.T = self.T
      u.Vm = array.array('f',[u.V0]) * steps
      for ch in u.channel_list:
        if ch.Recording:
          ch.trace = array.array('f',[0]) * steps

    for u in self.Units_vClamped:
      clamper = u.vClamper
      clamper.T = self.T
      clamper.J_ion  = array.array('f', [0]) * steps
      clamper.J_injs = array.array('f', [0]) * steps
      if clamper.Waveform != None:
        clamper._get_command(steps, self.dt)
        clamper.dot = [0.] + [ (clamper.Command[i]-clamper.Command[i-1])/self.dt for i in range(1,steps) ]
      else:
        clamper.Command = [0.0 + clamper.Baseline for i in range(steps)] # added May 16, 2023
        clamper.dot     = [0.0 for i in range(steps) ]              # added May 16, 2023

    for u in self.Units_iClamped: # generate command from specified waveform
      clamper = u.iClamper
      if clamper.Waveform != None:
        clamper._get_command(steps, self.dt)
      else:
        clamper.Command = [0.0 for i in range(steps)] # added May 16, 2023

    for clamper in self.C_Clampers: # We access the clamper through Experiment
      clamper.T = self.T
      if clamper.Waveform != None:
        clamper._get_command(steps, self.dt)
      else:
        clamper.Command = [0.0 for i in range(steps)] # added May 20, 2023


  def run(self, t, dt=None, verbose=True, reset_clock=True):
    if reset_clock: self.Clock = 0.0
    if dt != None: self.dt = dt
    steps = int(t/self.dt)
    self.T = [i*self.dt + self.Clock for i in range(steps)]
    if verbose:
      print( "Number of steps {}, with step length of {}".format(steps, self.dt) )
      print ("Integration in process ...")

    self._prep(steps)

    ###-------------------------------###

    t0 = time.time()
    for step_num in range(steps): # looping begins

      self.calc_J_axial()
      self.calc_J_total(step_num)

      for c_clamper in self.C_Clampers: # We access the clamper through Experiment
        c_clamper.Ligand.Conc = c_clamper.Command[step_num]

      for var_rate in self.VarRates:
        var_rate.update_rate()

      # Update_Start ========
      for u in self.Units_Others:
        ##print("reached??")
        u._update_Vm(self.dt, step_num)

        for ch in u.channel_list:
          ch._update_gate(u.Vi, self.dt)

      """
      for u in Units_CurrentRecorder: # Current Recording
        i_recorder = u.iRecorder
        i_recorder.J_ion[step_num]  = u.get_Jion(u.Vi)
        i_recorder.J_injs[step_num] = u.J_total
        i_recorder.J_cap[step_num]  = u.Cm * u.Vdot
      """

      for u in self.Units_vClamped: # Voltage-Clamp block
        clamper = u.vClamper
        u.Vi = clamper.Command[step_num]
        for ch in u.channel_list:
          ch._update_gate(u.Vi, self.dt)

        clamper.J_ion[step_num] = u.get_Jion(clamper.Command[step_num], step_num)
        clamper.J_injs[step_num] = u.J_total # renamed

      # Update_End ========
      for u in self.UNITS: # Store membrane potentials
        u.Vm[step_num] = u.Vi

    t1 = time.time()
    print ('Integration done in %f seconds.  ^-^\n' % (t1-t0))
    self.Clock = self.T[-1] + self.dt

    ##########################

  def plot_Vm(self, ylim=[-90,60]):
    xlim = (self.T[0], self.T[-1]+self.dt)
    n = 4
    m = n-1
    fig = plt.figure()
    ax1 = plt.subplot2grid((n,1), (0,0), rowspan=m)
    ax2 = plt.subplot2grid((n,1), (m,0))
    for u in self.UNITS:
      ax1.plot(self.T, u.Vm, linewidth=1.0)
      if u.iClamper:
        ax2.plot(self.T, u.iClamper.Command, linewidth=1.0)

    ax1.set_ylim(ylim)
    ax1.set_xlim(xlim)
    ax1.set_ylabel('mV', color='tab:blue')
    #ax1.tick_params(axis='y', labelcolor='tab:blue')
    ax1.tick_params(axis='x', colors='tab:blue')
    ax1.tick_params(axis='y', colors='tab:blue')
    ax1.spines['left'].set_color('tab:blue')
    ax1.spines['right'].set_color('tab:gray')
    ax1.spines['top'].set_color('tab:gray')
    ax1.spines['bottom'].set_color('tab:gray')

    ax2.set_xlim(xlim)
    ax2.set_xlabel('ms', color='tab:blue')
    ax2.set_ylabel(r'pA/$\mu$$m^2$', color='tab:blue')
    #ax2.tick_params(axis='y', labelcolor='tab:blue')
    ax2.tick_params(axis='x', colors='tab:blue')
    ax2.tick_params(axis='y', colors='tab:blue')
    ax2.spines['left'].set_color('tab:blue')
    ax2.spines['right'].set_color('tab:gray')
    ax2.spines['top'].set_color('tab:gray')
    ax2.spines['bottom'].set_color('tab:gray')

    plt.subplots_adjust(hspace=0.5)
    plt.show()

  def plot_lgic(self):
    xlim = (self.T[0], self.T[-1]+self.dt)
    n = 4
    m = n-1
    fig = plt.figure()
    ax1 = plt.subplot2grid((n,1), (0,0), rowspan=m)
    ax2 = plt.subplot2grid((n,1), (m,0))

    for clamper in self.V_Clampers:
        ax1.plot(clamper.T, clamper.J_ion, linewidth=1.0)
    ax1.set_xlim(xlim)
    ax1.set_ylabel('mV', color='tab:blue')
    ax1.tick_params(axis='y', colors='tab:blue')
    ax1.spines['left'].set_color('tab:blue')
    ax1.spines['right'].set_color('tab:gray')
    ax1.spines['top'].set_color('tab:gray')
    ax1.spines['bottom'].set_color('tab:gray')

    for clamper in self.C_Clampers:
        ax2.plot(clamper.T, clamper.Command, linewidth=1.0)
    ax2.set_xlim(xlim)
    ax2.set_xlabel('ms')
    ax2.set_ylabel('mV', color='tab:blue')
    ax2.tick_params(axis='y', colors='tab:blue')
    ax2.spines['left'].set_color('tab:blue')
    ax2.spines['right'].set_color('tab:gray')
    ax2.spines['top'].set_color('tab:gray')
    ax2.spines['bottom'].set_color('tab:gray')

    plt.subplots_adjust(hspace=0.5)
    plt.show()


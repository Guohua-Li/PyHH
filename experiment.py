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
 (6). length:        um
 (7). area:          um2
 (8). concentration: mM
 (9). Specific membrane conductance (Gm):       nS/um2
 (10).Specific membrance capacitance (Cm):      pF/um2
 (11).Specific axial conductance (Ga):          nS/um

Symbols:
  transmembrane current density (J_ion): pA/um2
  total injected current density (J_injs): pA /um2
  maximum conductances per unit area (gMax): nS/um2
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

"""

import array
import time

from math import exp, log10
import matplotlib.pyplot as plt#


def isCompartment(x):
  return hasattr(x,'is_Compartment')


class Experiment:
  """
  Integration
  """
  def __init__(self, *prep):
    self.UNITS = []
    for a in prep:
      if type (a) is list:
        for b in a:
          if isCompartment(b):
            self.UNITS.append(b)
          else:
            print ("-------------------------------------")
            print ("| Check the input to the experiment |")
            print ("-------------------------------------")
            raise Exception("Input error")
      else:
        if isCompartment(a):
          self.UNITS.append(a)
        else:
          print ("-------------------------------------")
          print ("| Check the input to the experiment |")
          print ("-------------------------------------")
          raise Exception("Input error")

    self.Clock = 0.0
    for cp in self.UNITS:
      #cp.Vi = cp.V0
      cp._calc_surface()
      cp._calc_gx()

    N = len(self.UNITS)
    print ('_________________________________________\n')
    print ('%i compartment(s) supplied'%(N))

  def save(self, filename):
    N = len(self.T)
    f = open(filename,'w')
    for cp in self.UNITS:
      f.write('# ID = %i\n' %(cp.ID))
      for k in range(N):
        f.write('%8.6f\n' %(cp.Vm[k]))
    f.close()

  def endpoint_info(self):
    for cp in self.UNITS:
      cp.show()
      for ch in cp.channel_list: ch.show()

  def calc_axial(self): # get J_axial (axial current density)
    for cp in self.UNITS:
      if cp.Parent:
        cp.J_axial = (cp.Vi - cp.Parent.Vi) * cp.gx/cp.Surface
      else:
        cp.J_axial = 0.0

  def calc_total(self, i): # net current to this compartment, excluding transmembrane ionic and capacitive currents
    for cp in self.UNITS:
      inj = 0.0 if cp.iClamper == None else cp.iClamper.Command[i] # clamper current density
      cp.J_total = inj - cp.J_axial # current leaving the compartment
      for chd in cp.Children:       # current coming from the child compartments
        cp.J_total += chd.J_axial

  def run(self, t, dt=0.005, verbose=True):
    steps = int(t/dt)
    self.T = [i*dt + self.Clock for i in range(steps)]
    if verbose:
      print( "Number of steps {}, with step length of {}".format(steps, dt) )
      print ("Integration in process ...")

    C_Clampers = []
    for cp in self.UNITS: # check for missing values
      for ch in cp.channel_list:
        if ch.Tag == 'LGIC':
          if ch.Ligand.Clamper: C_Clampers.append(ch.Ligand.Clamper)

    self.VoltGatings = []
    for cp in self.UNITS:
      for volt_gate in cp.RateList: # Only VoltageGate
        self.VoltGatings.append(volt_gate)

    ### prepare to run
    vcUnits = [i for i in self.UNITS if i.vClamper]
    icUnits = [i for i in self.UNITS if i.iClamper]
    imUnits = [i for i in self.UNITS if i.iMonitor]
    ncUnits = list(set(self.UNITS)-set(vcUnits))

    for cp in self.UNITS: # check for missing values
      cp.dt = dt
      cp.T = self.T

    for cp in self.UNITS: # prepare the storages based on run parameters
      cp.Vm = array.array('f',[cp.V0]) * steps
      for ch in cp.channel_list:
        ch.gM = array.array('f',[0]) * steps

    for cp in vcUnits:
      cp.vClamper.T = self.T
      cp.vClamper.J_ion  = array.array('f', [0]) * steps
      cp.vClamper.J_injs = array.array('f', [0]) * steps

    for cp in imUnits:
      i_monitor = cp.iMonitor
      i_monitor.J_ion  = array.array('f', [0]) * steps
      i_monitor.J_injs = array.array('f', [0]) * steps
      i_monitor.J_cap  = array.array('f', [0]) * steps # Capacitive current density

    for cp in icUnits: # generate command from specified waveform
      clamper = cp.iClamper
      if clamper.Waveform:
        clamper.Command = [clamper.Waveform._func(i*dt) for i in range(steps)]
      else:
        clamper.Command = [0.0 for i in range(steps)] # added May 16, 2023

    for clamper in C_Clampers: # We access the clamper through Experiment
      clamper.T = self.T

      if clamper.Waveform:
        clamper.Command = [clamper.Waveform._func(i*dt) for i in range(steps)]
      else:
        clamper.Command = [0.0 for i in range(steps)] # added May 20, 2023

    for cp in vcUnits: # generate command from specified waveform
      clamper = cp.vClamper
      if clamper.Waveform:
        clamper.Command = [clamper.Waveform._func(i*dt) + clamper.Baseline for i in range(steps)]
        clamper.dot = [0.] + [ (clamper.Command[i]-clamper.Command[i-1])/dt for i in range(1,steps) ]
      else:
        clamper.Command = [0.0 + clamper.Baseline for i in range(steps)] # added May 16, 2023
        clamper.dot     = [0.0 for i in range(steps) ]              # added May 16, 2023

    ###-------------------------------###

    t0 = time.time()
    for step_num in range(steps): # looping begins

      self.calc_axial()
      self.calc_total(step_num)

      for c_clamper in C_Clampers: # We access the clamper through Experiment
        c_clamper.Ligand.Conc = c_clamper.Command[step_num]

      for volt_gate in self.VoltGatings:
        volt_gate.update()

      # Update_Start ========
      for cp in ncUnits:
        cp._update_Vm(dt)

        for ch in cp.channel_list:
          ch._update_gate(cp.Vi, dt)

      for cp in imUnits: # Current Recording
        i_monitor = cp.iMonitor
        i_monitor.J_ion[step_num]  = cp.get_Jion(cp.Vi)
        i_monitor.J_injs[step_num] = cp.J_total
        i_monitor.J_cap[step_num]  = cp.Cm * cp.Vdot

      for cp in vcUnits: # Voltage-Clamp block
        clamper = cp.vClamper
        cp.Vi = clamper.Command[step_num]
        for ch in cp.channel_list:
          ch._update_gate(cp.Vi, dt)

        clamper.J_ion[step_num] = cp.get_Jion(clamper.Command[step_num])
        clamper.J_injs[step_num] = cp.J_total # renamed

      # Update_End ========
      for cp in self.UNITS: # Store membrane potentials
        cp.Vm[step_num] = cp.Vi
        #for ch in cp.channel_list:
        #  ch.gM[step_num] = ch.gMax

    t1 = time.time()
    print ('Integration done in %f seconds.  ^-^\n' % (t1-t0))
    self.Clock = self.T[-1] + dt

    ##########################

  def plot_Vm(self,vlimit=[-90,60]):
    n = 4
    m = n-1
    fig = plt.figure()
    ax1 = plt.subplot2grid((n,1), (0,0), rowspan=m)
    ax2 = plt.subplot2grid((n,1), (m,0))
    for cp in self.UNITS:
      ax1.plot(self.T, cp.Vm, linewidth=2.0)
      ax1.set_ylim(vlimit)
      if cp.iClamper:
        ax2.plot(self.T, cp.iClamper.Command, linewidth=1.0)

    ax1.set_ylabel('voltage (mV)')
    ax2.set_xlabel('time (ms)')
    ax2.set_ylabel('current (pA/um2)')
    plt.show()

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
 (9). Specific membrane conductance (Gm):       nS / um2
 (10).Specific membrance capacitance (Cm):      pF / um2
 (11).Specific axial conductance (Ga):          nS / um

Symbols:
  transmembrane current density (Jm): pA/um2
  cytosolic current density (Jn): pA ??
  maximum conductances per unit area (gMax): nS/um2
  axial current density (Jx?):
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
    if N > 1:
      print ('%i compartments supplied\n'%(N))
    else:
      print ('%i compartment supplied\n'%(N))

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

  def run(self, t, dt=0.005):
    steps = int(t/dt)
    self.T = [i*dt+self.Clock for i in range(steps)]
    self.Clock += t

    print ("Integration in process")
    print ('Number of steps: %i'% (steps))
    print ("Step length: %5.4f ms" % (dt))
    print ('Wait!')
    print ('...')

    self.C_Clampers = []
    for cp in self.UNITS: # check for missing values
      for ch in cp.channel_list:
        if ch.Tag == 'LGIC':
          if ch.Ligand.Clamper: self.C_Clampers.append(ch.Ligand.Clamper)

    self.Rates = []
    for cp in self.UNITS:
      for x in cp.RateList:
        self.Rates.append(x)

    ### prepare to run
    self.VC_UNITS = [i for i in self.UNITS if i.vClamper]
    self.CC_UNITS = [i for i in self.UNITS if i.cClamper]
    self.IC_UNITS = [i for i in self.UNITS if i.iClamper]
    self.IM_UNITS = [i for i in self.UNITS if i.iMonitor]
    self.NC_UNITS = list(set(self.UNITS)-set(self.VC_UNITS))

    for cp in self.UNITS: # check for missing values
      cp.dt = dt
      cp.T = self.T
      for ch in cp.channel_list:
        if ch.gMax == None:
          raise Exception('Ion channel %s in Compartment %i doen not have gMax value.')

        if ch.ER == None:
          raise Exception('Ion channel %s in Compartment %i doen not have ER value.')

    for cp in self.UNITS: # prepare the storages based on run parameters
      cp.Vm = array.array('f',[cp.V0]) * steps
      for ch in cp.channel_list:
        ch.gM = array.array('f',[0]) * steps

    for cp in self.VC_UNITS:
      cp.vClamper.T = self.T
      cp.vClamper.Jm = array.array('f',[0]) * steps
      cp.vClamper.Jn = array.array('f',[0]) * steps

    for cp in self.IM_UNITS:
      cp.iMonitor.Jm = array.array('f',[0]) * steps
      cp.iMonitor.Jn = array.array('f',[0]) * steps
      cp.iMonitor.Jc = array.array('f',[0]) * steps

    for cp in self.IC_UNITS: # generate command from specified waveform
      ic = cp.iClamper
      if ic.Waveform:
        ic.Command = [ic.Waveform._func(self.T[i]) for i in range(steps)]

    for cc in self.C_Clampers: # We access the clamper through Experiment
      cc.T = self.T
      cc.Command = [cc.Waveform._func(i*dt) for i in range(steps)]

    for cp in self.VC_UNITS: # generate command from specified waveform
      vc = cp.vClamper
      if vc.Waveform:
        vc.Command = [vc.Waveform._func(i*dt) + vc.Baseline for i in range(steps)]
        vc.dot = [0.] + [ (vc.Command[i]-vc.Command[i-1])/dt for i in range(1,steps) ]

    ###-------------------------------###
    t0 = time.time()
    for step_num in range(steps): # looping begins
      for cp in self.UNITS:
        if cp.iClamper:      # we don't directly access the clamper
          cp.Ja = cp.iClamper.Command[step_num]
        else:
          cp.Ja = 0.0

        if cp.Parent:        # get Jx (axial current density)
          cp.Jx = (cp.Vi - cp.Parent.Vi) * cp.gx/cp.Surface
        else:
          cp.Jx = 0.0

      for cp in self.UNITS:  # sum up non-transmembrane, non-capacitive current
        cp.Jn = cp.Ja - cp.Jx
        for chd in cp.Child:
          cp.Jn += chd.Jx

      for cc in self.C_Clampers: # We access the clamper through Experiment
        cc.Ligand.Conc = cc.Command[step_num]

      for x in self.Rates:
        x.update()

      # Update_Start ========
      for cp in self.NC_UNITS:
        cp._update_Vm(dt)

        for ch in cp.channel_list:
          ch._update_gate(cp.Vi,dt)
          #ch._update_gMax(dt)

      for cp in self.IM_UNITS: # Current Recording
        im = cp.iMonitor
        im.Jm[step_num] = cp.get_Jion(cp.Vi)
        im.Jn[step_num] = cp.Jn
        im.Jc[step_num] = cp.Cm * cp.Vdot

      for cp in self.VC_UNITS: # Voltage-Clamp block
        vc = cp.vClamper
        cp.Vi = vc.Command[step_num]
        for ch in cp.channel_list:
          ch._update_gate(cp.Vi,dt)

        vc.Jm[step_num] = cp.get_Jion(vc.Command[step_num])
        vc.Jn[step_num] = cp.Jn

      # Update_End ========
      for cp in self.UNITS: # Store membrane potentials
        cp.Vm[step_num] = cp.Vi
        #for ch in cp.channel_list:
        #  ch.gM[step_num] = ch.gMax

    t1 = time.time()
    print ('Integration done in %f seconds.  ^-^\n' % (t1-t0))
    ##########################

  """def plot(self,ylim=[-90,60]):
    fig = plt.figure()
    for cp in self.UNITS:
      plt.plot(self.T, cp.Vm, linewidth=2.0)
    plt.ylim(ylim)
    plt.xlabel('time (ms)')
    plt.ylabel('mV')
    plt.show()"""
  def plot(self,ylim=[-90,60]):
    n = 4
    m = n-1
    fig = plt.figure()
    ax1 = plt.subplot2grid((n,1), (0,0), rowspan=m)
    ax2 = plt.subplot2grid((n,1), (m,0))
    for cp in self.UNITS:
      ax1.plot(self.T, cp.Vm, linewidth=2.0)
      ax1.set_ylim(ylim)
      if cp.iClamper:
        ax2.plot(self.T, cp.iClamper.Command, linewidth=1.0)

    ax1.set_ylabel('voltage (mV)')
    ax2.set_xlabel('time (ms)')
    ax2.set_ylabel('current (pA/um2)')
    plt.show()


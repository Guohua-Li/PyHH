__doc__ = """
                      PyHH Neuron Simulator

PyHH is a small Python library for simulation of ion channel dynamics and the
electrical activity of cells. It provides two ways of simulating ion channels:
the Hodgkin-Huxley approach and the Markov approach.

Basic units used in PyHH are:
 (1). time:          ms           (5). capacitance:   pF
 (2). voltage:       mV           (6). length:        um
 (3). current:       pA           (7). area:          um2
 (4). conductance:   nS           (8). concentration: mM

 (9). Specific membrane conductance (Gm):       nS / um2
 (10).Specific membrance capacitance (Cm):      pF / um2
 (11).Specific axial conductance (Ga):          nS / um

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
import matplotlib.pyplot as pl
from math import exp, log10
import copy
import array
import time


def arange(a,b,step):
  N = int((b-a)/step)
  return [a+i*step for i in range(N)]

def isCompartment(x):
  return hasattr(x,'is_Compartment')


def Fm(V): # Na channel, m gate
  v = V + 35
  if v == 0: alpha = 1.0
  else:      alpha = 0.1 * v /(1-exp(-v/10.0))
  beta = 4.0 * exp(-(V+60)/18)
  tau = 1.0/(alpha + beta)
  return tau, alpha * tau

def Fh(V): # Na channel, h gate
  alpha = 0.07 * exp(-0.05*(V+60.0))
  beta  = 1. / (1+exp(-0.1*(V+30.0)))
  tau = 1. / (alpha + beta)
  return tau, alpha * tau

def Fn(V): # K channel, n gate
  if V != -50.0: alpha = 0.01 * (V + 50.0) / (1-exp(-(V + 50.0)/10))
  else:          alpha = 0.1
  beta = 0.125 * exp(-0.0125*(V+60))
  tau = 1./(alpha + beta)
  return tau, alpha * tau



def _plt(f1,f2=None,f3=None):
  Start, Finish = -100, 60
  V = arange(Start, Finish, 0.5)
  X1, tau1, inf1 = [], [], []
  for i in range(len(V)):
    try:
      x, y = f1(V[i])
      X1.append(V[i])
      tau1.append(x)
      inf1.append(y)
    except:
      pass
  if f2:
    X2, tau2, inf2 = [], [], []
    for i in range(len(V)):
      try:
        x, y = f2(V[i])
        X2.append(V[i])
        tau2.append(x)
        inf2.append(y)
      except:
        pass
  if f3:
    X3, tau3, inf3 = [], [], []
    for i in range(len(V)):
      try:
        x, y = f3(V[i])
        X3.append(V[i])
        tau3.append(x)
        inf3.append(y)
      except:
        pass
  fig = pl.figure()
  if f1 and f2 and f3:
    ax1 = fig.add_subplot(321)
    ax2 = fig.add_subplot(322)
    ax3 = fig.add_subplot(323)
    ax4 = fig.add_subplot(324)
    ax5 = fig.add_subplot(325)
    ax6 = fig.add_subplot(326)
    ax5.set_xlabel('mV')
    ax6.set_xlabel('mV')
  elif f1 and f2:
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    ax3.set_xlabel('mV')
    ax4.set_xlabel('mV')
  else:
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.set_xlabel('mV')
    ax2.set_xlabel('mV')

  ax1.plot(X1, tau1, label="line 1")
  ax2.plot(X1, inf1)
  ax1.set_xlim([Start,Finish])
  ax2.set_xlim([Start,Finish])

  if f2:
    ax3.plot(X2, tau2, label="line 2")
    ax4.plot(X2, inf2)
    ax3.set_xlim([Start,Finish])
    ax4.set_xlim([Start,Finish])

  if f3:
    ax5.plot(X3, tau3, label="line 3")
    ax6.plot(X3, inf3)
    ax5.set_xlim([Start,Finish])
    ax6.set_xlim([Start,Finish])

  #legend = ax1.legend(loc='upper right', shadow=True) #, fontsize='x-large'
  pl.show()


COUNT = 0

def increment():
    global COUNT
    COUNT = COUNT+1


class Alpha: 
  """
  alpha function for current injection
  """
  def __init__(self, delay, tau, amplitude):
    self.Delay = delay
    self.Tau = tau
    self.Amplitude = amplitude

  def _func(self, t):
    t = t - self.Delay
    if   t < 0: return 0.0
    else:       return 2.7183*self.Amplitude/self.Tau * t *exp(-t/self.Tau)

class Rect:
  """
  Rectangular function for current, voltage and concentration clampers
  """
  def __init__(self, delay, width, amplitude):
    self.Delay = delay
    self.Width = width
    self.Amplitude = amplitude

  def _func(self, t):
    t = t - self.Delay
    if   t < 0:          return 0.0
    elif t < self.Width: return self.Amplitude
    else:                return 0.0

class Train:
  """
  A train of Rectangular pulses
  """
  def __init__(self, delay = 1, width = 0.4, interval = 1.5, amplitude = 0.5, number = 5):
    self.Delay = delay
    self.Amplitude = amplitude
    self.Width = width
    self.interval = interval
    self.Number = number

  def _func(self, t):
    t = t - self.Delay
    P = (self.Width+self.interval)
    n = t // P
    if n >= self.Number: return 0.0
    if t < 0: return 0.0
    t = t % P
    if t < self.Width: return self.Amplitude
    else: return 0.0


class CClamper:
  """
  Concentration Clamper
  """
  def __init__(self, ligand):
    self.Command  = None
    self.Waveform = None
    self.Ligand   = ligand
    ligand.Clamper  = self    # added after v1.0
    #self.Ligand.Clamper = self
    self.Tag      = 'Ligand'

  """def connect(self, cmpt):
    if type(cmpt) is list:
      for cp in cmpt:
        if cp.cClamper:
          raise Exception('The compartment already has a CClamper')
        cp.cClamper = self
    else:
      if cmpt.cClamper:
        raise Exception('The compartment already has a CClamper')
      cmpt.cClamper = self"""

  def set_amplitude(self, val):
    self.Waveorm.Amplitude = val

  def set_width(self, val):
    self.Waveorm.Width = val



class IClamper:
  """
  Current Clamper
  """
  def __init__(self):
    self.Command  = None
    self.Waveform = None
    self.Tag      = 'Current'

  def connect(self, cmpt):
    if type(cmpt) is list:
      for cp in cmpt:
        if cp.iClamper:
          raise Exception('The compartment already has an IClamper')
        cp.iClamper = self
    else:
      if cmpt.iClamper:
        raise Exception('The compartment already has an IClamper')
      cmpt.iClamper = self

  def set_amplitude(self, val):
    self.Waveorm.Amplitude = val

  def set_width(self, val):
    self.Waveorm.Width = val


class VClamper:
  """
  Voltage clamper
  """
  def __init__(self):
    self.Command  = None
    self.Waveform = None
    self.Baseline = -60
    self.Tag      = 'Voltage'
    self.Jc     = None # to store capacity current
    self.Jn     = None # to store cytosolic current
    self.Jm     = None # to store transmembrane current

  def connect(self, cmpt):
    if type(cmpt) is list:
      for cp in cmpt:
        if cp.vClamper:
          raise Exception('The compartment already has a VClamper')
        cp.vClamper = self
    else:
      if cmpt.vClamper:
        raise Exception('The compartment already has a VClamper')
      cmpt.vClamper = self

  def set_baseline(self, val):
    self.Baseline = val

  def set_amplitude(self, val):
    self.Waveorm.Amplitude = val

  def set_width(self, val):
    self.Waveorm.Width = val

  def calc_Jp(self):
    Jp = [i+j+k for i,j,k in zip(self.Jm, self.Jc, self.Jn)]
    return Jp

  def save(self, filename):
    N = len(self.Jm)
    f = open(filename,'w')
    for k in range(N):
      s = '%7.5f %7.5f %7.5f\n'%(self.Jm[k],self.Jn[k],self.Jc[k])
      f.write(s)
    f.close()

class IMonitor:
  """
  Current monitor
  """
  def __init__(self):
    self.Jc     = None # to store capacity current
    self.Jn     = None # to store cytosolic current
    self.Jm     = None # to store transmembrane current

  def save(self, filename):
    N = len(self.Jm)
    f = open(filename,'w')
    for k in range(N):
      s = '%7.5f %7.5f %7.5f\n'%(self.Jm[k],self.Jn[k],self.Jc[k])
      f.write(s)
    f.close()



class LeakChannel:
  """
  Leak conductance
  """
  def __init__(self, gMax=0.03, ER=-60):
    self.gMax = gMax
    self.ER   = ER
    self.Tag = 'Leak'
    self.ptr = None

  def _ProbOpen(self): return 1.0
  def _update_gate(self, V, dt): return
  def _update_gMax(self,dt): return
  def store_gates(self, mem=0): return
  def load_gates(self, mem=0): return
  def show(self): pass


class KChannel:
  """
  Classical potassium channel
  """
  def __init__(self, f=Fn, gMax=0.18, ER=-90):
    self._f = f
    self.a0 = 0.299
    self.a  = self.a0
    self.gMax = gMax
    self.ER   = ER
    self.Fg   = gMax
    self.Rg   = 0.01
    self.gM   = None  # for storage of the gMax
    self.ptr = None
    self.Tag  = 'KDR'

  def _ProbOpen(self):
    return self.a**4

  def _update_gate(self, V, dt):
    k, f = self._f(V)
    self.a = exp(-dt/k) * (self.a-f) + f

  def _update_gMax(self, dt):
    gf = self.Fg * self.a**4
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

  def store_gates(self, mem=0):
    self.a0 = self.a

  def load_gates(self, mem = 0):
    self.a = self.a0

  def show(self):
    print(self.Tag)
    print('a =  %4.3f'% (self.a))

  def plot(self):
    _plt(self._f)


class NaChannel:
  """
  Classical sodium channel
  """
  def __init__(self, fa, fb, gMax, ER):
    self._fa = fa
    self._fb = fb
    self.a0 = 0.046
    self.b0 = 0.638
    self.a  = self.a0
    self.b  = self.b0
    self.gMax = gMax
    self.Fg   = gMax
    self.ER   = ER
    self.Rg   = 0.01
    self.gM   = None   # for storage of the gMax
    self.Tag  = 'NaC'
    self.ptr = None

  def _ProbOpen(self):
    return (self.a**3) * self.b

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb

  def show(self):
    print(self.Tag)
    print('a =  %4.3f' % (self.a))
    print('b =  %4.3f' % (self.b))

  def _update_gMax(self, dt):
    gf = self.Fg * self.a**3 * self.b
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

  def store_gates(self, mem = 0):
    if mem == 0: self.a0, self.b0 = self.a, self.b
    else: self.a1, self.b1 = self.a, self.b

  def load_gates(self, mem = 0):
    if mem == 0: self.a, self.b = self.a0, self.b0
    else: self.a, self.b = self.a1, self.b1

  def plot(self):
    _plt(self._fa, self._fb)


class Compartment:
  """
  Compartment on which neurons are built.
  """
  is_Compartment = 1
  def __init__(self, diameter = 20, length = 50):
    self.V0    = -61.237
    self.Vi    = self.V0
    self.Vdot  = 0.
    self.Cm    = 0.01
    self.Gx    = 700.0
    self.Length = length
    self.Diameter = diameter
    self.iClamper = None
    self.vClamper = None
    self.cClamper = None
    self.iMonitor = None
    self.Parent = None
    self.Child = []
    self.ID = COUNT
    increment()
    self.Vm    = None # for storage of the membrane potential
    self.channel_list = [] 

  def get_Jion(self, V):
    Jion = 0.0
    for ch in self.channel_list:
      g = ch.gMax * ch._ProbOpen()
      Jion += g * (V - ch.ER)
    return Jion

  def add_channel(self, Channel):
    if type(Channel) is list:
      raise Exception('add_channel does not accept a list')

    if Channel.Tag == 'LGIC':
      Channel.binding = copy.deepcopy(Channel.binding)
      ch = LGIC(Channel.transit,Channel.binding,Channel.gMax,Channel.ER) # python3.5
    else:
      ch = copy.deepcopy(Channel)
    self.channel_list.append(ch)
    ch.ptr = self
    return ch

  def add_iclamper(self):
    clamper = IClamper()
    self.iClamper = clamper
    return clamper

  def add_vclamper(self):
    clamper = VClamper()
    self.vClamper = clamper
    return clamper

  """def add_cclamper(self, lig):
    clamper = CClamper(lig)
    self.cClamper = clamper
    return clamper"""

  def add_imonitor(self):
    monitor = IMonitor()
    self.iMonitor = monitor
    return monitor

  def add_channels(self, Channels):
    if type(Channels) is not list:
      raise Exception('add_channels accept a list')

    for ch in Channels:
      self.add_channel(ch)

  def _calc_surface(self):
    if self.Length == None:
      self.Surface = 3.1415926 * (self.Diameter**2) # sphere
      return

    self.Surface = 3.1415926 * self.Diameter * self.Length
    S = 3.1415926 * self.Diameter * self.Diameter / 4.
    if self.Parent:
      if self.Child == []:
        self.Surface = self.Surface + S
    else:
      self.Surface = self.Surface + S
      if self.Child == []:
        self.Surface = self.Surface + S

  def _calc_gx(self):
    if self.Length:
      Cross = 3.1415926 * self.Diameter * self.Diameter / 4.
      self.gx = self.Gx * Cross / self.Length
    else:
      self.gx = 0

  def connect(self, cmpt):
    if not cmpt.Length:
      raise Exception('Do not connect a cylinder to a sphere.')
    self.Child.append(cmpt)
    cmpt.Parent = self

  def attached_to(self, cmpt):
    if not self.Length:
      raise Exception('Do not attach a sphere to a cylinder.')
    cmpt.Child.append(self)
    self.Parent = cmpt

  def _update_Vm(self, dt):
    P = self.Vi
    if self.channel_list == []:
      self.Vi = P + dt * self.Jn / self.Cm
      return
    sigma_g  = sigma_gE = 0.0
    for ch in self.channel_list:
      g = ch.gMax * ch._ProbOpen()
      sigma_g += g
      sigma_gE += g * ch.ER
    """if sigma_g == 0:
      self.Vi = P + dt * self.Jn / self.Cm
      return"""
    k = sigma_g / self.Cm
    Vf = (self.Jn + sigma_gE) / sigma_g
    self.Vi = exp(-k*dt) * (P-Vf) + Vf
    self.Vdot = (self.Vi-P)/dt

  def save(self, filename):
    N = len(self.Vm)
    f = open(filename,'w')
    for k in range(N):
      f.write('%8.6f\n' %(self.Vm[k]))
    f.close()

  def set_V(self, val):
    self.V0    = val
    self.Vi    = val

  def show(self):
    print ('Vm =%f' % (self.Vi))


class Experiment:
  """
  Integration
  """
  def __init__(self, prep):
    if type(prep) is list:
      self.UNITS = []
      for a in prep:
        if isCompartment(a):
          self.UNITS
        else:
          print ("-------------------------------------")
          print ("| Check the input to the experiment |")
          print ("-------------------------------------")
          raise Exception("Input error")
      self.UNITS = prep
    else:
      if isCompartment(prep):
        self.UNITS = [prep]
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

  def run(self, t, dt):
    steps = int(t/dt)
    self.T = [i*dt+self.Clock for i in range(steps)]
    self.Clock += t

    print ("Integration in process")
    print ('Number of steps: %i'% (steps))
    print ("Step length: %5.4f ms" % (dt))
    print ('Wait!')
    print ('...')

    self.C_Clampers = [cp.cClamper for cp in self.UNITS if cp.cClamper]
    self.C_Clampers = []
    for cp in self.UNITS: # check for missing values
      for ch in cp.channel_list:
        if ch.Tag == 'LGIC':
          if ch.Ligand.Clamper: self.C_Clampers.append(ch.Ligand.Clamper)

    ### prepare to run
    self.VC_UNITS = [i for i in self.UNITS if i.vClamper]
    self.CC_UNITS = [i for i in self.UNITS if i.cClamper]
    self.IC_UNITS = [i for i in self.UNITS if i.iClamper]
    self.IM_UNITS = [i for i in self.UNITS if i.iMonitor]
    self.NC_UNITS = list(set(self.UNITS)-set(self.VC_UNITS))

    for cp in self.UNITS: # check for missing values
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
      cp.vClamper.Jm = array.array('f',[0]) * steps
      cp.vClamper.Jn = array.array('f',[0]) * steps
      cp.vClamper.Jc = array.array('f',[0]) * steps

    for cp in self.IM_UNITS:
      cp.iMonitor.Jm = array.array('f',[0]) * steps
      cp.iMonitor.Jn = array.array('f',[0]) * steps
      cp.iMonitor.Jc = array.array('f',[0]) * steps

    for cp in self.IC_UNITS: # generate command from specified waveform
      ic = cp.iClamper
      if ic.Waveform:
        ic.Command = [ic.Waveform._func(self.T[i]) for i in range(steps)]

    #for cp in self.CC_UNITS: # generate command from specified waveform
    #  cc = cp.cClamper
    for cc in self.C_Clampers: # We access the clamper through Experiment
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

      # Update_Start ========
      for cp in self.NC_UNITS:
        cp._update_Vm(dt)

        for ch in cp.channel_list:
          ch._update_gate(cp.Vi,dt)
          #ch._update_gMax(dt)

      for cp in self.IM_UNITS: # Current Recording
        im = cp.iMonitor
        im.Jc[step_num] = cp.Cm * cp.Vdot
        im.Jm[step_num] = cp.get_Jion(cp.Vi)
        im.Jn[step_num] = cp.Jn

      for cp in self.VC_UNITS: # Voltage-Clamp block
        vc = cp.vClamper
        cp.Vi = vc.Command[step_num]
        for ch in cp.channel_list:
          ch._update_gate(cp.Vi,dt)

        vc.Jc[step_num] = cp.Cm * vc.dot[step_num]
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



class Ligand:
  """
  This class defines ligand objects.
  """
  def __init__(self, concentration = 0):
    self.Conc = concentration
    self.Clamper = None

class LGIC:
  """
  This class defines LGICs according to Markov transition schemes.
  """
  def __init__(self, transit, binding, gMax, ER):
    self.transit = transit # keep it for copying the channel
    self.binding = binding # keep it for copying the channel
    self.State = 0
    self.do_parsing(transit, binding)
    self.N = len(self.RateMatrix)
    self.StateProb = [0] * self.N    # create and initialize the state probability vector
    self.StateProb[0] = 1.0
    self.gMax = gMax
    self.ER   = ER
    self.Ligands = list(binding.keys())
    if len(self.Ligands) == 1: self.Ligand = self.Ligands[0]
    self.Tag = 'LGIC'

  def do_parsing(self, transit, binding):
    States = list(transit.keys()) # 2.7 and 3.5 are different
    States.sort()
    N = len(States)
    StateIndexDict = {}
    for k in range(N):
      s = States[k]
      StateIndexDict[s] = k

    self.RateMatrix = [[1]*N for i in range(N)] # like ones((N, N), float)
    for i in range(N):
      si = States[i]
      for j in range(N):
        sj = States[j]
        if sj in transit[si].keys():
          self.RateMatrix[i][j] = transit[si][sj]
        else:
          self.RateMatrix[i][j] = 0.0

    self.GateDict = {}
    for lg in binding.keys():
      paths = []
      scheme = binding[lg]
      for a in scheme.keys():
        b = scheme[a]
        i = StateIndexDict[a]
        j = StateIndexDict[b]
        paths.append((i,j))
      self.GateDict[lg] = paths

    self.OpenIndex = []
    for k in range(N):
      if States[k][0] == 'O': self.OpenIndex.append(k)

  def _ProbOpen(self):
    gate = 0.0
    for k in self.OpenIndex: gate += self.StateProb[k]
    return gate

  def _update_gate(self, V, h):
    N = self.N
    rate = [[j for j in row] for row in self.RateMatrix] # copy self.RateMatrix
    for lg, paths in self.GateDict.items():
      for s, t in paths:
        rate[s][t]   = rate[s][t] * lg.Conc
    A = [sum(i) for i in rate]           # sum of each row
    m = [ [j*self.StateProb[i] for j in rate[i]] for i in range(len(rate)) ] # rate * self.XC
    """
    m = []
    for i in range(len(rate)):
      row = rate[i]
      k = self.StateProb[i]
      n = [j*k for j in row]
      m.append(n)"""
    B = list(map(sum,zip(*m))) # sum of each vol, faster
    for k in range(N):
      if A[k] == 0:
        self.StateProb[k] = self.StateProb[k] + h * B[k]
      else:
        C = B[k]/A[k]
        self.StateProb[k] = (self.StateProb[k] - C) * exp(-A[k]*h) + C
    S = sum(self.StateProb)
    if S != 0:
      self.StateProb = [x/S for x in self.StateProb]

  def show(self):
    N = len(self.StateProb)
    print ("State Probabilities:")
    for k in range(N):
      print ("%5.4f"%(self.StateProb[k]))


class CaPool:
  """
  Calcium microdomain
  """
  def __init__(self, func):
    self._f = func
    self.Conc = 0.0001
    self.Target = None

  def update(self, J, dt):
    tau, inf = self._f(J)
    self.Conc = exp(-dt/tau) * (self.Conc-inf) + inf


NaC = NaChannel(Fm, Fh, gMax=0.6, ER=50) # the unit for gMax is nS/um2
KDR = KChannel(Fn, gMax=0.18, ER=-90)
gL  = LeakChannel(gMax=0.03,ER=-60)



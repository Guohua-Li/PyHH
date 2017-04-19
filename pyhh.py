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

def isLGIC(x):
  return hasattr(x,'is_LGIC')


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
  alpha function for Waveform setting
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
    ligand.Clamper  = self
    self.Tag      = 'Ligand'

  def set_amplitude(self, val):
    self.Waveform.Amplitude = val

  def set_width(self, val):
    self.Waveform.Width = val

  def plot(self):
    fig = pl.figure()
    pl.plot(self.T, self.Command, linewidth=2.0)
    pl.ylim()
    pl.xlabel('time (ms)')
    pl.ylabel('mV')
    pl.show()


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
    self.Waveform.Amplitude = val

  def set_width(self, val):
    self.Waveform.Width = val


class VClamper:
  """
  Voltage clamper
  """
  def __init__(self):
    self.Command  = None
    self.Waveform = None
    self.Baseline = -60
    self.Tag      = 'Voltage'
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
    self.Waveform.Amplitude = val

  def set_width(self, val):
    self.Waveform.Width = val

  def calc_Jp(self):
    Jp = [i+j for i,j in zip(self.Jm, self.Jn)]
    return Jp

  def save(self, filename):
    N = len(self.Jm)
    f = open(filename,'w')
    for k in range(N):
      s = '%7.5f %7.5f\n'%(self.Jm[k],self.Jn[k])
      f.write(s)
    f.close()

  """def plot(self, show_Jm=1,show_Jn=1,show_Jp=0):
    fig = pl.figure()
    pl.subplot(2,1,1)
    if show_Jm==1: pl.plot(self.T, self.Jm, linewidth=3.0)
    if show_Jn==1: pl.plot(self.T, self.Jn, linewidth=2.0)
    if show_Jp==1:
      Jp = self.calc_Jp()
      pl.plot(self.T, Jp, linewidth=1.0)

    a = max(self.Jm)
    b = max(self.Jn)
    c = min(self.Jm)
    d = min(self.Jn)
    M = max(a,b)
    m = min(c,d)
    R = (M-m)/10.
    pl.ylim([m-R,M+R])
    pl.ylabel('Jm, Jn')
    pl.subplot(2,1,2)
    pl.plot(self.T, self.Command, linewidth=1.0)
    pl.xlabel('time (ms)')
    pl.ylim([self.Baseline-5,self.Waveform.Amplitude+self.Baseline+5])
    pl.show()"""

  def plot(self, show_Jm=1,show_Jn=1,show_Jp=0):
    fig = pl.figure()
    ax1 = pl.subplot2grid((3,1), (0,0), rowspan=2)
    ax2 = pl.subplot2grid((3,1), (2,0))

    if show_Jm==1: ax1.plot(self.T, self.Jm, linewidth=3.0)
    if show_Jn==1: ax1.plot(self.T, self.Jn, linewidth=2.0)
    if show_Jp==1:
      Jp = self.calc_Jp()
      ax1.plot(self.T, Jp, linewidth=1.0)

    a = max(self.Jm)
    b = max(self.Jn)
    c = min(self.Jm)
    d = min(self.Jn)
    M = max(a,b)
    m = min(c,d)
    R = (M-m)/10.
    pl.ylim([m-R,M+R])
    pl.ylabel('Jm, Jn')

    ax2.plot(self.T, self.Command, linewidth=1.0)
    ax2.set_xlabel('time (ms)')
    ax2.set_ylim([self.Baseline-5,self.Waveform.Amplitude+self.Baseline+5])
    pl.show()


class IMonitor:
  """
  Current monitor
  """
  def __init__(self):
    self.Jn     = None # to store cytosolic current
    self.Jm     = None # to store transmembrane current
    self.Jc     = None # to store capacitive current

  def save(self, filename):
    N = len(self.Jm)
    f = open(filename,'w')
    for k in range(N):
      s = '%7.5f %7.5f %7.5f\n'%(self.Jm[k], self.Jn[k], self.Jc[k])
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
  def _set_gMax(self, val): self.gMax = val
  def _set_ER(self, val): self.gMax = val



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

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

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

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

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
  def __init__(self, diameter = 20, length = 50, channel_list=None):
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
    if channel_list != None:
      self.add_channels(channel_list)
    else:
      self.channel_list = []

  def get_Jion(self, V):
    Jion = 0.0
    for ch in self.channel_list:
      g = ch.gMax * ch._ProbOpen()
      Jion += g * (V - ch.ER)
    return Jion

  def add_channels(self,*channels):
    L = []
    for i in channels:
      if type (i) is list:
        for b in i:
          ch = self._add_channel(b)
          L.append(ch)
      else:
        ch = self._add_channel(i)
        L.append(ch)
    if len(L)==1: return L[0]
    return tuple(L)

  def _add_channel(self, Channel):
    if isLGIC(Channel):
      Channel.binding = copy.deepcopy(Channel.binding)
      ch = LGIC(Channel.transit,Channel.binding,Channel.gMax,Channel.ER) # python3.5
    else:
      ch = copy.deepcopy(Channel)
    self.channel_list.append(ch)
    ch.ptr = self
    return ch

  def add_vgic(self, gMax, ER, a = None, b = None, c = None):
    if a == None and b == None and c == None:
      ch = LeakChannel(gMax, ER)
      self.channel_list.append(ch)
      ch.ptr = self
      return ch

    if b == None and c == None:
      Fa, A = a[0], a[1]
      if A == 1:
        ch = HH_a(Fa, gMax, ER)#
        self.channel_list.append(ch)
        ch.ptr = self
        return ch

      if A > 1:
        if A == 4: # special case
          ch = KChannel(Fa, gMax, ER)
          self.channel_list.append(ch)
          ch.ptr = self
        else:
          ch = HH_aA(Fa, A, gMax, ER)#
          self.channel_list.append(ch)
          ch.ptr = self
      return ch

    if c == None:
      Fa, A = a[0], a[1]
      Fb, B = b[0], b[1]
      if B == 1 and A == 1:
        ch = HH_ab(Fa, Fb, gMax, ER)#
        self.channel_list.append(ch)
        ch.ptr = self
        return ch

      if B == 1 and A >  1:
        if A == 3: # special case
          ch = NaChannel(Fa, Fb, gMax, ER)
        else:
          ch = HH_aAb(Fa, A, Fb, gMax, ER)

        self.channel_list.append(ch)
        ch.ptr = self
        return ch

      if B >  1 and A >  1:
        ch = HH_aAbB(Fa, A, Fb, B, gMax, ER)
        self.channel_list.append(ch)
        ch.ptr = self
        return ch

    Fa, A = a[0], a[1]
    Fb, B = b[0], b[1]
    Fc, C = c[0], c[1]
    if C == 1 and B == 1 and A == 1: ch = HH_abc(Fa,Fb,Fc, gMax, ER)
    if C == 1 and B == 1 and A >  1:
      if A == 3:
        ch = NaChannel_MHMJ(Fa, Fb, Fc, gMax, ER)
      else:
        ch = HH_aAbc(Fa, A, Fb, Fc, gMax, ER)
    if C == 1 and B  > 1 and A >  1: ch = HH_aAbBc(Fa,A, Fb, B, Fc, gMax, ER)
    if C >  1 and B  > 1 and A >  1: ch = HH_aAbBcC(Fa, A, Fb, B, Fc, C, gMax, ER)
    self.channel_list.append(ch)
    ch.ptr = self
    return ch

  def add_lgic(self, transit, binding, gMax, ER):
    t = copy.deepcopy(transit)
    b = copy.deepcopy(binding)
    ch = LGIC(t, b, gMax, ER)
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

  def add_imonitor(self):
    monitor = IMonitor()
    self.iMonitor = monitor
    return monitor

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

  def plot(self):
    fig = pl.figure()
    if self.iClamper:
      ax1 = pl.subplot2grid((3,1), (0,0), rowspan=2)
      ax2 = pl.subplot2grid((3,1), (2,0))
      ax1.plot(self.T, self.Vm, linewidth=1.0)
      ax2.plot(self.T, self.iClamper.Command, linewidth=1.0)
      ax2.set_xlabel('time (ms)')
      #ax2.set_ylim([self.Baseline-5,self.Waveform.Amplitude+self.Baseline+5])
    else:
      N = len(self.Vm)
      if N == 0: return
      pl.plot(self.T, self.Vm, linewidth=2.0)
      #pl.ylim(ylim)
      pl.xlabel('time (ms)')
      pl.ylabel('mV')
    pl.show()




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

  def run(self, t, dt):
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
    fig = pl.figure()
    for cp in self.UNITS:
      pl.plot(self.T, cp.Vm, linewidth=2.0)
    pl.ylim(ylim)
    pl.xlabel('time (ms)')
    pl.ylabel('mV')
    pl.show()"""
  def plot(self,ylim=[-90,60]):
    fig = pl.figure()
    ax1 = pl.subplot2grid((3,1), (0,0), rowspan=2)
    ax2 = pl.subplot2grid((3,1), (2,0))
    for cp in self.UNITS:
      ax1.plot(self.T, cp.Vm, linewidth=2.0)
      if cp.iClamper:
        ax2.plot(self.T, cp.iClamper.Command, linewidth=1.0)

    ax1.set_ylabel('voltage (mV)')
    ax2.set_xlabel('time (ms)')
    ax2.set_ylabel('current (pA/um2)')
    pl.show()



class Ligand:
  """
  This class defines ligand objects.
  """
  def __init__(self, concentration = 0):
    self.Conc = concentration
    self.Clamper = None

class LGIC:
  """
  LGICs based on Markov transition schemes.
  """
  is_LGIC = 1
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

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

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




###################################################################################
class HH_a:
  """
  gating: a
  """
  def __init__(self, f, gMax, ER):
    self._f = f
    self.a0 = 0.046
    self.a  = self.a0
    self.gMax = gMax
    self.ER = ER
    self.Pools = []
    self.Tag = 'HH_a'
    self.ptr = None

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d
  def _ProbOpen(self):
    return self.a

  def _update_gate(self, V, dt):
    tau, inf = self._f(V)
    self.a = exp(-dt/tau) * (self.a-inf) + inf
    if len(self.Pools) > 0:     #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def set_a(self,a):
    self.a = a

  def show(self):
    print(self.Tag)
    print('a', self.a)

  def plot(self):
    _plt(self._f)


class HH_aAb:
  """
  gating: (a**A) * b
  """
  def __init__(self, fa, A, fb, gMax, ER):
    self._fa = fa
    self._fb = fb
    self.A = A
    self.a0 = 0.046
    self.b0 = 0.638
    self.a  = self.a0
    self.b  = self.b0
    self.gMax = gMax
    self.ER = ER
    self.Pools = []
    self.Tag = 'HH_a(A)b'
    self.ptr = None

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

  def couple(self, p):
    self.pool_pointer = p
    p.target = self

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return (self.a**self.A) * self.b

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    if len(self.Pools) > 0:
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)

  def plot(self):
    _plt(self._fa, self._fb)


class HH_aAbB:
  """
  gating: (a**N) * (b**M)
  """
  def __init__(self, fa, A, fb, B, gMax, ER):
    self._fa = fa
    self._fb = fb
    self.A = A
    self.B = B
    self.a0 = 0.046
    self.b0 = 0.638
    self.a  = 0.046
    self.b  = 0.638
    self.gMax = gMax
    self.ER = ER
    self.Pools = []
    self.Tag = 'a(A)b(B)'
    self.ptr = None

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return (self.a**self.A) * (self.b**self.B)

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb

    if len(self.Pools) > 0:
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)

  def plot(self):
    _plt(self._fa, self._fb)


class HH_aA:
  """
  gating: a**A
  """
  def __init__(self, f, A, gMax, ER):
    self._Func = f
    self.N = A
    self.a0 = 0.299
    self.a = self.a0
    self.gMax = gMax
    self.ER = ER
    self.Pools = []
    self.Tag = 'HH_a(A)'
    self.ptr = None

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return self.a**self.N

  def _update_gate(self, V, dt):
    k, f = self._Func(V)
    self.a = exp(-dt/k) * (self.a-f) + f

    if len(self.Pools) > 0:    #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def show(self):
    print(self.Tag)
    print('a', self.a)

  def plot(self):
    _plt(self._fa)


class HH_ab:
  """
  gating: a*b
  """
  def __init__(self, fa, fb, gMax, ER):
    self._fa = fa
    self._fb= fb
    self.a0 = 0.046
    self.b0 = 0.046
    self.a  = 0.046
    self.b  = 0.046
    self.gMax = gMax
    self.ER = ER
    self.Pools = []
    self.Tag = 'HH_ab'
    self.ptr = None

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return self.a * self.b

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb

    if len(self.Pools) > 0:    #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def set_a(self,a):
    self.a = a

  def set_b(self,b):
    self.b = b

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)

  def plot(self):
    _plt(self._fa, self._fb)


class HH_abc:
  """
  gating: a*b*c
  """
  def __init__(self, fa, fb, fc, gMax, ER):
    self._fa = fa
    self._fb = fb
    self._fc = fc
    self.a0 = 0.046
    self.b0 = 0.638
    self.c0 = 0.638
    self.a  = 0.046
    self.b  = 0.638
    self.c  = 0.638
    self.gMax = gMax
    self.Fg   = gMax
    self.ER   = ER
    self.Rg   = 0.01
    self.gM   = None
    self.Pools = []
    self.Tag  = 'HH_abc'
    self.ptr = None

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return self.a * self.b * self.c

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    tc, fc = self._fc(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    self.c = exp(-dt/tc) * (self.c-fc) + fc

    if len(self.Pools) > 0:    #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def _update_gMax(self, dt):
    gf = self.Fg * self.a**3 * self.b * self.c
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

  def set_m(self, a):
    self.a = a

  def set_h(self, b):
    self.b = b

  def set_i(self, c):
    self.c = c

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)
    print('c', self.c)

  def plot(self):
    _plt(self._fa, self._fb, self._fc)


class HH_aAbc:
  """
  gating: (a**A) * b * c
  """
  def __init__(self, fa, A, fb, fc, gMax, ER):
    self._fa = fa
    self._fb = fb
    self._fc = fc
    self.A = A
    self.a0 = 0.046
    self.b0 = 0.638
    self.c0 = 0.638
    self.a  = 0.046
    self.b  = 0.638
    self.c  = 0.638
    self.gMax = gMax
    self.Fg   = gMax
    self.ER   = ER
    self.Rg   = 0.01
    self.gM   = None
    self.Pools = []
    self.Tag  = 'HH_aAbc'
    self.ptr = None

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return (self.a**self.A) * self.b * self.c

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    tc, fc = self._fc(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    self.c = exp(-dt/tc) * (self.c-fc) + fc

    if len(self.Pools) > 0:    #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def _update_gMax(self, dt):
    gf = self.Fg * self.a**3 * self.b * self.c
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

  def set_a(self, a):
    self.a = a

  def set_b(self, b):
    self.b = b

  def set_c(self, c):
    self.c = c

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)
    print('c', self.c)

  def plot(self):
    _plt(self._fa, self._fb, self._fc)


class HH_aAbBc:
  """
  gating: (a**A) * (b**B) * c
  """
  def __init__(self, fa, A, fb, B, fc, gMax, ER):
    self._fa = fa
    self._fb = fb
    self._fc = fc
    self.A = A
    self.B = B
    self.a0 = 0.046
    self.b0 = 0.638
    self.c0 = 0.638
    self.a  = self.a0
    self.b  = self.b0
    self.c  = self.c0
    self.gMax = gMax
    self.Fg   = gMax
    self.ER   = ER
    self.Rg   = 0.01
    self.gM   = None
    self.Pools = []
    self.Tag  = 'HH_aAbBc'
    self.ptr = None

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return (self.a**self.A) * (self.b**self.B) * self.c

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    tc, fc = self._fc(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    self.c = exp(-dt/tc) * (self.c-fc) + fc

    if len(self.Pools) > 0:    #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def _update_gMax(self, dt):
    gf = self.Fg * self.a**3 * self.b * self.c
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

  def set_a(self, a):
    self.a = a

  def set_b(self, b):
    self.b = b

  def set_c(self, c):
    self.c = c

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)
    print('c', self.c)

  def plot(self):
    _plt(self._fa, self._fb, self._fc)


class HH_aAbBcC:
  """
  gating: (a**A) * (b**B) * (c**C)
  """
  def __init__(self, fa, A, fb, B, fc, C, gMax, ER):
    self._fa = fa
    self._fb = fb
    self._fc = fc
    self.A = A
    self.B = B
    self.B = C
    self.a0 = 0.046
    self.b0 = 0.638
    self.c0 = 0.638
    self.a  = self.a0
    self.b  = self.b0
    self.c  = self.c0
    self.gMax = gMax
    self.Fg   = gMax
    self.ER   = ER
    self.Rg   = 0.01
    self.gM   = None
    self.Pools = []
    self.Tag  = 'HH_a(A)b(B)c(C)'
    self.ptr = None

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return (self.a**self.A) * (self.b**self.B) * (self.c**self.C)

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    tc, fc = self._fc(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    self.c = exp(-dt/tc) * (self.c-fc) + fc

    if len(self.Pools) > 0:    #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def _update_gMax(self, dt):
    gf = self.Fg * self.a**3 * self.b * self.c
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

  def set_a(self, a):
    self.a = a

  def set_b(self, b):
    self.b = b

  def set_c(self, c):
    self.c = c

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)
    print('c', self.c)

  def plot(self):
    _plt(self._fa, self._fb, self._fc)


class NaChannel_MHMJ:
  """
  MHMJ sodium channel
  """
  def __init__(self, fm, fh, fi, gMax, ER):
    self._fa = fm
    self._fb = fh
    self._fc = fi
    self.a  = 0.046
    self.b  = 0.638
    self.c  = 0.638
    self.a0 = 0.046
    self.b0 = 0.638
    self.c0 = 0.638
    self.gMax = gMax
    self.Fg   = gMax
    self.ER   = ER
    self.Rg   = 0.01
    self.gM   = None
    self.Tag  = 'MHMJ'

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return (self.a**3) * self.b * self.c

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    tc, fc = self._fc(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    self.c = exp(-dt/tc) * (self.c-fc) + fc

  def _update_gMax(self, dt):
    gf = self.Fg * self.a**3 * self.b * self.c
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

  def set_m(self, m):
    self.a = m
  def set_h(self, h):
    self.b = h
  def set_i(self, i):
    self.c = i

  def show(self):
    print(self.Tag)
    print('m', self.a)
    print('h', self.b)
    print('i', self.c)

  def plot(self):
    _plt(self._fa, self._fb, self._fc)


def voltage_gated(gMax, ER, a = None, b = None, c = None):
  if a == None and b == None and c == None:
    return LeakChannel(gMax, ER)

  if b == None and c == None:
    Fa, A = a[0], a[1]
    if A == 1: return HH_a(Fa, gMax, ER)#
    if A >  1:
      if A == 4: # special case
        return KChannel(Fa, gMax, ER)
      else:
        return HH_aA(Fa, A, gMax, ER)#

  if c == None:
    Fa, A = a[0], a[1]
    Fb, B = b[0], b[1]
    if B == 1 and A == 1: return HH_ab(Fa, Fb, gMax, ER)#
    if B == 1 and A >  1:
      if A == 3: # special case
        return NaChannel(Fa, Fb, gMax, ER)
      else:
        return HH_aAb(Fa, A, Fb, gMax, ER)
    if B >  1 and A >  1: return HH_aAbB(Fa, A, Fb, B, gMax, ER)

  Fa, A = a[0], a[1]
  Fb, B = b[0], b[1]
  Fc, C = c[0], c[1]
  if C == 1 and B == 1 and A == 1: return HH_abc(Fa,Fb,Fc, gMax, ER)
  if C == 1 and B == 1 and A >  1: return HH_aAbc(Fa,A,Fb,Fc, gMax, ER)
  if C == 1 and B  > 1 and A >  1: return HH_aAbBc(Fa,A, Fb, B, Fc, gMax, ER)
  if C >  1 and B  > 1 and A >  1: return HH_aAbBcC(Fa, A, Fb, B, Fc, C, gMax, ER)

###############################################################################
class HHC_a:
  """
  gating: a
  """
  def __init__(self, f, gMax, ER):
    self._f = f
    self.a  = 0.046
    self.a0 = 0.046
    self.gMax = gMax
    self.ER = ER
    self.Pools = []
    self.Tag = 'HHC_a'
    self.ptr = None

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

  def couple(self,p):
    self.pool_pointer = p
    p.Target = self

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return self.a

  def _update_gate(self, V, dt):
    if self.pool_pointer.Conc == 0: C = 0.0001
    else: C = self.pool_pointer.Conc
    tau, inf = self._f(V, C)
    self.a = exp(-dt/tau) * (self.a-inf) + inf
    if len(self.Pools) > 0:
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def set_a(self,a):
    self.a = a

  def show(self):
    print(self.Tag)
    print('a', self.a)

  def plot(self):
    _plt(self._f)


class HHC_aA:
  """
  gating: (a**A)
  """
  def __init__(self, fa, A, gMax, ER):
    self._fa = fa
    self.A = A
    self.a0 = 0.046
    self.a  = 0.046
    self.gMax = gMax
    self.ER = ER
    self.Pools = []
    self.Tag = 'HHC_a(A)'
    self.ptr = None

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

  def couple(self, p):
    self.pool_pointer = p
    p.target = self

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return (self.a**self.A)

  def _update_gate(self, V, dt):
    if self.pool_pointer.Conc == 0: C = 0.0001
    else: C = self.pool_pointer.Conc
    ta, fa = self._fa(V,C)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    if len(self.Pools) > 0:
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def plot(self):
    _plt(self._fa)


class HHC_aAb:
  """
  gating: (a**A) * b
  """
  def __init__(self, fa, A, fb, gMax, ER):
    self._fa = fa
    self._fb = fb
    self.A = A
    self.a0 = 0.046
    self.b0 = 0.638
    self.a  = 0.046
    self.b  = 0.638
    self.gMax = gMax
    self.ER = ER
    self.Pools = []
    self.Tag = 'HHC_a(A)b'
    self.ptr = None

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

  def couple(self, p):
    self.pool_pointer = p
    p.target = self

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return (self.a**self.A) * self.b

  def _update_gate(self, V, dt):
    if self.pool_pointer.Conc == 0: C = 0.0001
    else: C = self.pool_pointer.Conc
    ta, fa = self._fa(V,C)
    tb, fb = self._fb(V,V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    if len(self.Pools) > 0:
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)

  def plot(self):
    _plt(self._fa, self._fb)


class HHC_ab:
  """
  gating: a * b
  """
  def __init__(self, fa, fb, gMax, ER):
    self._fa = fa
    self._fb= fb
    self.a0 = 0.046
    self.b0 = 0.046
    self.a  = 0.046
    self.b  = 0.046
    self.gMax = gMax
    self.ER = ER
    self.Pools = []
    self.Tag = 'HHC_ab'
    self.ptr = None

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

  def couple(self,p):
    self.pool_pointer = p
    p.Target = self

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return self.a * self.b

  def _update_gate(self, V, dt):
    if self.pool_pointer.Conc == 0: C = 0.0001
    else: C = self.pool_pointer.Conc
    ta, fa = self._fa(V,C)
    tb, fb = self._fb(V,C)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb

    if len(self.Pools) > 0:    #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def set_a(self,a):
    self.a = a

  def set_b(self,b):
    self.b = b

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)

  def plot(self):
    _plt(self._fa, self._fb)

class HHC_aAbB:
  """
  (a**N) * (b**M)
  """
  def __init__(self, fa, A, fb, B, gMax, ER):
    self._fa = fa
    self._fb = fb
    self.A = A
    self.B = B
    self.a0 = 0.046
    self.b0 = 0.638
    self.a  = 0.046
    self.b  = 0.638
    self.gMax = gMax
    self.ER = ER
    self.Pools = []
    self.Tag = 'a(A)b(B)'
    self.ptr = None

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

  def couple(self, p):
    self.pool_pointer = p
    p.target = self

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return (self.a**self.A) * (self.b**self.B)

  def _update_gate(self, V, dt):
    if self.pool_pointer.Conc == 0: C = 0.0001
    else: C = self.pool_pointer.Conc
    ta, fa = self._fa(V,C)
    tb, fb = self._fb(V,C)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb

    if len(self.Pools) > 0:    #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)

  def plot(self):
    _plt(self._fa, self._fb)



class HHC_abc:
  """
  gating: a * b * c
  """
  def __init__(self, fa, fb, fc, gMax, ER):
    self._fa = fa
    self._fb = fb
    self._fc = fc
    self.a0 = 0.046
    self.b0 = 0.638
    self.c0 = 0.638
    self.a  = 0.046
    self.b  = 0.638
    self.c  = 0.638
    self.gMax = gMax
    self.Fg   = gMax
    self.ER   = ER
    self.Rg   = 0.01
    self.gM   = None
    self.Pools = []
    self.Tag  = 'HHC_abc'
    self.ptr = None

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

  def couple(self, p):
    self.pool_pointer = p
    p.target = self

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return self.a * self.b * self.c

  def _update_gate(self, V, dt):
    if self.pool_pointer.Conc == 0: C = 0.0001
    else: C = self.pool_pointer.Conc
    ta, fa = self._fa(V,C)
    tb, fb = self._fb(V,C)
    tc, fc = self._fc(V,C)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    self.c = exp(-dt/tc) * (self.c-fc) + fc

    if len(self.Pools) > 0:    #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def _update_gMax(self, dt):
    gf = self.Fg * self.a**3 * self.b * self.c
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

  def set_m(self, a):
    self.a = a

  def set_h(self, b):
    self.b = b

  def set_i(self, c):
    self.c = c

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)
    print('c', self.c)

  def plot(self):
    _plt(self._fa, self._fb, self._fc)

class HHC_aAbc:
  """
  gating: (a**A) * b * c
  """
  def __init__(self, fa, A, fb, fc, gMax, ER):
    self._fa = fa
    self._fb = fb
    self._fc = fc
    self.A = A
    self.a0 = 0.046
    self.b0 = 0.638
    self.c0 = 0.638
    self.a  = 0.046
    self.b  = 0.638
    self.c  = 0.638
    self.gMax = gMax
    self.Fg   = gMax
    self.ER   = ER
    self.Rg   = 0.01
    self.gM   = None
    self.Pools = []
    self.Tag  = 'HHC_aAbc'
    self.ptr = None

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

  def couple(self, p):
    self.pool_pointer = p
    p.target = self

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return (self.a**self.A) * self.b * self.c

  def _update_gate(self, V, dt):
    if self.pool_pointer.Conc == 0: C = 0.0001
    else: C = self.pool_pointer.Conc
    ta, fa = self._fa(V,C)
    tb, fb = self._fb(V,C)
    tc, fc = self._fc(V,C)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    self.c = exp(-dt/tc) * (self.c-fc) + fc

    if len(self.Pools) > 0:    #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def _update_gMax(self, dt):
    gf = self.Fg * self.a**3 * self.b * self.c
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

  def set_a(self, a):
    self.a = a

  def set_b(self, b):
    self.b = b

  def set_c(self, c):
    self.c = c

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)
    print('c', self.c)

  def plot(self):
    _plt(self._fa, self._fb, self._fc)




class HHC_aAbBc:
  """
  gating: (a**A) * (b**B) * c
  """
  def __init__(self, fa, A, fb, B, fc, gMax, ER):
    self._fa = fa
    self._fb = fb
    self._fc = fc
    self.A = A
    self.B = B
    self.a0 = 0.046
    self.b0 = 0.638
    self.c0 = 0.638
    self.a  = self.a0
    self.b  = self.b0
    self.c  = self.c0
    self.gMax = gMax
    self.Fg   = gMax
    self.ER   = ER
    self.Rg   = 0.01
    self.gM   = None
    self.Pools = []
    self.Tag  = 'HH_a(A)b(B)c'
    self.ptr = None

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

  def couple(self, p):
    self.pool_pointer = p
    p.target = self

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return (self.a**self.A) * (self.b**self.B) * self.c

  def _update_gate(self, V, dt):
    if self.pool_pointer.Conc == 0: C = 0.0001
    else: C = self.pool_pointer.Conc
    ta, fa = self._fa(V,C)
    tb, fb = self._fb(V,C)
    tc, fc = self._fc(V,C)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    self.c = exp(-dt/tc) * (self.c-fc) + fc

    if len(self.Pools) > 0:    #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def _update_gMax(self, dt):
    gf = self.Fg * self.a**3 * self.b * self.c
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

  def set_a(self, a):
    self.a = a

  def set_b(self, b):
    self.b = b

  def set_c(self, c):
    self.c = c

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)
    print('c', self.c)

  def plot(self):
    _plt(self._fa, self._fb, self._fc)


class HHC_aAbBcC:
  """
  (a**A) * (b**B) * (c**C)
  """
  def __init__(self, fa, A, fb, B, fc, C, gMax, ER):
    self._fa = fa
    self._fb = fb
    self._fc = fc
    self.A = A
    self.B = B
    self.B = C
    self.a0 = 0.046
    self.b0 = 0.638
    self.c0 = 0.638
    self.a  = self.a0
    self.b  = self.b0
    self.c  = self.c0
    self.gMax = gMax
    self.Fg   = gMax
    self.ER   = ER
    self.Rg   = 0.01
    self.gM   = None
    self.Pools = []
    self.Tag  = 'HH_a(A)b(B)c(C)'
    self.ptr = None

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

  def couple(self, p):
    self.pool_pointer = p
    p.target = self

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return (self.a**self.A) * (self.b**self.B) * (self.c**self.C)

  def _update_gate(self, V, dt):
    if self.pool_pointer.Conc == 0: C = 0.0001
    else: C = self.pool_pointer.Conc
    ta, fa = self._fa(V,C)
    tb, fb = self._fb(V,C)
    tc, fc = self._fc(V,C)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    self.c = exp(-dt/tc) * (self.c-fc) + fc

    if len(self.Pools) > 0:    #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def _update_gMax(self, dt):
    gf = self.Fg * self.a**3 * self.b * self.c
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

  def set_a(self, a):
    self.a = a

  def set_b(self, b):
    self.b = b

  def set_c(self, c):
    self.c = c

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)
    print('c', self.c)

  def plot(self):
    _plt(self._fa, self._fb, self._fc)

def calcium_activated(gMax, ER, a, b = None, c = None):
  if b == None and c == None:
    Fa, A = a[0], a[1]
    if A == 1: return HHC_a(Fa, gMax, ER)#
    if A >  1: return HHC_aA(Fa, A, gMax, ER)#

  if c == None:
    Fa, A = a[0], a[1]
    Fb, B = b[0], b[1]
    if B == 1 and A == 1: return HHC_ab(Fa, Fb, gMax, ER)#
    if B == 1 and A >  1: return HHC_aAb(Fa, A, Fb, gMax, ER)
    if B >  1 and A >  1: return HHC_aAbB(Fa, A, Fb, B, gMax, ER)

  Fa, A = a[0], a[1]
  Fb, B = b[0], b[1]
  Fc, C = c[0], c[1]
  
  if C == 1 and B == 1 and A == 1: return HHC_abc(Fa,Fb,Fc, gMax, ER)
  if C == 1 and B == 1 and A >  1: return HHC_aAbc(Fa,A,Fb,Fc, gMax, ER)
  if C == 1 and B  > 1 and A >  1: return HHC_aAbBc(Fa,A,Fb,B,Fc, gMax, ER)
  if C >  1 and B  > 1 and A >  1: return HHC_aAbBcC(Fa,A,Fb,B,Fc,C, gMax, ER)


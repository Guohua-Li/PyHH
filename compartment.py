import copy
from lgic import LGIC, TransRate #
from clampers import VClamper, IClamper#
from math import exp#
import matplotlib.pyplot as plt #

COUNT = 0

def increment():
    global COUNT
    COUNT = COUNT+1

def isLGIC(x):
  return hasattr(x,'is_LGIC')

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
    self.RateList = []
    self.ID = COUNT
    increment()
    self.Vm    = None # for storage of the membrane potential
    self.channel_list = []
    if channel_list != None:
      self.add_channels(channel_list)

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
      Channel.transit = copy.deepcopy(Channel.transit)
      Channel.binding = copy.deepcopy(Channel.binding)
      ch = LGIC(Channel.transit, Channel.binding, Channel.gMax, Channel.ER) # python3.5
      """for i in ch.Transit:
        if isinstance(i, TransRate):
          i.Ptr = self
          self.RateList.append(i)"""
      for i in range(ch.nStates):
        for j in range(ch.nStates):
          x = ch.Transit[i][j]
          if isinstance(x, TransRate):
            x.Compt = self
            self.RateList.append(x)

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

  def add_vclamper(self, baseline):
    clamper = VClamper(baseline)
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
    S = 3.1415926 * self.Diameter * self.Diameter / 4.0
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
      self.gx = self.Gx * Cross / self.Length # gx is the axial conductance (nS)
    else:
      self.gx = 0

  def connect(self, child_cmpt): # branch
    if not child_cmpt.Length:
      raise Exception('Do not connect a cylinder to a sphere.')
    self.Child.append(child_cmpt)
    child_cmpt.Parent = self

  def attached_to(self, parent_cmpt):
    if not self.Length:
      raise Exception('Do not attach a sphere to a cylinder.')
    parent_cmpt.Child.append(self)
    self.Parent = parent_cmpt

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
    M,m = max(self.Vm), min(self.Vm)
    R = (M - m)/10
    fig = plt.figure()
    if self.iClamper:
      ax1 = plt.subplot2grid((4,1), (0,0), rowspan=3)
      ax2 = plt.subplot2grid((4,1), (3,0))
      ax1.plot(self.T, self.Vm, linewidth=1.0)
      ax2.plot(self.T, self.iClamper.Command, linewidth=1.0)
      ax1.set_ylim([m-R,M+R])
      ax2.set_xlabel('time (ms)')
      ax2.set_ylim([0,1.1*self.iClamper.Waveform.Amplitude])
    else:
      N = len(self.Vm)
      if N == 0: return
      plt.plot(self.T, self.Vm, linewidth=2.0)
      plt.xlabel('time (ms)')
      plt.ylabel('mV')
      plt.set_ylim([m-R,M+R])
    plt.show()



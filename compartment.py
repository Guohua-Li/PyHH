import copy

from lgic import LGIC, Transition
from clampers import VClamper, IClamper
from math import exp
import matplotlib.pyplot as plt

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
  def __init__(self, diameter = 20, length = 50, channels=None, parent=None):
    self.V0    = -61.237
    self.Vi    = self.V0 # Vi: intracellular potential of the compartment, just a single value
    self.Vm    = None # an array.array to store membrane, modified in the Experiment class
    self.Vdot  = 0.0
    self.Cm    = 0.01 # Specific membrance capacitance:      pF / um2
    self.Gx    = 700.0
    # cytoplasm conductivity: the reciprocal of electrical resistivity
    # or speciﬁc electrical resistance for the cytoplasm modelled as a resistor
    # The new value is taken from 
    # Wang, K., Zhao, Y., Chen, D. et al. Specific membrane capacitance, cytoplasm conductivity and
    # instantaneous Young’s modulus of single tumour cells. Sci Data 4, 170015
    # (2017). https://doi.org/10.1038/sdata.2017.15
    #self.Gx    = 958.0 # average of values from 5 types of cell, the unit is ns/um
    self.Length = length
    self.Diameter = diameter
    self.iClamper = None
    self.vClamper = None
    self.cClamper = None
    self.iMonitor = None

    self.Parent = None
    if parent != None:
      self.attached_to(parent)

    self.Children = []
    self.RateList = []
    self.ID = COUNT
    increment()

    # newly define the make the program clear
    # J_axial is not used in the current class
    self.J_axial = 0.0 # axial current density, direction is from the child to parent, see Experiment class

    # J_total is used in the current class for updating Vi
    self.J_total = 0.0 # all current injected into the segment from neighboring compartments and IClamper

    self.Density = {} # testing
    self.Reversal = {} #testing

    self.channel_list = []
    if channels != None:
      for ch_obj in channels.keys():
        copied = self._add_channel_copy(ch_obj) # copy and add
        self.Density[copied]  = channels[ch_obj][0]
        self.Reversal[copied] = channels[ch_obj][1]

  def add_channels(self, channels): # testing, seem ok
    L = []
    for ch_obj in channels.keys():
      copied = self._add_channel_copy(ch_obj) # copy and add
      self.Density[copied]  = channels[ch_obj][0]
      self.Reversal[copied] = channels[ch_obj][1]
      L.append(copied)

    if len(L)==1: return L[0]
    return tuple(L)

  def _add_channel_copy(self, Channel):
    if isLGIC(Channel):
      Channel.transit = copy.deepcopy(Channel.transit)
      Channel.binding = copy.deepcopy(Channel.binding)
      ch_obj = LGIC(Channel.transit, Channel.binding)
      """for i in ch_obj.Transit:
        if isinstance(i, Transition):
          i.Ptr = self
          self.RateList.append(i)"""

      state_numbers = range(ch_obj.nStates)
      for state_i in state_numbers:
        for state_j in state_numbers:
          x = ch_obj.Transit[state_i][state_j] # x is a number or a Transition object
          if isinstance(x, Transition):
            x.Compt = self
            self.RateList.append(x) # only the voltageGate
    else: # non-lgic
      ch_obj = copy.deepcopy(Channel)

    self.channel_list.append(ch_obj)
    ch_obj.ptr = self
    return ch_obj

  def add_vgic(self, a = None, b = None, c = None):
    if a == None and b == None and c == None:
      ch_obj = LeakChannel()
      self.channel_list.append(ch_obj)
      ch_obj.ptr = self
      return ch_obj

    if b == None and c == None:
      Fa, A = a[0], a[1]
      if A == 1:
        ch = HH_a(Fa)#
        self.channel_list.append(ch)
        ch.ptr = self
        return ch

      if A > 1:
        if A == 4: # special case
          ch = KChannel(Fa)
          self.channel_list.append(ch)
          ch.ptr = self
        else:
          ch = HH_aA(Fa, A)#
          self.channel_list.append(ch)
          ch.ptr = self
      return ch

    if c == None:
      Fa, A = a[0], a[1]
      Fb, B = b[0], b[1]
      if B == 1 and A == 1:
        ch = HH_ab(Fa, Fb)#
        self.channel_list.append(ch)
        ch.ptr = self
        return ch

      if B == 1 and A >  1:
        if A == 3: # special case
          ch = NaChannel(Fa, Fb)
        else:
          ch = HH_aAb(Fa, A, Fb)

        self.channel_list.append(ch)
        ch.ptr = self
        return ch

      if B >  1 and A >  1:
        ch = HH_aAbB(Fa, A, Fb, B)
        self.channel_list.append(ch)
        ch.ptr = self
        return ch

    Fa, A = a[0], a[1]
    Fb, B = b[0], b[1]
    Fc, C = c[0], c[1]
    if C == 1 and B == 1 and A == 1: ch = HH_abc(Fa,Fb,Fc)
    if C == 1 and B == 1 and A >  1:
      if A == 3:
        ch = NaChannel_MHMJ(Fa, Fb, Fc)
      else:
        ch = HH_aAbc(Fa, A, Fb, Fc)
    if C == 1 and B  > 1 and A >  1: ch = HH_aAbBc(Fa,A, Fb, B, Fc)
    if C >  1 and B  > 1 and A >  1: ch = HH_aAbBcC(Fa, A, Fb, B, Fc, C)
    self.channel_list.append(ch)
    ch.ptr = self
    return ch

  def add_lgic(self, transit, binding, gMax, ER):
    t = copy.deepcopy(transit)
    b = copy.deepcopy(binding)
    ch = LGIC(t, b)
    self.Density[ch] = gMax
    self.Reversal[ch] = ER
    self.channel_list.append(ch)
    ch.ptr = self
    return ch

  def get_Jion(self, V):
    # this is required by voltage clamp,
    # not current clamp
    Jion = 0.0
    for ch_obj in self.channel_list:
      g = self.Density[ch_obj] * ch_obj._ProbOpen() # ch_obj.gMax
      Jion += g * (V - self.Reversal[ch_obj]) # ch_obj.ER
    return Jion

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
    S = 3.1415926 * self.Diameter * self.Diameter / 4.0 # this is the circle area
    if self.Parent:
      if self.Children == []:
        self.Surface = self.Surface + S # terminal segment?
    else:
      self.Surface = self.Surface + S # one end is closed
      if self.Children == []:
        self.Surface = self.Surface + S # the other end is closed

  def _calc_gx(self):
    if self.Length:
      Cross = 3.1415926 * self.Diameter * self.Diameter / 4.0
      self.gx = self.Gx * Cross / self.Length # gx is the axial conductance (nS)
    else:
      self.gx = 0

  def connect(self, child_cmpt): # branch
    if not child_cmpt.Length:
      raise Exception('Do not connect a cylinder to a sphere.')
    self.Children.append(child_cmpt)
    child_cmpt.Parent = self

  def attached_to(self, parent_cmpt):
    #if not self.Length:
    if self.Length == None:
      raise Exception('Do not attach a sphere to a cylinder.')
    parent_cmpt.Children.append(self)
    self.Parent = parent_cmpt

  def _update_Vm(self, dt):
    if self.channel_list == []: # passive circuit
      self.Vi = self.Vi + dt * self.J_total / self.Cm # Specific membrance capacitance: pF/um2
      return

    sigma_g  = sigma_gE = 0.0
    for ch_obj in self.channel_list:
      g = self.Density[ch_obj] * ch_obj._ProbOpen() # g: unit area conductance, ch_obj.gMax
      sigma_g += g
      sigma_gE += g * self.Reversal[ch_obj] #ch_obj.ER

    V0 = self.Vi
    """if sigma_g == 0:
      self.Vi = V0 + dt * (self.J_total+sigma_gE) / self.Cm
      return"""

    k = sigma_g / self.Cm # Specific membrance capacitance: pF/um2
    V_inf = (self.J_total + sigma_gE) / sigma_g
    self.Vi = exp(-k*dt) * (V0-V_inf) + V_inf
    self.Vdot = (self.Vi-V0)/dt

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

  def plot_Vm(self, vlimit=None):
    M, m = max(self.Vm), min(self.Vm)
    if vlimit == None and M > m:
      R = (M - m)/10
      vlimit = [m-R,M+R]

    fig = plt.figure()
    if self.iClamper: # we need to show the command
      ax1 = plt.subplot2grid((4,1), (0,0), rowspan=3)
      ax2 = plt.subplot2grid((4,1), (3,0))
      ax1.plot(self.T, self.Vm, linewidth=1.0)
      ax2.plot(self.T, self.iClamper.Command, linewidth=1.0)
      if vlimit != None:
        ax1.set_ylim(vlimit)
      ax2.set_xlabel('time (ms)')
      ax2.set_ylim([0,1.1*self.iClamper.Waveform.Amplitude])
    else:
      N = len(self.Vm)
      if N == 0: return
      plt.plot(self.T, self.Vm, linewidth=2.0)
      plt.xlabel('time (ms)')
      plt.ylabel('mV')
      if vlimit != None:
        plt.ylim(vlimit)
    plt.show()

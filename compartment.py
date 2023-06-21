# modified on May 19, 2023

import copy
import array
import math

from lgic import LGIC, VariableRate
from clampers import VClamper, IClamper
from math import exp
import matplotlib.pyplot as plt

COUNT = 0

def increment():
  global COUNT
  COUNT += 1


def isLGIC(x):
  return hasattr(x,'is_LGIC')

class Compartment:
  """
  Compartment on which neurons are built.
  """
  is_Compartment = 1
  def __init__(self, diameter = 20, length = None, channels=None, parent=None):
    self.V0    = -61.237
    self.Vi    = self.V0 # Vi: intracellular potential of the compartment, just a single value
    self.Vm    = None # an array.array to store membrane, modified in the Experiment class
    self.Vdot  = 0.0
    self.Cm    = 0.01 # Specific membrance capacitance: pF/um2
    self.Gx    = 700.0

    self.Length = length
    self.Diameter = diameter
    self.iClamper = None
    self.vClamper = None
    self.cClamper = None
    self.iRecorder = None

    self.Parent = None
    if parent != None:
      self.attached_to(parent)

    self.Children = []
    self.VoltGatedRates = [] #RateList
    self.ID = COUNT
    increment()

    # J_axial is not used in the current class
    self.J_axial = 0.0 # axial current density, direction is from the child to parent, see Experiment class

    # J_total is used in the current class for updating Vi
    self.J_total = 0.0 # all current injected into the segment from neighboring compartments and IClamper

    self.Density = {}
    self.Reversal = {}

    self.channel_list = []
    if channels != None:
      for ch_obj in channels.keys():
        copied = self._add_channel_copy(ch_obj) # copy and add to channel_list
        self.Density[copied]  = channels[ch_obj][0]
        self.Reversal[copied] = channels[ch_obj][1]

  def add_channels(self, channels): # testing, seem ok
    L = []
    for ch_obj in channels.keys():
      copied = self._add_channel_copy(ch_obj) # copy and add to channel_list
      self.Density[copied]  = channels[ch_obj][0]
      self.Reversal[copied] = channels[ch_obj][1]
      L.append(copied)

    if len(L)==1: return L[0]
    return tuple(L)

  def get_channel(self, tag):
    for chnn in self.channel_list:
      if chnn.Tag == tag:
        return chnn

  def _add_channel_copy(self, Channel):
    if isLGIC(Channel):
      Channel.transit = copy.deepcopy(Channel.transit)
      Channel.binding = copy.deepcopy(Channel.binding)
      ch_obj = LGIC(Channel.transit, Channel.binding)
      """for i in ch_obj.Transit:
        if isinstance(i, VariableRate):
          i.Ptr = self
          self.VoltGatedRates.append(i) # RateList
          """
      state_numbers = range(ch_obj.nStates)
      for state_i in state_numbers:
        for state_j in state_numbers:
          x = ch_obj.Transit[state_i][state_j] # x is a number or a VariableRate object
          if isinstance(x, VariableRate):
            x.Compt = self
            self.VoltGatedRates.append(x) # RateList
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

  def get_Jion(self, V, step):
    # this is required only by voltage clamp
    Jion = 0.0
    for ch_obj in self.channel_list:
      p = ch_obj._ProbOpen(step)
      g = self.Density[ch_obj] * p # ch_obj.gMax
      Jion += g * (V - self.Reversal[ch_obj]) # ch_obj.ER
    return Jion

  def add_iclamper(self, waveform=None, delay=None, amplitude=None, width=None, tau=None, rlist=None, interval=None, number=None): #
    clamper = IClamper()
    self.iClamper = clamper
    if waveform != None:
      clamper.set_waveform(waveform, delay, amplitude, width, tau, rlist, interval, number)
    return clamper

  def add_vclamper(self, baseline, waveform=None, delay=None, amplitude=None, width=None, tau=None, rlist=[]):
    clamper = VClamper(self, baseline)
    if waveform != None: clamper.set_waveform(waveform, delay, amplitude, width, tau, rlist=[])
    self.vClamper = clamper
    return clamper

  def _calc_surface(self):
    if self.Length == None: # sphere
      self.Surface = 3.1415926 * (self.Diameter**2) # sphere area
      for child in self.Children:
        self.Surface -= 3.1415926 * child.Diameter * child.Diameter / 4.0
    else:
      self.Surface = 3.1415926 * self.Diameter * self.Length # curved surface area
      base_area = 3.1415926 * self.Diameter**2 / 4.0
      if self.Parent:
        if self.Children == []:
          self.Surface += base_area # This is a terminal segment.
      else:
        self.Surface += base_area # one end is closed
        if self.Children == []:
          self.Surface += base_area # the other end is closed too

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

  def _update_Vm(self, dt, step):
    if self.channel_list == []: # passive circuit
      self.Vi = self.Vi + dt * self.J_total / self.Cm # Specific membrance capacitance: pF/um2
      return

    sigma_g  = sigma_gE = 0.0
    for ch_obj in self.channel_list:
      p = ch_obj._ProbOpen(step)
      g = self.Density[ch_obj] * p # g: unit area conductance, ch_obj.gMax
      sigma_g += g
      sigma_gE += g * self.Reversal[ch_obj] #ch_obj.ER
      """
      if ch_obj.record_mode=="none": continue
      if   ch_obj.record_mode == "g":
        ch_obj.trace[step] = g
        ptint("g=",g)
      elif ch_obj.record_mode == "p":
        ch_obj.trace[step] = p
      """
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

  def _plot_Vm(self, xlim, ylim, show_density=True):
    N = len(self.Vm)
    if N == 0: return

    fig = plt.figure()
    """plt.plot(self.T, self.Vm, linewidth=1.0)
    plt.xlabel('ms')
    plt.ylabel('mV')
    plt.xlim(xlim)
    plt.ylim(ylim)"""

    ax = fig.add_subplot(111)
    ax.tick_params(axis='y', labelcolor='tab:blue')
    ax.plot(self.T, self.Vm, linewidth=1.0)
    ax.set_xlabel('ms')
    ax.set_ylabel('mV', color='tab:blue')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.tick_params(axis='y', colors='tab:blue')
    ax.spines['left'].set_color('tab:blue')

    #ax.spines['bottom'].set_color('red')
    #ax.spines['top'].set_color('red')
    #ax.xaxis.label.set_color('red')
    #ax.tick_params(axis='x', colors='red')

    plt.show()

  def _plot_Vm_clamper(self, xlim, ylim, show_density=True):
    fig = plt.figure()
    if show_density:
      i_cmd = self.iClamper.Command
    else:
      N = len(self.iClamper.Command)
      i_cmd = array.array('f', [0]) * N
      for i in range(N):
        i_cmd[i] = self.iClamper.Command[i] * self.Surface/1000.

    ax1 = plt.subplot2grid((4,1), (0,0), rowspan=3)
    ax1.plot(self.T, self.Vm, linewidth=1.0)
    ax1.set_ylim(ylim)
    ax1.set_xlim(xlim)
    ax1.set_ylabel('mV')
    ax1.yaxis.label.set_color('tab:blue')
    #ax1.tick_params(axis='y', labelcolor='tab:blue')
    ax1.tick_params(axis='y', colors='tab:blue')
    ax1.spines['left'].set_color('tab:blue')
    ax1.spines['right'].set_color('tab:gray')
    ax1.spines['top'].set_color('tab:gray')
    ax1.spines['bottom'].set_color('tab:gray')

    ax2 = plt.subplot2grid((4,1), (3,0))
    ax2.plot(self.T, i_cmd, linewidth=1.0)
    ax2.set_xlim(xlim)
    ax2.set_xlabel('ms')
    ax2.yaxis.label.set_color('tab:blue')
    #ax2.tick_params(axis='y', labelcolor='tab:blue')
    ax2.tick_params(axis='y', colors='tab:blue')
    ax2.spines['left'].set_color('tab:blue')
    ax2.spines['right'].set_color('tab:gray')
    ax2.spines['top'].set_color('tab:gray')
    ax2.spines['bottom'].set_color('tab:gray')

    if show_density:
      amp = self.iClamper.Amplitude
      ylabel = r'pA/$\mu$$m^2$'
    else:
      amp = self.iClamper.Amplitude*self.Surface/1000.
      ylabel = "nA"

    if amp < 0:
      y0, ym = 1.1*amp, -0.1*amp
    elif amp < 0.6:
      y0 = -0.05*amp
      ym = 1.1*amp
    else:
      y0 = -0.1*amp
      ym = 1.1*math.ceil(amp)
    ax2.set_ylim([y0, ym])
    ax2.set_ylabel(ylabel)
    plt.subplots_adjust(hspace=0.5)
    plt.show()

  def plot_Vm(self, ylim=None, show_density=True):
    dt = self.T[-1]/(len(self.T)-1)
    xlim = (self.T[0],self.T[-1]+dt)
    M, m = max(self.Vm), min(self.Vm)
    if ylim == None and M > m:
      R = (M - m)/10
      ylim = [m-R,M+R]

    if self.iClamper:
      self._plot_Vm_clamper(xlim, ylim, show_density)
    else:
      self._plot_Vm(xlim, ylim, show_density)

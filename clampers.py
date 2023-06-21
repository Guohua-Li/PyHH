from math import exp
from array import array
import matplotlib.pyplot as plt

class Clamper:

  def __init__(self):
    self.Baseline = 0.0
    self.Waveform  = None
    self.Command  = None
    self.Amplitude = None
    self.Delay = None
    self.Width = None
    self.Tau   = None
    self.rect_list = None
    self.Interval = None

  def set_amplitude(self, val):
    self.Amplitude = val

  def set_width(self, val):
    self.Width = val

  def _rect(self, N, dt):
    self.Command = array('f', [0]) * N
    for i in range(N):
      t = i * dt
      t -= self.Delay
      if t < 0:
        self.Command[i] = self.Baseline
      elif t < self.Width:
        self.Command[i] = self.Baseline + self.Amplitude
      else:
        self.Command[i] = self.Baseline

  def _alpha(self, N, dt):
    self.Command = array('f', [0]) * N
    for i in range(N):
      t = i * dt
      t = t - self.Delay
      if   t < 0: self.Command[i] = 0.0
      else:       self.Command[i] = 2.7183*self.Amplitude/self.Tau * t *exp(-t/self.Tau)

  def _rects(self, N, dt):
    self.Command = array('f',[0]) * N
    T0 = 0
    for st in self.rlist:
        int_st = tuple([round(x/dt) for x in st])
        t1 = int_st[0]
        t2 = t1 + int_st[1]
        for i in range(T0, T0+t1):
            self.Command[i] = 0
        for i in range(T0+t1, T0+t2):
            self.Command[i] = st[2]
        T0 += t2
        if T0 >= steps: break

  def _train(self, N, dt):
    self.Command = array('f',[0]) * N
    for i in range(N):
      t = i * dt - self.Delay
      if t < 0:
        self.Command[i] = 0.0
        continue

      P = (self.Width+self.Interval)
      if (t // P) >= self.Number:
        self.Command[i] = 0.0
        continue

      if t % P < self.Width:
        self.Command[i] = self.Amplitude
      else:
        self.Command[i] =  0.0

  def _get_command(self, N, dt):
    if   self.Waveform == 'rect':
      self._rect(N, dt)
    elif self.Waveform == 'alpha':
      self._alpha(N, dt)
    elif self.Waveform == 'rects':
      self._rects(N, dt)
    elif self.Waveform == 'train':
      self._train(N, dt)
    else:
      self.Command = array('f',[0]) * N
      print("Waveform not defined")
      
  def set_waveform(self, waveform, delay, amplitude=0.0, width=None, tau=None, rlist=None, interval=None, number=None):
    self.Waveform = waveform
    self.Delay = delay
    self.Amplitude = amplitude
    self.Width = width
    self.Tau = tau
    self.rlist = rlist
    self.Interval = interval
    self.Number = number



class CClamper(Clamper):
  """
  Concentration Clamper
  """
  def __init__(self, ligand):
    Clamper.__init__(self)
    self.Ligand    = ligand
    ligand.Clamper = self
    self.Tag       = 'Ligand'

  def plot(self):
    fig = plt.figure()
    plt.plot(self.T, self.Command, linewidth=2.0)
    plt.ylim()
    plt.xlabel('time (ms)')
    plt.ylabel('mV')
    plt.show()


class IClamper(Clamper):
  """
  Current Clamper
  """
  def __init__(self, compartment=None):
    Clamper.__init__(self)
    self.Tag      = 'Current'
    if compartment != None: self.clamp(compartment)

  def clamp(self, cmpt):
    if type(cmpt) is list:
      for cp in cmpt:
        if cp.iClamper:
          raise Exception('The compartment already has an IClamper')
        cp.iClamper = self
    else:
      if cmpt.iClamper:
        raise Exception('The compartment already has an IClamper')
      cmpt.iClamper = self

  def plot(self, T):
    fig = plt.figure()
    plt.plot(T, self.Command, linewidth=2.0)
    plt.ylim()
    plt.xlabel('time (ms)')
    plt.ylabel('mV')
    plt.show()



class VClamper(Clamper):
  """
  Voltage clamper
  """
  def __init__(self, compartment=None, baseline=-60.0):
    #if not isinstance(compartment, Compartment): raise Exception('A compartment instance needed')
    #if !isinstance(baseline
    Clamper.__init__(self)
    self.Baseline = baseline
    self.Tag      = 'Voltage'
    self.J_injs   = None # to store non-transmembrane currents (injected plus axial)
    self.J_ion    = None # to store total transmembrane ionic currents
    if compartment != None: self.clamp(compartment)

  def clamp(self, cmpt):
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

  def save(self, filename):
    N = len(self.J_ion)
    f = open(filename,'w')
    for k in range(N):
      s = '%7.5f %7.5f\n'%(self.J_ion[k],self.J_injs[k])
      f.write(s)
    f.close()

  def plot_current(self, show_J_ion=True, show_J_injs=True, show_Jp=False):
    dt = self.T[-1]/(len(self.T)-1)
    xlim = (self.T[0],self.T[-1]+dt)

    n = 4
    m = n-1
    fig = plt.figure()
    ax1 = plt.subplot2grid((n,1), (0,0), rowspan=m)
    ax2 = plt.subplot2grid((n,1), (m,0))

    if show_J_ion:
      ax1.plot(self.T, self.J_ion, linewidth=1.0)

    if show_J_injs:
      ax1.plot(self.T, self.J_injs, linewidth=1.0)

    if show_Jp:
      #Jp = self.calc_Jp()
      Jp = [j_ion+j_injs for j_ion,j_injs in zip(self.J_ion, self.J_injs)]
      ax1.plot(self.T, Jp, linewidth=1.0)

    ax1.set_ylabel(r'pA/$\mu$$m^2$', color='tab:blue')
    ax1.set_xlim(xlim)
    ax1.set_xticklabels([])
    #ax1.tick_params(axis='x', colors='tab:blue')
    ax1.tick_params(axis='y', colors='tab:blue')
    ax1.spines['left'].set_color('tab:blue')
    ax1.spines['right'].set_color('tab:gray')
    ax1.spines['top'].set_color('tab:gray')
    ax1.spines['bottom'].set_color('tab:gray')

    a = max(self.J_ion)
    b = max(self.J_injs)
    c = min(self.J_ion)
    d = min(self.J_injs)
    M = max(a,b)
    m = min(c,d)
    R = (M-m)/10.
    ax1.set_ylim([m-R,M+R])

    ax2.set_xlim(xlim)
    ax2.plot(self.T, self.Command, linewidth=1.0)
    ax2.set_xlabel('time (ms)')
    y1 = self.Baseline + 5
    if self.Amplitude != None:
      y1 += self.Amplitude
    ax2.set_ylim([self.Baseline-5, y1])
    #ax2.tick_params(axis='x', colors='tab:blue')
    ax2.tick_params(axis='y', colors='tab:blue')
    ax2.spines['left'].set_color('tab:blue')
    ax2.spines['right'].set_color('tab:gray')
    ax2.spines['top'].set_color('tab:gray')
    ax2.spines['bottom'].set_color('tab:gray')

    plt.subplots_adjust(hspace=0.2)
    plt.show()

"""
class IRecorder:
  def __init__(self):
    self.J_ion  = None # to store transmembrane ionic current
    self.J_cap  = None # to store capacitive current density
    self.J_injs = None # none transmembrane current, namely, other net current injected into this compartment

  def save(self, filename):
    N = len(self.J_ion)
    f = open(filename,'w')
    for k in range(N):
      s = '%7.5f %7.5f %7.5f\n'%(self.J_ion[k], self.J_injs[k], self.J_cap[k])
      f.write(s)
    f.close()
"""

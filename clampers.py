from math import exp
import matplotlib.pyplot as plt


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
    if (t // P) >= self.Number: return 0.0
    if t < 0: return 0.0
    t = t % P
    if t < self.Width: return self.Amplitude
    else: return 0.0

class Clamper:

  def __init__(self):
    self.Command  = None
    self.Waveform = None

  def set_amplitude(self, val):
    self.Waveform.Amplitude = val

  def set_width(self, val):
    self.Waveform.Width = val


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
  def __init__(self):
    Clamper.__init__(self)
    self.Tag      = 'Current'

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


class VClamper(Clamper):
  """
  Voltage clamper
  """
  def __init__(self, baseline):
    Clamper.__init__(self)
    self.Baseline = baseline
    self.Tag      = 'Voltage'
    self.J_injs   = None # to store cytosolic current? see Experiment class for meaning
    self.J_ion    = None # to store transmembrane current? see Experiment class for meaning

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

  def calc_Jp(self):
    Jp = [i+j for i,j in zip(self.J_ion, self.J_injs)]
    return Jp

  def save(self, filename):
    N = len(self.J_ion)
    f = open(filename,'w')
    for k in range(N):
      s = '%7.5f %7.5f\n'%(self.J_ion[k],self.J_injs[k])
      f.write(s)
    f.close()

  def plot(self, show_J_ion=True, show_J_injs=True, show_Jp=False):
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
      Jp = self.calc_Jp()
      ax1.plot(self.T, Jp, linewidth=1.0)

    a = max(self.J_ion)
    b = max(self.J_injs)
    c = min(self.J_ion)
    d = min(self.J_injs)
    M = max(a,b)
    m = min(c,d)
    R = (M-m)/10.
    plt.ylim([m-R,M+R])
    plt.ylabel('J_ion, J_injs')

    ax2.plot(self.T, self.Command, linewidth=1.0)
    ax2.set_xlabel('time (ms)')
    ax2.set_ylim([self.Baseline-5,self.Waveform.Amplitude+self.Baseline+5])
    plt.show()


class IMonitor:
  """
  Current monitor
  """
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

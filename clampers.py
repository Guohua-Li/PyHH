from math import exp#
import matplotlib.pyplot as plt#

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
    fig = plt.figure()
    plt.plot(self.T, self.Command, linewidth=2.0)
    plt.ylim()
    plt.xlabel('time (ms)')
    plt.ylabel('mV')
    plt.show()


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
  def __init__(self, baseline):
    self.Command  = None
    self.Waveform = None
    self.Baseline = baseline
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
    n = 4
    m = n-1
    fig = plt.figure()
    ax1 = plt.subplot2grid((n,1), (0,0), rowspan=m)
    ax2 = plt.subplot2grid((n,1), (m,0))

    if show_Jm==1:
      ax1.plot(self.T, self.Jm, linewidth=1.0)

    if show_Jn==1:
      ax1.plot(self.T, self.Jn, linewidth=1.0)

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
    plt.ylim([m-R,M+R])
    plt.ylabel('Jm, Jn')

    ax2.plot(self.T, self.Command, linewidth=1.0)
    ax2.set_xlabel('time (ms)')
    ax2.set_ylim([self.Baseline-5,self.Waveform.Amplitude+self.Baseline+5])
    plt.show()


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


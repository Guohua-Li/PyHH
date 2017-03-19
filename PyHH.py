import matplotlib.pyplot as pl
from math import exp
import copy
import array
import time
import sys

def arange(a,b,step):
  N = int((b-a)/step)
  return [a+i*step for i in range(N)]


def Fm(V): # Func for Na channel
  if V != -35.0: alpha = 0.1 * (V + 35.0) /(1-exp(-(V + 35.0)/10.0))
  else:          alpha = 1.0
  beta = 4 * exp(-(V+60)/18)
  tau = 1./(alpha + beta)
  return tau, alpha * tau

def Fh(V): # Func for Na channel
  alpha = 0.07 * exp(-0.05*(V+60.0))
  beta  = 1. / (1+exp(-0.1*(V+30.0)))
  tau = 1. / (alpha + beta)
  return tau, alpha * tau

def Fn(V): # Func for K channel
  if V != -50.0: alpha = 0.01 * (V + 50.0) / (1-exp(-(V + 50.0)/10))
  else:          alpha = 0.1
  beta = 0.125 * exp(-0.0125*(V+60))
  tau = 1./(alpha + beta)
  return tau, alpha * tau




COUNT = 0

def increment():
    global COUNT
    COUNT = COUNT+1

def add(a,b): # fastest
  return [x+y for x,y in zip(a,b)]



class Alpha: # alpha function for current injection
  def __init__(self, delay, tau, amplitude):
    self.Delay = delay
    self.Tau = tau
    self.Amplitude = amplitude

  def _func(self, t):
    t = t - self.Delay
    if   t < 0: return 0.0
    else:       return 2.7183*self.Amplitude/self.Tau * t *exp(-t/self.Tau)

class Rect: # Rectangular function for current, voltage and concentration clampers

  def __init__(self, delay, width, amplitude):
    self.Delay = delay
    self.Width = width
    self.Amplitude = amplitude
    self.Tau = 0.01

  def _func(self, t):
    t = t - self.Delay
    if   t < 0:          return 0.0
    elif t < self.Width: return self.Amplitude
    else:                return 0.0

class Train:
  """ a train of Rectangular pulses"""
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
  """Concentration Clamper"""
  def __init__(self, ligand):
    self.Command  = None
    self.Waveform = None
    self.Ligand   = ligand
    self.Tag      = 'Ligand'

  def connect(self, cmpt):
    if type(cmpt) is list:
      for cp in cmpt:
        cp.lClamp = self
    else:
      cmpt.lClamp = self

  def set_amplitude(self, val):
    self.Waveorm.Amplitude = val

  def set_width(self, val):
    self.Waveorm.Width = val



class IClamper:
  """Current Clamper"""
  def __init__(self):
    self.Command  = None
    self.Waveform = None
    self.Tag      = 'Current'

  def connect(self, cmpt):
    if type(cmpt) is list:
      for cp in cmpt:
        cp.iClamp = self
    else:
      cmpt.iClamp = self

  def set_amplitude(self, val):
    self.Waveorm.Amplitude = val

  def set_width(self, val):
    self.Waveorm.Width = val


class VClamper:
  """Voltage clamper"""
  def __init__(self):
    self.Command  = None
    self.Waveform = None
    self.Baseline = -60
    self.Tag      = 'Voltage'
    self.Jc     = None # to store capacity current
    self.Jn     = None # to store cytosolic current
    self.Jm     = None # to store transmembrane current
    self.Jp     = None # Jc + Jm + Jn

  def connect(self, cmpt):
    if type(cmpt) is list:
      for cp in cmpt:
        cp.vClamp = self
    else:
      cmpt.vClamp = self

  def set_baseline(self, val):
    self.Waveorm.Baseline = val

  def set_amplitude(self, val):
    self.Waveorm.Amplitude = val

  def set_width(self, val):
    self.Waveorm.Width = val

  def save(self, filename):
    N = len(self.Jp)
    f = open(filename,'w')
    for k in range(N):
      s = '%7.5f %7.5f %7.5f %7.5f\n'%(self.Jp[k],self.Jm[k],self.Jn[k],self.Jc[k])
      f.write(s)
    f.close()

"""
class alphaConductance:"""

class LeakChannel:
  """Leak conductance"""
  def __init__(self, gMax=0.03, ER=-60):
    self.gMax = gMax
    self.ER   = ER
    self.Tag = 'Leak'

  def _gating(self): return 1.0
  def _update_gate(self, V, dt): return
  def _update_gMax(self,dt): return
  def store_gates(self, mem=0): return
  def load_gates(self, mem=0): return
  def show(self): pass


class KChannel:
  """Classical potassium channel"""
  def __init__(self, f=Fn, gMax=0.18, ER=-90):
    self._Func = f
    self.n0 = 0.299
    self.n  = self.n0
    self.gMax = gMax
    self.ER   = ER
    self.Fg   = gMax
    self.Rg   = 0.01
    self.gM   = None  # for storage of the gMax
    self.Tag  = 'KDR'

  def _gating(self):
    return self.n**4

  def _update_gate(self, V, dt):
    k, f = self._Func(V)
    self.n = exp(-dt/k) * (self.n-f) + f

  def _update_gMax(self, dt):
    gf = self.Fg * self.n**4
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

  def store_gates(self, mem=0):
    self.n0 = self.n

  def load_gates(self, mem = 0):
    self.n = self.n0

  def show(self):
    print(self.Tag)
    print('n =  %4.3f'% (self.n))

  def plot(self):
    V = arange(-90, 50, 0.5)
    N = len(V)
    tau = [0.0] * N
    inf = [0.0] * N
    for v in V:
      tau[i], inf[i] = self._Func(v)
    fig = pl.figure()
    ax1 = fig.add_subplot(121)
    ax1.plot(V, tau)
    ax2 = fig.add_subplot(122)
    ax2.plot(V, inf)
    #xlabel('mV')


class NaChannel:
  """Classical sodium channel"""
  def __init__(self, fm, fh, gMax, ER):
    self._Func1 = fm
    self._Func2 = fh
    self.m  = 0.046
    self.h  = 0.638
    self.m0 = 0.046
    self.h0 = 0.638
    self.gMax = gMax
    self.Fg   = gMax
    self.ER   = ER
    self.Rg   = 0.01
    self.gM   = None   # for storage of the gMax
    self.Tag  = 'NaC'

  def _gating(self):
    return (self.m**3) * self.h

  def _update_gate(self, V, dt):
    km, fm = self._Func1(V)
    kh, fh = self._Func2(V)
    self.m = exp(-dt/km) * (self.m-fm) + fm
    self.h = exp(-dt/kh) * (self.h-fh) + fh

  def show(self):
    print(self.Tag)
    print('m =  %4.3f' % (self.m))
    print('h =  %4.3f' % (self.h))

  def _update_gMax(self, dt):
    gf = self.Fg * self.m**3 * self.h
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

  def store_gates(self, mem = 0):
    if mem == 0: self.m0, self.h0 = self.m, self.h
    else: self.m1, self.h1 = self.m, self.h

  def load_gates(self, mem = 0):
    if mem == 0: self.m, self.h = self.m0, self.h0
    else: self.m, self.h = self.m1, self.h1



class Compartment:
  """
  Compartment on which neurons are built.
  """
  is_Compartment = 1

  def __init__(self, diameter = 20, length = 50):
    self.V0    = -61.237
    self.Vi    = self.V0
    self.Cm    = 0.01
    self.Gx    = 700.0
    self.Length = length
    self.Diameter = diameter
    self.iClamp = None
    self.vClamp = None
    self.lClamp = None
    self.Parent = None
    self.Child = []
    self.ID = COUNT
    increment()

    self.Vm    = None # for storage of the membrane potential
    self.channel_list = [] 

  def get_Jion(self, V):
    Jion = 0.0
    for ch in self.channel_list:
      g = ch.gMax * ch._gating()
      Jion += g * (V - ch.ER)
    return Jion

  def add_channel(self, Channel):
    if Channel.Tag == 'LGIC':
      Channel.binding = copy.deepcopy(Channel.binding)
      ch = LGIC(Channel.transit,Channel.binding,Channel.gMax,Channel.ER) # deal with python3.5
    else:
      ch = copy.deepcopy(Channel) #Python3.5???
    self.channel_list.append(ch)
    ch.ptr = self
    if   ch.Tag == 'NaC': self.NaC = ch
    elif ch.Tag == 'KDR': self.KDR = ch
    elif ch.Tag == 'Leak': self.Leak = ch
    elif ch.Tag == 'LGIC': self.LGIC = ch
    return ch

  def add_channels(self, Channels):
    if type(Channels) is not list:
      print ('Channel not added')
      return
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
      print ('Error:')
      print ('Connection failed.\n Do not connect a cylinder to a sphere.')
      sys.exit(1)
    self.Child.append(cmpt)
    cmpt.Parent = self

  def attached_to(self, cmpt):
    if not self.Length:
      print ('Error:')
      print ('Attachment failed.\n Do not attach a sphere to a cylinder.')
      sys.exit(1)
    cmpt.Child.append(self)
    self.Parent = cmpt

  def _update_Vm(self, dt):
    if self.channel_list == []: return self.Vi + dt * self.Jn / self.Cm
    sigma_g  = sigma_gE = 0.0
    for ch in self.channel_list:
      g = ch.gMax * ch._gating()
      sigma_g += g
      sigma_gE += g * ch.ER
    k = sigma_g / self.Cm
    Vf = (self.Jn + sigma_gE) / sigma_g
    return exp(-k*dt) * (self.Vi-Vf) + Vf

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

  def __init__(self, prep):
    if type(prep) is list:
      self.cmptList = prep
    else:
      self.cmptList = [prep]

    self.TotalTime = 0.0
    for cp in self.cmptList:
      cp.Vi = cp.V0
      cp._calc_surface()
      cp._calc_gx()

    N = len(self.cmptList)
    print ('_________________________________________\n')
    if N > 1:
      print ('%i compartments supplied\n'%(N))
    else:
      print ('%i compartment supplied\n'%(N))

  def save(self, filename):
    N = len(self.T)
    f = open(filename,'w')
    for cp in self.cmptList:
      f.write('# ID = %i\n' %(cp.ID))
      for k in range(N):
        f.write('%8.6f\n' %(cp.Vm[k]))
    f.close()

  def endpoint_info(self):
    for cp in self.cmptList:
      cp.show()
      for ch in cp.channel_list: ch.show()

  def run(self, t, dt):
    steps = int(t/dt)
    self.T = [i*dt+self.TotalTime for i in range(steps)]
    self.TotalTime += t

    print ("Integration in process")
    print ('Number of steps: %i'% (steps))
    if dt < 0.001:
      print ("Step length: %5.4f ms" % (dt))
    else:
      print ("Step length: %4.3f ms" % (dt))
    print ('Wait!')
    print ('...')

    self.L_Clampers = []
    for cp in self.cmptList:
      if cp.lClamp: self.L_Clampers.append(cp.lClamp)

    ### prepare to run
    self.V_Compartments = []
    self.L_Compartments = []
    self.I_Compartments = []
    self.F_Compartments = [] # free, no clampers attached
    for cp in self.cmptList:
      if   cp.vClamp and cp.vClamp not in self.V_Compartments:
        self.V_Compartments.append(cp)
      elif cp.lClamp and cp.lClamp not in self.L_Compartments:
        self.L_Compartments.append(cp)
      elif cp.iClamp and cp.iClamp not in self.I_Compartments:
        self.I_Compartments.append(cp)
      else:
        self.F_Compartments.append(cp)

    self.N_Compartments = self.I_Compartments + self.L_Compartments + self.F_Compartments

    for cp in self.cmptList: # check for missing values
      for ch in cp.channel_list:
        if ch.gMax == None:
          print ('Ion channel %s in Compartment %i doen not have gMax value')
          sys.exit(1)
        if ch.ER == None:
          print ('Ion channel %s in Compartment %i doen not have ER value')
          sys.exit(1)

    for cp in self.cmptList: # prepare the storages based on run parameters
      cp.Vm = array.array('f',[cp.V0]) * steps
      for ch in cp.channel_list:
        ch.gM = array.array('f',[0]) * steps

    for cp in self.cmptList:
      if cp.vClamp:
        cp.vClamp.Jm = array.array('f',[0]) * steps
        cp.vClamp.Jn = array.array('f',[0]) * steps
        cp.vClamp.Jc = array.array('f',[0]) * steps
        cp.vClamp.Jp = array.array('f',[0]) * steps

    for cp in self.cmptList: # generate command from specified waveform
      ic = cp.iClamp
      if ic:
        ic.Command = array.array('f',[0]) * steps
        if ic.Waveform != None:
          for i in range(steps):
            ic.Command[i] = ic.Waveform._func(self.T[i])

      lc = cp.lClamp
      if lc:
        lc.Command = array.array('f',[0]) * steps
        if lc.Waveform != None:
          for i in range(steps):
            lc.Command[i] = lc.Waveform._func(i*dt)

      vClmp = cp.vClamp
      if vClmp:
        vClmp.Command = array.array('f',[0]) * steps
        if vClmp.Waveform != None:
          for i in range(steps):
            vClmp.Command[i] = vClmp.Waveform._func(i*dt) + vClmp.Baseline
        vClmp.dot = [0.] + [ (vClmp.Command[i]-vClmp.Command[i-1])/dt for i in range(1,steps) ]

    ###

    t0 = time.time() # looping begins
    for step_num in range(steps):
      for cp in self.cmptList:
        if cp.iClamp:      # we don't directly access the clamper
          cp.Ja = cp.iClamp.Command[step_num]
        else:
          cp.Ja = 0.0

        if cp.Parent:      # get Jx (axial current density)
          cp.Jx = (cp.Vi - cp.Parent.Vi) * cp.gx/cp.Surface
        else:
          cp.Jx = 0.0

      for cp in self.cmptList:       # sum up non-Im current, or cytoplasmuc current
        cp.Jn = cp.Ja - cp.Jx
        for chd in cp.Child:
          cp.Jn += chd.Jx

      for lc in self.L_Clampers: # We access the clamper through Experiment
        lc.Ligand.Conc = lc.Command[step_num]

      # Update_Start ========
      for cp in self.N_Compartments:
        cp.Vi = cp._update_Vm(dt)

        for ch in cp.channel_list:
          ch._update_gate(cp.Vi,dt)
          #ch._update_gMax(dt)

      for cp in self.V_Compartments: # Voltage-Clamp block
        vClmp = cp.vClamp
        cp.Vi = vClmp.Command[step_num] # ???
        for ch in cp.channel_list:
          ch._update_gate(cp.Vi,dt)

        vClmp.Jc[step_num] = cp.Cm * vClmp.dot[step_num]
        vClmp.Jm[step_num] = cp.get_Jion(vClmp.Command[step_num])
        vClmp.Jn[step_num] = cp.Jn # here I corrected a mistake in the submitted version
        vClmp.Jp[step_num] = cp.Jn + vClmp.Jc[step_num] + vClmp.Jm[step_num]

      # Update_End ========
      for cp in self.cmptList: # Store membrane potentials
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
  This class defines receptor-operated channels (LGIC) according to Markov schemes.
  """

  def __init__(self, transit, binding, gMax, ER):
    self.transit = transit # keep it for copying the channel
    self.binding = binding # keep it for copying the channel
    self.State = 0
    self.LigandList = binding.keys()
    self.do_parsing(transit, binding)
    self.N = len(self.RateMatrix)    #self.N = len(self.RM)
    self.StateProb = [0] * self.N # create and initiate state probability vector
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

  def _gating(self):
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



class NaChannel_MHMJ:
  """MHMJ sodium channel"""
  def __init__(self, fm, fh, fi, gMax, ER):
    self._Func1 = fm
    self._Func2 = fh
    self._Func3 = fi
    self.m  = 0.046
    self.h  = 0.638
    self.i  = 0.638
    self.m0 = 0.046
    self.h0 = 0.638
    self.i0 = 0.638
    self.gMax = gMax
    self.Fg   = gMax
    self.ER   = ER
    self.Rg   = 0.01
    self.gM   = None
    self.Tag  = 'MHMJ'

  def _gating(self):
    return (self.m**3) * self.h * self.i

  def _update_gate(self, V, dt):
    km, fm = self._Func1(V)
    kh, fh = self._Func2(V)
    ki, fi = self._Func3(V)
    self.m = exp(-dt/km) * (self.m-fm) + fm
    self.h = exp(-dt/kh) * (self.h-fh) + fh
    self.i = exp(-dt/ki) * (self.i-fi) + fi

  def _update_gMax(self, dt):
    gf = self.Fg * self.m**3 * self.h * self.i
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

  def set_m(self, m):
    self.m = m
  def set_h(self, h):
    self.h = h
  def set_i(self, i):
    self.i = i

  def show(self):
    print(self.Tag)
    print('m', self.m)
    print('h', self.h)
    print('i', self.i)

  def plot(self):
    V = arange(-90, 50, 0.5)
    N = len(V)
    tau_a, inf_a = [0.0] * N, [0.0] * N
    tau_b, inf_b = [0.0] * N, [0.0] * N
    tau_i, inf_i = [0.0] * N, [0.0] * N
    for i in range(N):
      tau_a[i], inf_a[i] = self._Func1(V[i])
      tau_b[i], inf_b[i] = self._Func2(V[i])
      tau_i[i], inf_i[i] = self._Func3(V[i])

    fig = pl.figure()
    ax1 = fig.add_subplot(321)
    ax1.plot(V, inf_a)
    ax2 = fig.add_subplot(322)
    ax2.plot(V, tau_a)
    ax3 = fig.add_subplot(323)
    ax3.plot(V, inf_b)
    ax4 = fig.add_subplot(324)
    ax4.plot(V, tau_b)
    ax5 = fig.add_subplot(325)
    ax5.plot(V, inf_i)
    ax6 = fig.add_subplot(326)
    ax6.plot(V, tau_i)
    ax1.set_xlim([-90,50])
    ax2.set_xlim([-90,50])
    ax3.set_xlim([-90,50])
    ax4.set_xlim([-90,50])
    ax5.set_xlim([-90,50])
    ax6.set_xlim([-90,50])
    ax5.set_xlabel('mV')
    ax6.set_xlabel('mV')
    ax1.set_ylabel('inf')
    ax3.set_ylabel('inf')
    ax5.set_ylabel('inf')
    ax2.set_ylabel('tau')
    ax4.set_ylabel('tau')
    ax6.set_ylabel('tau')


class HH_aAbB:
  """Ion channels with (a**N) * (b**M) gating"""
  def __init__(self, fa, Na, fb, Nb, gMax, ER):
    self._Func1 = fa
    self._Func2 = fb
    self.Na = Na
    self.Nb = Nb
    self.a  = 0.046
    self.b  = 0.638
    self.a0 = 0.046
    self.b0 = 0.638
    self.gMax = gMax
    self.ER = ER
    self.Tag = 'NaYb'

  def _gating(self):
    return (self.a**self.Na) * (self.b**self.Nb)

  def _update_gate(self, V, dt):
    ka, fa = self._Func1(V)
    kb, fb = self._Func2(V)
    self.a = exp(-dt/ka) * (self.a-fa) + fa
    self.b = exp(-dt/kb) * (self.b-fb) + fb

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)

  def plot(self):
    V = arange(-90, 50, 0.5)
    N = len(V)
    tau_a = [0.0] * N
    inf_a = [0.0] * N
    tau_b = [0.0] * N
    inf_b = [0.0] * N
    for i in range(N):
      tau_a[i], inf_a[i] = self._Func1(V[i])
      tau_b[i], inf_b[i] = self._Func2(V[i])

    fig = pl.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(V, tau_a)
    ax1.plot(V, tau_b)
    ax2 = fig.add_subplot(212)
    ax2.plot(V, inf_a)
    ax2.plot(V, inf_b)
    #xlabel('mV')

class HH_aA:
  """Ion channels with a**N gating"""
  def __init__(self, f, A, gMax, ER):
    self._Func = f
    self.N = A
    self.a0 = 0.299
    self.a = self.a0
    self.gMax = gMax
    self.ER = ER
    self.Tag = 'Na0b'

  def _gating(self):
    return self.a**self.N

  def _update_gate(self, V, dt):
    k, f = self._Func(V)
    self.a = exp(-dt/k) * (self.a-f) + f

  def show(self):
    print(self.Tag)
    print('a', self.a)

  def plot(self):
    V = arange(-90, 50, 0.5)
    N = len(V)
    tau_a = [0.0] * N
    inf_a = [0.0] * N
    for i in range(N):
      tau_a[i], inf_a[i] = self._Func1(V[i])

    fig = pl.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(V, tau_a)
    ax2 = fig.add_subplot(212)
    ax2.plot(V, inf_a)
    #xlabel('mV')

class HH_aAb:
  """Ion channels with (a**A) * b gating"""
  def __init__(self, fa, A, fb, gMax, ER):
    self._Func1 = fa
    self._Func2 = fb
    self.Na = A
    self.a  = 0.046
    self.b  = 0.638
    self.a0 = 0.046
    self.b0 = 0.638
    self.gMax = gMax
    self.ER = ER
    self.Tag = 'NaYb'

  def _gating(self):
    return (self.a**self.Na) * self.b

  def _update_gate(self, V, dt):
    ka, fa = self._Func1(V)
    kb, fb = self._Func2(V)
    self.a = exp(-dt/ka) * (self.a-fa) + fa
    self.b = exp(-dt/kb) * (self.b-fb) + fb

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)

  def plot(self):
    V = arange(-90, 50, 0.5)
    N = len(V)
    tau_a = [0.0] * N
    inf_a = [0.0] * N
    tau_b = [0.0] * N
    inf_b = [0.0] * N
    for i in range(N):
      tau_a[i], inf_a[i] = self._Func1(V[i])
      tau_b[i], inf_b[i] = self._Func2(V[i])

    fig = pl.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(V, tau_a)
    ax1.plot(V, tau_b)
    ax2 = fig.add_subplot(212)
    ax2.plot(V, inf_a)
    ax2.plot(V, inf_b)
    #xlabel('mV')

class HH_a:
  """Ion channels with one gate type"""
  def __init__(self, f, gMax, ER):
    self._Func = f
    self.a  = 0.046
    self.a0 = 0.046
    self.gMax = gMax
    self.ER = ER
    self.Tag = '1a0b'

  def _gating(self):
    return self.a

  def _update_gate(self, V, dt):
    ka, fa = self._Func(V)
    self.a = exp(-dt/ka) * (self.a-fa) + fa

  def set_a(self,a):
    self.a = a

  def show(self):
    print(self.Tag)
    print('a', self.a)

  def plot(self):
    V = arange(-90, 50, 0.5)
    N = len(V)
    tau = [0.0] * N
    inf = [0.0] * N
    for i in range(N):
      tau[i], inf[i] = self._Func(V[i])
    fig = pl.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(V, tau)
    ax2 = fig.add_subplot(212)
    ax2.plot(V, inf)
    #xlabel('mV')

class HH_ab:
  """Ion channels with a gating"""
  def __init__(self, fa, fb, gMax, ER):
    self._Func1 = fa
    self._Func2= fb
    self.a  = 0.046
    self.a0 = 0.046
    self.b  = 0.046
    self.b0 = 0.046
    self.gMax = gMax
    self.ER = ER
    self.Tag = '1a0b'

  def _gating(self):
    return self.a * self.b

  def _update_gate(self, V, dt):
    ka, fa = self._Func1(V)
    kb, fb = self._Func2(V)
    self.a = exp(-dt/ka) * (self.a-fa) + fa
    self.b = exp(-dt/kb) * (self.b-fb) + fb

  def set_a(self,a):
    self.a = a

  def set_b(self,b):
    self.b = b

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)

  def plot(self):
    V = arange(-90, 50, 0.5)
    N = len(V)
    tau_a = [0.0] * N
    inf_a = [0.0] * N
    tau_b = [0.0] * N
    inf_b = [0.0] * N
    for i in range(N):
      tau_a[i], inf_a[i] = self._Func1(V[i])
      tau_b[i], inf_b[i] = self._Func2(V[i])

    fig = pl.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(V, tau_a)
    ax1.plot(V, tau_b)
    ax2 = fig.add_subplot(212)
    ax2.plot(V, inf_a)
    ax2.plot(V, inf_b)
    #xlabel('mV')

class HH_abc:
  """MHMJ sodium channel"""
  def __init__(self, fa, fb, fc, gMax, ER):
    self._Func1 = fa
    self._Func2 = fb
    self._Func3 = fc
    self.m  = 0.046
    self.h  = 0.638
    self.i  = 0.638
    self.m0 = 0.046
    self.h0 = 0.638
    self.i0 = 0.638
    self.gMax = gMax
    self.Fg   = gMax
    self.ER   = ER
    self.Rg   = 0.01
    self.gM   = None
    self.Tag  = 'MHMJ'

  def _gating(self):
    return self.m * self.h * self.i

  def _update_gate(self, V, dt):
    km, fm = self._Func1(V)
    kh, fh = self._Func2(V)
    ki, fi = self._Func3(V)
    self.m = exp(-dt/km) * (self.m-fm) + fm
    self.h = exp(-dt/kh) * (self.h-fh) + fh
    self.i = exp(-dt/ki) * (self.i-fi) + fi

  def _update_gMax(self, dt):
    gf = self.Fg * self.m**3 * self.h * self.i
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

  def set_m(self, m):
    self.m = m

  def set_h(self, h):
    self.h = h

  def set_i(self, i):
    self.i = i

  def show(self):
    print(self.Tag)
    print('m', self.m)
    print('h', self.h)
    print('i', self.i)

  def plot(self):
    V = arange(-90, 50, 0.5)
    N = len(V)
    tau_a, inf_a = [0.0] * N, [0.0] * N
    tau_b, inf_b = [0.0] * N, [0.0] * N
    tau_i, inf_i = [0.0] * N, [0.0] * N
    for i in range(N):
      tau_a[i], inf_a[i] = self._Func1(V[i])
      tau_b[i], inf_b[i] = self._Func2(V[i])
      tau_i[i], inf_i[i] = self._Func3(V[i])

    fig = pl.figure()
    ax1 = fig.add_subplot(321)
    ax1.plot(V, tau_a)
    ax2 = fig.add_subplot(322)
    ax2.plot(V, inf_a)
    ax3 = fig.add_subplot(323)
    ax3.plot(V, tau_b)
    ax4 = fig.add_subplot(324)
    ax4.plot(V, inf_b)
    ax5 = fig.add_subplot(325)
    ax5.plot(V, tau_i)
    ax6 = fig.add_subplot(326)
    ax6.plot(V, inf_i)
    #xlabel('mV')




GABA = Ligand()
Glu = Ligand()

gabar_transit = {
                 'C0': {'C1': 20.0},
                 'C1': {'C0': 4.6, 'O1': 3.3, 'C2': 10.0},
                 'C2': {'C1': 9.2, 'O2': 10.6},
                 'O1': {'C1': 9.8},
                 'O2': {'C2': 0.41}
                }

gabar_binding = {
                 GABA:
                 {
                   'C0': 'C1',
                   'C1': 'C2'
                 }
                }


nmdar_transit = {
                 'C0': {'C1': 5.0},
                 'C1': {'C0': 0.0129, 'C2': 5.0},
                 'C2': {'C1': 0.0129, 'D': 0.0084, 'O': 0.0465},
                 'D':  {'C2': 0.0068},
                 'O':  {'C2': 0.0738}
                }

nmdar_binding = {Glu:
                 {'C0': 'C1',
                  'C1': 'C2'}
                }


# gate_a = (Fa, A), (Fm,3), (Fn,4)
# gate_b = (Fb, B)
# gate_c = (Fc, C), (None,0)

def HHChannel(gate_a, gate_b, gate_c, gMax=None, ER=None):
  Fa, A = gate_a[0], gate_a[1]
  Fb, B = gate_b[0], gate_b[1]
  Fc, B = gate_c[0], gate_c[1]
  if C == 0 and B == 0 and A == 1: return HH_a(Fa, gMax, ER) #
  if C == 0 and B == 0 and A >  1: return HH_aA(Fa, A, gMax, ER)#
  if C == 0 and B == 1 and A == 1: return HH_ab(Fa, Fb, gMax, ER)#
  if C == 0 and B == 1 and A >  1: return HH_aAb(Fa, A, Fb, gMax, ER)
  if C == 0 and B >  1 and A >  1: return HH_aAbB(Fa, A, Fb, B, gMax, ER)
  if C == 1 and B == 1 and A == 1: return HH_abc(Fa, Fb, Fc, gMax, ER)
"""
  if C == 1 and B == 1 and A >  1: return HH_aAbc(Fa,Fb)
  if C == 1 and B  > 1 and A >  1: return HH_aAbBc(Fa,Fb)
  if C >  1 and B  > 1 and A >  1: return HH_aAbBcC(Fa,Fb)"""



NaC = NaChannel(Fm, Fh, gMax=0.6, ER=50)    # the unit for gMax is ?S/um2
KDR = KChannel(Fn, gMax=0.18, ER=-90)
gL  = LeakChannel(gMax=0.03,ER=-60)



from math import exp
import copy
import array
import time
import sys

def Fm(V): # Func for Na channel
  if V != -35.0: am = 0.1 * (V + 35.0) /(1-exp(-(V + 35.0)/10.0))
  else:          am = 1.0
  bm = 4 * exp(-(V+60)/18)
  return am, bm

def Fh(V): # Func for Na channel
  ah = 0.07 * exp(-0.05*(V+60.0))
  bh = 1 / (1+exp(-0.1*(V+30.0)))
  return ah, bh

def Fn(V): # Func for K channel
  if V != -50.0: an = 0.01 * (V + 50.0) / (1-exp(-(V + 50.0)/10))
  else:          an = 0.1
  bn = 0.125 * exp(-0.0125*(V+60))
  return an, bn


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
    if   t < 0:          return 0.0
    else:                return 2.7183*self.Amplitude/self.Tau * t *exp(-t/self.Tau)

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


class CClamper: # Concentration Clamper

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



class IClamper: # Current Clamper

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


class VClamper: #Voltage clamper

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


class LeakChannel:
  """Leak conductance"""
  def __init__(self, gMax=0.03, ER=-60):
    self.gMax = gMax
    self.ER   = ER
    self.Tag = 'Leak'

  def _gating(self): return 1.0
  def _update(self, V, h): return
  def store_gates(self, mem=0): return
  def load_gates(self, mem=0): return
  def _update_gMax(self,h): return


class KChannel:
  """Classical potassium channel"""
  def __init__(self, gMax=0.18, ER=-90, func=Fn):
    self._Func = func
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

  def _update(self, V, h):
    alpha, beta = self._Func(V)
    k = alpha + beta
    nf = alpha / k
    self.n = exp(-h*k) * (self.n-nf) + nf

  def _update_gMax(self, dt):
    gf = self.Fg * self.n**4
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

  def store_gates(self, mem=0):
    self.n0 = self.n

  def load_gates(self, mem = 0):
    self.n = self.n0

class NaChannel:
  """Classical sodium channel"""

  def __init__(self, gMax=0.6, ER=50, fm=Fm, fh=Fh):
    self._Func1 = fm
    self._Func2 = fh
    self.m  = 0.046
    self.h  = 0.638
    self.m0 = 0.046
    self.h0 = 0.638
    self.m1 = 0.046
    self.h1 = 0.638
    self.gMax = gMax
    self.Fg   = gMax
    self.ER   = ER
    self.Rg   = 0.01
    self.gM   = None   # for storage of the gMax
    self.Tag  = 'NaC'

  def _gating(self):
    return (self.m**3) * self.h

  def _update(self, V, dt):
    am, bm = self._Func1(V)
    ah, bh = self._Func2(V)
    km = am + bm
    fm = am / km
    kh = ah + bh
    fh = ah / kh
    self.m = exp(-dt*km) * (self.m-fm) + fm
    self.h = exp(-dt*kh) * (self.h-fh) + fh

  def _update_gMax(self, dt):
    gf = self.Fg * self.m**3 * self.h
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

  def store_gates(self, mem = 0):
    if mem == 0: self.m0, self.h0 = self.m, self.h
    else: self.m1, self.h1 = self.m, self.h

  def load_gates(self, mem = 0):
    if mem == 0: self.m, self.h = self.m0, self.h0
    else: self.m, self.h = self.m1, self.h1



class HH1V:

  def __init__(self, func, N):
    self._Func = func
    self.N = N
    self.n0 = 0.299
    self.n = self.n0
    self.Tag = 'HH1V'
    self.gMax = 0.18
    self.ER = -90

  def _gating(self):
    return self.n**self.N

  def _update(self, V, dt):
    alpha, beta = self._Func(V)
    tau = 1.0/(alpha+beta)
    ff = alpha * tau
    self.n = exp(-dt/tau) * (self.n-ff) + ff

  def store_gates(self, mem=0):
    self.n0 = self.n

  def load_gates(self, mem = 0):
    self.n = self.n0


class HH2V:
  is_Channel = 1
  is_HHChannel = 1
  is_HHV = 1

  def __init__(self, fm, Nm, fh, Nh):
    self._Func1 = fm
    self._Func2 = fh
    self.Nm = Nm
    self.Nh = Nh
    self.m  = 0.046
    self.h  = 0.638
    self.m0 = 0.046
    self.h0 = 0.638
    self.m1 = 0.046
    self.h1 = 0.638
    self.gMax = 0.6 # ???
    self.ER = 50    # ???
    self.Tag = 'HH2V'

  def _gating(self):
    return (self.m**self.Nm) * (self.h**self.Nh)

  def _update(self, V, h):
    am, bm = self._Func1(V)
    ah, bh = self._Func2(V)
    tm = 1.0/(am+bm)
    fm = am * tm
    th = 1.0/(ah+bh)
    fh = ah * th
    self.m = exp(-h/tm) * (self.m-fm) + fm
    self.h = exp(-h/th) * (self.h-fh) + fh

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


class Experiment:

  def __init__(self, prep):
    if type(prep) is list:
      self.cmptList = prep
    else:
      self.cmptList = [prep]

    aNumb = []
    for cp in self.cmptList:
      if cp.ID in aNumb:
        print ('\nWarning: Compartments share same ID = %i.\n'%(cp.ID))
      else:
        aNumb.append(cp.ID)

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

  def run(self, t, dt):
    steps = self.nStep = int(t/dt)
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

    self.T = [i*dt for i in range(steps)]
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

    for cp in self.cmptList: # generate Ja pulse
      ic = cp.iClamp
      if ic:
        ic.Command = array.array('f',[0]) * steps
        if ic.Waveform != None:
          for i in range(steps):
            ic.Command[i] = ic.Waveform._func(i*dt) # why +=

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
          ch._update(cp.Vi,dt)
          #ch._update_gMax(dt)

      for cp in self.V_Compartments: # Voltage-Clamp block
        vClmp = cp.vClamp
        cp.Vi = vClmp.Command[step_num] # ???
        for ch in cp.channel_list:
          ch._update(cp.Vi,dt)

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


gL = LeakChannel()
NaC = NaChannel()
KDR = KChannel()

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

  def _update(self, V, h):
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

  def _gating(self):
    gate = 0.0
    for k in self.OpenIndex: gate += self.StateProb[k]
    return gate




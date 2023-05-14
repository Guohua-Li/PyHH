from math import exp#

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



NaC = NaChannel(Fm, Fh, gMax=0.6, ER=50) # the unit for gMax is nS/um2
KDR = KChannel(Fn, gMax=0.18, ER=-90)
gL  = LeakChannel(gMax=0.03,ER=-60)


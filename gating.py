from math import exp


import matplotlib.pyplot as plt

def Fm(V): # Na channel, m gate
  v = V + 35 # 35 -> 40
  alpha = 1.0 if v == 0 else 0.1 * v /(1-exp(-0.1*v))
  beta = 4.0 * exp(-(V+60)/18) # 60->65
  tau = 1.0/(alpha + beta)
  return tau, alpha * tau

def Fh(V): # Na channel, h gate
  alpha = 0.07 * exp(-0.05*(V+60.0)) # 60 -> 65
  beta  = 1.0 / (1.0+exp(-0.1*(V+30.0))) # 30 -> 35
  tau = 1.0 / (alpha + beta)
  return tau, alpha * tau

def Fn(V): # K channel, n gate
  v = V + 50.0
  alpha = 0.1 if v == 0 else 0.01 * v / (1.0-exp(-0.1*v)) # 50 -> 55
  beta = 0.125 * exp(-0.0125*(V+60)) # 60 -> 65
  tau = 1.0/(alpha + beta)
  return tau, alpha * tau



def arange(a,b,step):
  N = int((b-a)/step)
  return [a+i*step for i in range(N)]


def _plt(f1, f2=None, f3=None, Vi=None):
  if Vi == None:
    Start, Finish = -100, 80
  else:
    Start, Finish = Vi[0], Vi[1]
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
  fig = plt.figure()
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
  plt.show()


class LeakChannel:
  """
  Leak conductance
  """
  def __init__(self):
    self.Tag = 'Leak'
    self.ptr = None

  def _ProbOpen(self): return 1.0
  def _update_gate(self, V, dt): return
  def store_gates(self, mem=0): return
  def load_gates(self, mem=0): return
  def show(self): pass

class NaChannel:
  """
  Classical sodium channel
  """
  def __init__(self, fa, fb):
    self._fa = fa
    self._fb = fb
    self.a0 = 0.046
    self.b0 = 0.638
    self.a  = self.a0
    self.b  = self.b0
    self.Tag  = 'NaT'
    self.ptr = None

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

  def plot(self):
    _plt(self._fa, self._fb)


class KChannel:
  """
  Classical potassium channel
  """
  def __init__(self, f = Fn):
    self._f = f
    self.a0 = 0.299
    self.a  = self.a0
    self.ptr = None
    self.Tag  = 'KDR'

  def _ProbOpen(self):
    return self.a**4

  def _update_gate(self, V, dt):
    k, f = self._f(V)
    self.a = exp(-dt/k) * (self.a-f) + f

  def show(self):
    print(self.Tag)
    print('a =  %4.3f'% (self.a))

  def plot(self):
    _plt(self._f)

NaT = NaChannel(Fm, Fh) # the unit for gMax is nS/um2
KDR = KChannel(Fn)
gL  = LeakChannel()


class HH_a:
  """
  gating: a
  """
  def __init__(self, f):
    self._f = f
    self.a0 = 0.046
    self.a  = self.a0
    self.Tag = 'HH_a'
    self.ptr = None

  def _ProbOpen(self):
    return self.a

  def _update_gate(self, V, dt):
    tau, inf = self._f(V)
    self.a = exp(-dt/tau) * (self.a-inf) + inf

  def set_a(self,a):
    self.a = a

  def show(self):
    print(self.Tag)
    print('a', self.a)

  def plot(self, Vm=None):
    _plt(self._f, Vi=Vm)


class HH_aAb:
  """
  gating: (a**A) * b
  """
  def __init__(self, fa, A, fb):
    self._fa = fa
    self._fb = fb
    self.A = A
    self.a0 = 0.046
    self.b0 = 0.638
    self.a  = self.a0
    self.b  = self.b0
    self.Tag = 'HH_a(A)b'
    self.ptr = None

  def couple(self, p):
    self.pool_pointer = p
    p.target = self

  def _ProbOpen(self):
    return (self.a**self.A) * self.b

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)


  def plot(self, Vm=None):
    _plt(self._fa, self._fb, Vi = Vm)


class HH_aAbB:
  """
  gating: (a**N) * (b**M)
  """
  def __init__(self, fa, A, fb, B):
    self._fa = fa
    self._fb = fb
    self.A = A
    self.B = B
    self.a0 = 0.046
    self.b0 = 0.638
    self.a  = 0.046
    self.b  = 0.638
    self.Tag = 'a(A)b(B)'
    self.ptr = None

  def _ProbOpen(self):
    return (self.a**self.A) * (self.b**self.B)

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)

  def plot(self, Vm=None):
    _plt(self._fa, self._fb, Vi=Vm)


class HH_aA:
  """
  gating: a**A
  """
  def __init__(self, f, A):
    self._Func = f
    self.N = A
    self.a0 = 0.299
    self.a = self.a0
    self.Tag = 'HH_a(A)'
    self.ptr = None

  def _ProbOpen(self):
    return self.a**self.N

  def _update_gate(self, V, dt):
    k, f = self._Func(V)
    self.a = exp(-dt/k) * (self.a-f) + f

  def show(self):
    print(self.Tag)
    print('a', self.a)

  def plot(self, Vm=None):
    _plt(self._fa, Vi=Vm)


class HH_ab:
  """
  gating: a*b
  """
  def __init__(self, fa, fb):
    self._fa = fa
    self._fb= fb
    self.a0 = 0.046
    self.b0 = 0.046
    self.a  = 0.046
    self.b  = 0.046
    self.Tag = 'HH_ab'
    self.ptr = None

  def _ProbOpen(self):
    return self.a * self.b

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb

  def set_a(self,a):
    self.a = a

  def set_b(self,b):
    self.b = b

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)

  def plot(self, Vm=None):
    _plt(self._fa, self._fb, Vi = Vm)



class HH_abc:
  """
  gating: a*b*c
  """
  def __init__(self, fa, fb, fc):
    self._fa = fa
    self._fb = fb
    self._fc = fc
    self.a0 = 0.046
    self.b0 = 0.638
    self.c0 = 0.638
    self.a  = 0.046
    self.b  = 0.638
    self.c  = 0.638
    self.Tag  = 'HH_abc'
    self.ptr = None

  def _ProbOpen(self):
    return self.a * self.b * self.c

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    tc, fc = self._fc(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    self.c = exp(-dt/tc) * (self.c-fc) + fc

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

  def plot(self, Vm=None):
    _plt(self._fa, self._fb, self._fc, Vi=Vm)


class HH_aAbc:
  """
  gating: (a**A) * b * c
  """
  def __init__(self, fa, A, fb, fc):
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
    self.Tag  = 'HH_aAbc'
    self.ptr = None

  def _ProbOpen(self):
    return (self.a**self.A) * self.b * self.c

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    tc, fc = self._fc(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    self.c = exp(-dt/tc) * (self.c-fc) + fc

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

  def plot(self, Vm=None):
    _plt(self._fa, self._fb, self._fc, Vi = Vm)


class HH_aAbBc:
  """
  gating: (a**A) * (b**B) * c
  """
  def __init__(self, fa, A, fb, B, fc):
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
    self.Tag  = 'HH_aAbBc'
    self.ptr = None

  def _ProbOpen(self):
    return (self.a**self.A) * (self.b**self.B) * self.c

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    tc, fc = self._fc(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    self.c = exp(-dt/tc) * (self.c-fc) + fc

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

  def plot(self, Vm=None):
    _plt(self._fa, self._fb, self._fc, Vi=Vm)


class HH_aAbBcC:
  """
  gating: (a**A) * (b**B) * (c**C)
  """
  def __init__(self, fa, A, fb, B, fc, C):
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
    self.Tag  = 'HH_a(A)b(B)c(C)'
    self.ptr = None

  def _ProbOpen(self):
    return (self.a**self.A) * (self.b**self.B) * (self.c**self.C)

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    tc, fc = self._fc(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    self.c = exp(-dt/tc) * (self.c-fc) + fc

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

  def plot(self, Vm=None):
    _plt(self._fa, self._fb, self._fc, Vi=Vm)



class NaChannel_MHMJ:
  """
  MHMJ sodium channel
  """
  def __init__(self, fm, fh, fi):
    self._fa = fm
    self._fb = fh
    self._fc = fi
    self.a  = 0.046
    self.b  = 0.638
    self.c  = 0.638
    self.a0 = 0.046
    self.b0 = 0.638
    self.c0 = 0.638
    self.Tag  = 'MHMJ'

  def _ProbOpen(self):
    return (self.a**3) * self.b * self.c

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    tc, fc = self._fc(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    self.c = exp(-dt/tc) * (self.c-fc) + fc

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

  def plot(self, Vm=None):
    _plt(self._fa, self._fb, self._fc, Vi=Vm)



def voltage_gated(gate_a = None, gate_b = None, gate_c = None):
  """
    if gate_a <> None, gate_a is gate_a tuple of two elements,
    the first one is the function, the second one is the power.
  """

  if gate_a == None and gate_b == None and gate_c == None:
    return LeakChannel()

  if gate_b == None and gate_c == None:
    Fa, A = gate_a[0], gate_a[1]
    if A == 1: return HH_a(Fa)
    if A >  1:
      if A == 4: # special case
        return KChannel(Fa)
      else:
        return HH_aA(Fa, A)#

  if gate_c == None:
    Fa, A = gate_a[0], gate_a[1]
    Fb, B = gate_b[0], gate_b[1]
    if B == 1 and A == 1: return HH_ab(Fa, Fb)#
    if B == 1 and A >  1:
      if A == 3: # special case
        return NaChannel(Fa, Fb)
      else:
        return HH_aAb(Fa, A, Fb)
    if B >  1 and A >  1: return HH_aAbB(Fa, A, Fb, B)

  Fa, A = gate_a[0], gate_a[1]
  Fb, B = gate_b[0], gate_b[1]
  Fc, C = gate_c[0], gate_c[1]
  if C == 1 and B == 1 and A == 1: return HH_abc(Fa,Fb,Fc)
  if C == 1 and B == 1 and A >  1: return HH_aAbc(Fa,A,Fb,Fc)
  if C == 1 and B  > 1 and A >  1: return HH_aAbBc(Fa,A, Fb, B, Fc)
  if C >  1 and B  > 1 and A >  1: return HH_aAbBcC(Fa, A, Fb, B, Fc, C)

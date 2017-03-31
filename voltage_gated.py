from math import exp
from pyhh import CaPool, LeakChannel, NaChannel, KChannel, _plt


class HH_a:
  """
  gating: a
  """
  def __init__(self, f, gMax, ER):
    self._f = f
    self.a0 = 0.046
    self.a  = self.a0
    self.gMax = gMax
    self.ER = ER
    self.Pools = []
    self.Tag = 'HH_a'
    self.ptr = None

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d
  def _ProbOpen(self):
    return self.a

  def _update_gate(self, V, dt):
    tau, inf = self._f(V)
    self.a = exp(-dt/tau) * (self.a-inf) + inf
    if len(self.Pools) > 0:     #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def set_a(self,a):
    self.a = a

  def show(self):
    print(self.Tag)
    print('a', self.a)

  def plot(self):
    _plt(self._f)


class HH_aAb:
  """
  gating: (a**A) * b
  """
  def __init__(self, fa, A, fb, gMax, ER):
    self._fa = fa
    self._fb = fb
    self.A = A
    self.a0 = 0.046
    self.b0 = 0.638
    self.a  = self.a0
    self.b  = self.b0
    self.gMax = gMax
    self.ER = ER
    self.Pools = []
    self.Tag = 'HH_a(A)b'
    self.ptr = None

  def couple(self, p):
    self.pool_pointer = p
    p.target = self

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return (self.a**self.A) * self.b

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    if len(self.Pools) > 0:
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)

  def plot(self):
    _plt(self._fa, self._fb)


class HH_aAbB:
  """
  gating: (a**N) * (b**M)
  """
  def __init__(self, fa, A, fb, B, gMax, ER):
    self._fa = fa
    self._fb = fb
    self.A = A
    self.B = B
    self.a0 = 0.046
    self.b0 = 0.638
    self.a  = 0.046
    self.b  = 0.638
    self.gMax = gMax
    self.ER = ER
    self.Pools = []
    self.Tag = 'a(A)b(B)'
    self.ptr = None

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d
  def _ProbOpen(self):
    return (self.a**self.A) * (self.b**self.B)

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb

    if len(self.Pools) > 0:    #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)

  def plot(self):
    _plt(self._fa, self._fb)


class HH_aA:
  """
  gating: a**A
  """
  def __init__(self, f, A, gMax, ER):
    self._Func = f
    self.N = A
    self.a0 = 0.299
    self.a = self.a0
    self.gMax = gMax
    self.ER = ER
    self.Pools = []
    self.Tag = 'HH_a(A)'
    self.ptr = None

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return self.a**self.N

  def _update_gate(self, V, dt):
    k, f = self._Func(V)
    self.a = exp(-dt/k) * (self.a-f) + f

    if len(self.Pools) > 0:    #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def show(self):
    print(self.Tag)
    print('a', self.a)

  def plot(self):
    _plt(self._fa)


class HH_ab:
  """
  gating: a*b
  """
  def __init__(self, fa, fb, gMax, ER):
    self._fa = fa
    self._fb= fb
    self.a0 = 0.046
    self.b0 = 0.046
    self.a  = 0.046
    self.b  = 0.046
    self.gMax = gMax
    self.ER = ER
    self.Pools = []
    self.Tag = 'HH_ab'
    self.ptr = None

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return self.a * self.b

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb

    if len(self.Pools) > 0:    #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def set_a(self,a):
    self.a = a

  def set_b(self,b):
    self.b = b

  def show(self):
    print(self.Tag)
    print('a', self.a)
    print('b', self.b)

  def plot(self):
    _plt(self._fa, self._fb)


class HH_abc:
  """
  gating: a*b*c
  """
  def __init__(self, fa, fb, fc, gMax, ER):
    self._fa = fa
    self._fb = fb
    self._fc = fc
    self.a0 = 0.046
    self.b0 = 0.638
    self.c0 = 0.638
    self.a  = 0.046
    self.b  = 0.638
    self.c  = 0.638
    self.gMax = gMax
    self.Fg   = gMax
    self.ER   = ER
    self.Rg   = 0.01
    self.gM   = None
    self.Pools = []
    self.Tag  = 'HH_abc'
    self.ptr = None

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return self.a * self.b * self.c

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    tc, fc = self._fc(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    self.c = exp(-dt/tc) * (self.c-fc) + fc

    if len(self.Pools) > 0:    #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def _update_gMax(self, dt):
    gf = self.Fg * self.a**3 * self.b * self.c
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

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

  def plot(self):
    _plt(self._fa, self._fb, self._fc)


class HH_aAbc:
  """
  gating: (a**A) * b * c
  """
  def __init__(self, fa, A, fb, fc, gMax, ER):
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
    self.gMax = gMax
    self.Fg   = gMax
    self.ER   = ER
    self.Rg   = 0.01
    self.gM   = None
    self.Pools = []
    self.Tag  = 'HH_aAbc'
    self.ptr = None

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return (self.a**self.A) * self.b * self.c

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    tc, fc = self._fc(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    self.c = exp(-dt/tc) * (self.c-fc) + fc

    if len(self.Pools) > 0:    #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def _update_gMax(self, dt):
    gf = self.Fg * self.a**3 * self.b * self.c
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

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

  def plot(self):
    _plt(self._fa, self._fb, self._fc)


class HH_aAbBc:
  """
  gating: (a**A) * (b**B) * c
  """
  def __init__(self, fa, A, fb, B, fc, gMax, ER):
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
    self.gMax = gMax
    self.Fg   = gMax
    self.ER   = ER
    self.Rg   = 0.01
    self.gM   = None
    self.Pools = []
    self.Tag  = 'HH_aAbBc'
    self.ptr = None

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return (self.a**self.A) * (self.b**self.B) * self.c

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    tc, fc = self._fc(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    self.c = exp(-dt/tc) * (self.c-fc) + fc

    if len(self.Pools) > 0:    #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def _update_gMax(self, dt):
    gf = self.Fg * self.a**3 * self.b * self.c
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

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

  def plot(self):
    _plt(self._fa, self._fb, self._fc)


class HH_aAbBcC:
  """
  gating: (a**A) * (b**B) * (c**C)
  """
  def __init__(self, fa, A, fb, B, fc, C, gMax, ER):
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
    self.gMax = gMax
    self.Fg   = gMax
    self.ER   = ER
    self.Rg   = 0.01
    self.gM   = None
    self.Pools = []
    self.Tag  = 'HH_a(A)b(B)c(C)'
    self.ptr = None

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return (self.a**self.A) * (self.b**self.B) * (self.c**self.C)

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    tc, fc = self._fc(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    self.c = exp(-dt/tc) * (self.c-fc) + fc

    if len(self.Pools) > 0:    #???
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def _update_gMax(self, dt):
    gf = self.Fg * self.a**3 * self.b * self.c
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

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

  def plot(self):
    _plt(self._fa, self._fb, self._fc)


class NaChannel_MHMJ:
  """
  MHMJ sodium channel
  """
  def __init__(self, fm, fh, fi, gMax, ER):
    self._fa = fm
    self._fb = fh
    self._fc = fi
    self.a  = 0.046
    self.b  = 0.638
    self.c  = 0.638
    self.a0 = 0.046
    self.b0 = 0.638
    self.c0 = 0.638
    self.gMax = gMax
    self.Fg   = gMax
    self.ER   = ER
    self.Rg   = 0.01
    self.gM   = None
    self.Tag  = 'MHMJ'

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return (self.a**3) * self.b * self.c

  def _update_gate(self, V, dt):
    ta, fa = self._fa(V)
    tb, fb = self._fb(V)
    tc, fc = self._fc(V)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    self.b = exp(-dt/tb) * (self.b-fb) + fb
    self.c = exp(-dt/tc) * (self.c-fc) + fc

  def _update_gMax(self, dt):
    gf = self.Fg * self.a**3 * self.b * self.c
    self.gMax = exp(-dt*self.Rg) * (self.gMax-gf) + gf

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

  def plot(self):
    _plt(self._fa, self._fb, self._fc)


def voltage_gated(gMax, ER, a = None, b = None, c = None):
  if a == None and b == None and c == None:
    return LeakChannel(gMax, ER)

  if b == None and c == None:
    Fa, A = a[0], a[1]
    if A == 1: return HH_a(Fa, gMax, ER)#
    if A >  1:
      if A == 4: # special case
        return KChannel(Fa, gMax, ER)
      else:
        return HH_aA(Fa, A, gMax, ER)#

  if c == None:
    Fa, A = a[0], a[1]
    Fb, B = b[0], b[1]
    if B == 1 and A == 1: return HH_ab(Fa, Fb, gMax, ER)#
    if B == 1 and A >  1:
      if A == 3: # special case
        return NaChannel(Fa, Fb, gMax, ER)
      else:
        return HH_aAb(Fa, A, Fb, gMax, ER)
    if B >  1 and A >  1: return HH_aAbB(Fa, A, Fb, B, gMax, ER)

  Fa, A = a[0], a[1]
  Fb, B = b[0], b[1]
  Fc, C = c[0], c[1]
  if C == 1 and B == 1 and A == 1: return HH_abc(Fa,Fb,Fc, gMax, ER)
  if C == 1 and B == 1 and A >  1: return HH_aAbc(Fa,A,Fb,Fc, gMax, ER)
  if C == 1 and B  > 1 and A >  1: return HH_aAbBc(Fa,A, Fb, B, Fc, gMax, ER)
  if C >  1 and B  > 1 and A >  1: return HH_aAbBcC(Fa, A, Fb, B, Fc, C, gMax, ER)


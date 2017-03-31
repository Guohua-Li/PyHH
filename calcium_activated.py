from math import exp
from pyhh import CaPool, _plt

def meshgrid(X,Y):
  M, N = len(X), len(Y)
  return [X for i in range(N)], [[Y[i]]*M for i in range(N)]

class HHC_a:
  """
  gating: a
  """
  def __init__(self, f, gMax, ER):
    self._f = f
    self.a  = 0.046
    self.a0 = 0.046
    self.gMax = gMax
    self.ER = ER
    self.Pools = []
    self.Tag = 'HHC_a'
    self.ptr = None

  def couple(self,p):
    self.pool_pointer = p
    p.Target = self

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return self.a

  def _update_gate(self, V, dt):
    if self.pool_pointer.Conc == 0: C = 0.0001
    else: C = self.pool_pointer.Conc
    tau, inf = self._f(V, C)
    self.a = exp(-dt/tau) * (self.a-inf) + inf
    if len(self.Pools) > 0:
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


class HHC_aA:
  """
  gating: (a**A)
  """
  def __init__(self, fa, A, gMax, ER):
    self._fa = fa
    self.A = A
    self.a0 = 0.046
    self.a  = 0.046
    self.gMax = gMax
    self.ER = ER
    self.Pools = []
    self.Tag = 'HHC_a(A)'
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
    return (self.a**self.A)

  def _update_gate(self, V, dt):
    if self.pool_pointer.Conc == 0: C = 0.0001
    else: C = self.pool_pointer.Conc
    ta, fa = self._fa(V,C)
    self.a = exp(-dt/ta) * (self.a-fa) + fa
    if len(self.Pools) > 0:
      J = self.gMax * self._ProbOpen() * (V - self.ER)
      for p in self.Pools:
        p.update(J, dt)

  def plot(self):
    _plt(self._fa)


class HHC_aAb:
  """
  gating: (a**A) * b
  """
  def __init__(self, fa, A, fb, gMax, ER):
    self._fa = fa
    self._fb = fb
    self.A = A
    self.a0 = 0.046
    self.b0 = 0.638
    self.a  = 0.046
    self.b  = 0.638
    self.gMax = gMax
    self.ER = ER
    self.Pools = []
    self.Tag = 'HHC_a(A)b'
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
    if self.pool_pointer.Conc == 0: C = 0.0001
    else: C = self.pool_pointer.Conc
    ta, fa = self._fa(V,C)
    tb, fb = self._fb(V,V)
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


class HHC_ab:
  """
  gating: a * b
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
    self.Tag = 'HHC_ab'
    self.ptr = None

  def couple(self,p):
    self.pool_pointer = p
    p.Target = self

  def add_CaMD(self,f):
    if self.ptr == None:
      raise Exception('Add CaMD after the channel is added to a compartment')
    d = CaPool(f)
    self.Pools.append(d)
    return d

  def _ProbOpen(self):
    return self.a * self.b

  def _update_gate(self, V, dt):
    if self.pool_pointer.Conc == 0: C = 0.0001
    else: C = self.pool_pointer.Conc
    ta, fa = self._fa(V,C)
    tb, fb = self._fb(V,C)
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

class HHC_aAbB:
  """
  (a**N) * (b**M)
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
    return (self.a**self.A) * (self.b**self.B)

  def _update_gate(self, V, dt):
    if self.pool_pointer.Conc == 0: C = 0.0001
    else: C = self.pool_pointer.Conc
    ta, fa = self._fa(V,C)
    tb, fb = self._fb(V,C)
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



class HHC_abc:
  """
  gating: a * b * c
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
    self.Tag  = 'HHC_abc'
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
    return self.a * self.b * self.c

  def _update_gate(self, V, dt):
    if self.pool_pointer.Conc == 0: C = 0.0001
    else: C = self.pool_pointer.Conc
    ta, fa = self._fa(V,C)
    tb, fb = self._fb(V,C)
    tc, fc = self._fc(V,C)
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

class HHC_aAbc:
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
    self.Tag  = 'HHC_aAbc'
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
    return (self.a**self.A) * self.b * self.c

  def _update_gate(self, V, dt):
    if self.pool_pointer.Conc == 0: C = 0.0001
    else: C = self.pool_pointer.Conc
    ta, fa = self._fa(V,C)
    tb, fb = self._fb(V,C)
    tc, fc = self._fc(V,C)
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




class HHC_aAbBc:
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
    self.Tag  = 'HH_a(A)b(B)c'
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
    return (self.a**self.A) * (self.b**self.B) * self.c

  def _update_gate(self, V, dt):
    if self.pool_pointer.Conc == 0: C = 0.0001
    else: C = self.pool_pointer.Conc
    ta, fa = self._fa(V,C)
    tb, fb = self._fb(V,C)
    tc, fc = self._fc(V,C)
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


class HHC_aAbBcC:
  """
  (a**A) * (b**B) * (c**C)
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
    return (self.a**self.A) * (self.b**self.B) * (self.c**self.C)

  def _update_gate(self, V, dt):
    if self.pool_pointer.Conc == 0: C = 0.0001
    else: C = self.pool_pointer.Conc
    ta, fa = self._fa(V,C)
    tb, fb = self._fb(V,C)
    tc, fc = self._fc(V,C)
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

def calcium_activated(gMax, ER, a, b = None, c = None):
  if b == None and c == None:
    Fa, A = a[0], a[1]
    if A == 1: return HHC_a(Fa, gMax, ER)#
    if A >  1: return HHC_aA(Fa, A, gMax, ER)#

  if c == None:
    Fa, A = a[0], a[1]
    Fb, B = b[0], b[1]
    if B == 1 and A == 1: return HHC_ab(Fa, Fb, gMax, ER)#
    if B == 1 and A >  1: return HHC_aAb(Fa, A, Fb, gMax, ER)
    if B >  1 and A >  1: return HHC_aAbB(Fa, A, Fb, B, gMax, ER)

  Fa, A = a[0], a[1]
  Fb, B = b[0], b[1]
  Fc, C = c[0], c[1]
  if C == 1 and B == 1 and A == 1: return HHC_abc(Fa,Fb)
  if C == 1 and B == 1 and A >  1: return HHC_aAbc(Fa,Fb)
  if C == 1 and B  > 1 and A >  1: return HHC_aAbBc(Fa,Fb)
  if C >  1 and B  > 1 and A >  1: return HHC_aAbBcC(Fa,Fb)



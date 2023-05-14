from math import exp

class Ligand:
  """
  This class defines ligand objects.
  """
  def __init__(self, concentration = 0):
    self.Conc = concentration
    self.Clamper = None

class TransRate:
  def __init__(self, func, category='V'):
    self.f = func
    self.Rate = 0
    self.Compt = None
    self.Category = category

  def update(self):
    V = self.Compt.Vi
    self.Rate = self.f(V)

class LGIC:
  """
  LGICs based on Markov transition schemes.
  """
  is_LGIC = 1
  def __init__(self, transit, binding, gMax, ER):
    self.transit = transit # keep it for copying the channel
    self.binding = binding # keep it for copying the channel
    self.gMax = gMax
    self.ER   = ER
    self.Tag = 'LGIC'
    self.ptr = None

    self.States = list(transit.keys()) # 2.x and 3.y are different
    self.States.sort()
    self.nStates = len(self.States)
    self.do_parsing(transit, binding)
    self._get_rate_matrix()
    self.N = len(self.RM)
    self.Prob = [0]*len(self.RM) # initialize the state probability vector
    self.Prob[0] = 1.0
    self.Ligands = list(binding.keys())
    if len(self.Ligands) == 1: self.Ligand = self.Ligands[0]

  def _set_gMax(self, val):
    self.gMax = val

  def _set_ER(self, val):
    self.gMax = val

  def do_parsing(self, transit, binding):
    StateIndexDict = {}
    for k in range(self.nStates):
      s = self.States[k]
      StateIndexDict[s] = k

    self.Transit = {}
    for i in range(self.nStates):
      si = self.States[i]
      for j in range(self.nStates):
        sj = self.States[j]
        if sj in transit[si].keys():
          if i in self.Transit:
            self.Transit[i].update({ j:transit[si][sj] })
          else:
            self.Transit[i] = { j:transit[si][sj] }
        else:
          if i in self.Transit:
            self.Transit[i].update({ j:0 })
          else:
            self.Transit[i] = { j:0 }

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
    for k in range(self.nStates):
      if self.States[k][0] == 'O': self.OpenIndex.append(k)

  def _get_rate_matrix(self):
    N = self.nStates
    states = range(N)
    self.RM = [[1]*N for i in states] # like ones((N, N), float)
    for i in states:
      for j in states:
        if isinstance(self.Transit[i][j], TransRate):
          self.RM[i][j] = self.Transit[i][j].Rate
        else:
          self.RM[i][j] = self.Transit[i][j]

  def _ProbOpen(self):
    gate = 0.0
    for k in self.OpenIndex: gate += self.Prob[k]
    return gate

  def _update_gate(self, V, h):
    self._get_rate_matrix()
    N = self.nStates
    states = range(N)
    for lg, paths in self.GateDict.items():
      for s, t in paths:
        self.RM[s][t] = self.RM[s][t] * lg.Conc
    A = [sum(i) for i in self.RM]           # sum of each row
    m = [ [j*self.Prob[i] for j in self.RM[i]] for i in states ] # self.RM * self.XC
    B = list(map(sum,zip(*m))) # sum of each vol, faster
    for k in states:
      if A[k] == 0:
        self.Prob[k] = self.Prob[k] + h * B[k]
      else:
        C = B[k]/A[k]
        self.Prob[k] = (self.Prob[k] - C) * exp(-A[k]*h) + C
    S = sum(self.Prob)
    if S != 0:
      self.Prob = [x/S for x in self.Prob]

  def show(self):
    print ("State Probabilities:")
    for k in range(self.nStates):
      print ("%5.4f"%(self.Prob[k]))



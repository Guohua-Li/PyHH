# modified on May 19, 2023
# modified on May 21, 2023
from math import exp

from gating import Gating

class Ligand:
  """
  This class defines ligand objects.
  """
  def __init__(self, concentration = 0):
    self.Conc = concentration
    self.Clamper = None


class VariableRate:
  def __init__(self, func, category = 'V'):
    self.Rate = 0
    self.f = func
    self.Compt = None
    self.Category = category
    self.Tag = "VariableRate"

  def update_rate(self):
    V = self.Compt.Vi
    self.Rate = self.f(V)


class LGIC(Gating):
  """
  LGICs based on Markov transition schemes.
  """
  is_LGIC = 1
  def __init__(self, transit, binding, recording="none"):
    Gating.__init__(self, recording)
    self.transit = transit # keep it for copying the channel
    self.binding = binding # keep it for copying the channel
    self.Tag = 'LGIC'
    self.ptr = None

    self.States = list(transit.keys()) # 2.x and 3.y are different
    self.States.sort()
    self.nStates = len(self.States)
    self.Transit = self.do_parsing(transit, binding)  # states in number, basically a sparse matrix in dic form

    self.OpenStateNumbers = [ k for k in range(self.nStates) if self.States[k][0] == 'O' ]
    self.LigandGatedDict = self.get_bind_paths(binding)
    """
    example of self.LigandGatedDict
    {<lgic.Ligand object at 0x7f156a967ee0>: [(0, 1), (1, 2)]}
    where 0->1, 1->2 are the paths for ligand binding
    """

    self.Prob = [0] * self.nStates # initialize the state probability vector
    self.Prob[0] = 1.0

    self.Ligands = list(binding.keys())
    if len(self.Ligands) == 1: self.Ligand = self.Ligands[0]

  def get_bind_paths(self, binding):
    StateMapIndex = {key: i for i, key in enumerate(self.States)}
    # dict(map(reversed, enumerate(self.States)))

    m_dict = {}
    for lg in binding.keys():
      paths = []
      scheme = binding[lg] # scheme is a dict
      for start in scheme.keys():
        end = scheme[start]
        i = StateMapIndex[start]
        j = StateMapIndex[end]
        paths.append((i,j))
      m_dict[lg] = paths
    return m_dict

  def do_parsing(self, transit, binding):
    trans_dict = {} # states in number, basically a sparse matrix in dic form
    state_numbers = range(self.nStates)
    for i in state_numbers:
      state_i = self.States[i]
      for j in state_numbers:
        state_j = self.States[j]
        if state_j in transit[state_i].keys():
          if i in trans_dict:
            trans_dict[i].update({ j:transit[state_i][state_j] })
          else:
            trans_dict[i] = { j:transit[state_i][state_j] }
        else:
          if i in trans_dict:
            trans_dict[i].update({ j: 0 })
          else:
            trans_dict[i] = { j: 0 }
    return trans_dict

  def _update_gate(self, V, dt): # used by xp, this is the solution of the ...
    state_numbers = range(self.nStates)
    rate_mat = [[0]*self.nStates for i in state_numbers] # like zeros/ones((N, N), float), N rows, N cols
    for i in state_numbers:
      for j in state_numbers:
        if isinstance(self.Transit[i][j], VariableRate): # NMDAR will use this
          rate_mat[i][j] = self.Transit[i][j].Rate # At the start, Rate = 0
        else:
          rate_mat[i][j] = self.Transit[i][j]


    for lg, paths in self.LigandGatedDict.items():
      for s, t in paths:
        rate_mat[s][t] = rate_mat[s][t] * lg.Conc

    A = [sum(i) for i in rate_mat]           # sum of each row
    m = [ [j*self.Prob[i] for j in rate_mat[i]] for i in state_numbers ] # self.RM * self.XC
    B = list(map(sum,zip(*m))) # sum of each vol, faster
    for k in state_numbers:
      if A[k] == 0:
        self.Prob[k] = self.Prob[k] + dt * B[k]
      else:
        C = B[k]/A[k]
        self.Prob[k] = (self.Prob[k] - C) * exp(-A[k]*dt) + C
    S = sum(self.Prob)
    if S != 0:
      self.Prob = [x/S for x in self.Prob]

  def _ProbOpen(self, step): # only used by compartment
    p = 0.0
    for k in self.OpenStateNumbers: p += self.Prob[k]
    if self.Recording: self.trace[step] = p
    return p

  def show(self):
    print ("State Probabilities:")
    for k in range(self.nStates):
      print ("%5.4f"%(self.Prob[k]))


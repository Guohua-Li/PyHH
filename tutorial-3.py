from pyhh import *

N = 50   # number of compartment to be created
D = 1.5  # diameter
L = 50   # length

P = [ Compartment(D, L, [NaC, KDR, gL]) for i in range(N) ]

for i in range(N-1):
  P[i].connect(P[i+1])

P[0].add_iclamper()
P[0].iClamper.Waveform = Rect(delay=1.5, width=1., amplitude=1.65)

xp = Experiment(P)
xp.run(30, 0.002)
xp.plot()


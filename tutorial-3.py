from pyhh import *


N = 50   # number of compartment to be created
D = 1.5  # diameter
L = 50   # length

P = [ Compartment(D, L) for i in range(N) ]

for i in range(N-1):
  P[i].connect(P[i+1])

for cpm in P:
  cpm.add_channels(NaC, KDR, gL)

P[0].add_iclamper()
P[0].iClamper.Waveform = Rect(delay=1.5, width=1., amplitude=1.65)

xp = Experiment(P)
xp.run(30, 0.002)
xp.plot()

"""
import matplotlib.pyplot as plt
fig, ax1 = plt.subplots()
for cpm in P:
  ax1.plot(xp.T, cpm.Vm, linewidth=2.0)

ax2 = ax1.twinx()
ax2.plot(xp.T, P[0].iClamper.Command, 'r-')
"""


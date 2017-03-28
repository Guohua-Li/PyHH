from pyhh import *
import pylab as plt

N = 50
D = 1.5
L = 50

P1 = []
for i in range(N):
  P1.append(Compartment(D, L))

#P1 = [ Compartment(diameter = 1.5, length = L) for i in range(1,N+1) ]

for cpm in P1:
  cpm.add_channels([NaC, KDR, gL])

for i in range(N-1):
  P1[i].connect(P1[i+1])

xp = Experiment(P1)

clp = IClamper()
clp.Waveform = Rect(delay=1.5, width=1.5, amplitude=1.4) # try amplitude = 1.3, 1.4, 1.5
#clp.Waveform = Alpha(delay=1, tau=0.5, amplitude=3.5)

clp.connect(P1[0])

xp.run(30, 0.002)

plt.figure()
plt.subplot(2,1,1)
for cpm in P1:
  plt.plot(xp.T, cpm.Vm, linewidth=2.0)
plt.ylim([-80,30])
plt.ylabel('V (mV)')

plt.subplot(2,1,2)
plt.plot(xp.T, clp.Command, linewidth=2.0)
plt.ylim([-0.2,1.5])
plt.xlabel('time (ms)')
plt.ylabel('current (pA/um2)')
plt.show()



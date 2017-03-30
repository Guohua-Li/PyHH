from pyhh import *
import pylab as plt



NaC = NaChannel(Fm, Fh, gMax=0.6, ER=50)    # the unit for gMax is nS/um2
KDR = KChannel(Fn, gMax=0.18, ER=-90)
gL  = LeakChannel(gMax=0.03, ER=-60)

NaL = NaChannel(Fm, Fh, gMax=0.2,ER=50)     # low density
NaH = NaChannel(Fm, Fh, 1.8, 50)            # high density

N = 30
D = 1.5
L = 50
Soma = Compartment(50,length=None)
#Soma.add_channels([NaC, KDR, gL])

P1 = [ Compartment(D, L) for i in range(N) ]

for cpm in P1:
  cpm.add_channels([NaC, KDR, gL])


P2 = [Compartment(D, L) for i in range(N)]
for cpm in P2:
  cpm.add_channels([NaL, KDR, gL])

P3 = [Compartment(D, L) for i in range(N)]
for cpm in P3:
  cpm.add_channels([NaH, KDR, gL])

for i in range(N-1):
  P1[i].connect(P1[i+1])
  P2[i].connect(P2[i+1])
  P3[i].connect(P3[i+1])

Soma.connect(P1[0])
Soma.connect(P2[0])
Soma.connect(P3[0])

xp = Experiment([Soma]+P1+P2+P3)
xp.run(5, 0.01)
xp.run(15, 0.02)
xp.Clock = 0
i_clp = IClamper()
i_clp.Waveform = Rect(delay=1, width=1.5, amplitude=3.0)
i_clp.connect(Soma)

xp.run(20, 0.005)
#xp.endpoint_info()

plt.figure()
plt.subplot(5,1,1)
plt.plot(xp.T, Soma.Vm, linewidth=2.0)
plt.ylim([-80,50])

plt.subplot(5,1,2)
for cpm in P1:
  plt.plot(xp.T, cpm.Vm, linewidth=2.0)
plt.ylim([-80,50])
plt.ylabel('V (mV)')

plt.subplot(5,1,3)
for cpm in P2:
  plt.plot(xp.T, cpm.Vm, linewidth=2.0)
plt.ylim([-80,50])
plt.ylabel('V (mV)')

plt.subplot(5,1,4)
for cpm in P3:
  plt.plot(xp.T, cpm.Vm, linewidth=2.0)
plt.ylim([-80,50])
plt.ylabel('V (mV)')

plt.subplot(5,1,5)
plt.plot(xp.T, i_clp.Command, linewidth=2.0)
plt.ylim([-0.5,4])
plt.ylabel('pA/um2')
plt.xlabel('time (ms)')
plt.show()


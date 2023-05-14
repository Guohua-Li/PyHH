from compartment import Compartment
from channels import Fm, Fh, NaChannel, Fn, KChannel, LeakChannel
from lgic import Ligand, LGIC
from clampers import CClamper, VClamper, IClamper, Rect
from experiment import Experiment
import matplotlib.pyplot as plt


NaC = NaChannel(Fm, Fh, gMax=0.6, ER=50)    # the unit for gMax is nS/um2
KDR = KChannel(Fn, gMax=0.18, ER=-90)
gL  = LeakChannel(gMax=0.03, ER=-60)

NaL = NaChannel(Fm, Fh, gMax=0.2, ER=50)     # low density
NaH = NaChannel(Fm, Fh, 1.0, 50)             # high density

N = 30

Soma = Compartment(50, length=None, channel_list=[NaH, KDR, gL])
P1 = [Compartment(2.0, 50.0, [NaH, KDR, gL]) for i in range(N)] # diameter = 2.0
P2 = [Compartment(1.5, 50.0, [NaH, KDR, gL]) for i in range(N)] # diameter = 1.5
P3 = [Compartment(1.0, 50.0, [NaH, KDR, gL]) for i in range(N)] # diameter = 1.0

for i in range(1, N):
  P1[i].attached_to(P1[i-1])
  P2[i].attached_to(P2[i-1])
  P3[i].attached_to(P3[i-1])

P1[0].attached_to(Soma)
P2[0].attached_to(Soma)
P3[0].attached_to(Soma)

xp = Experiment([Soma]+P1+P2+P3)
# relax to neuron
xp.run(5, 0.01)
xp.run(15, 0.02)


xp.Clock = 0
Soma.add_iclamper()
Soma.iClamper.Waveform = Rect(delay=1, width=1.5, amplitude=3.0)
xp.run(20, 0.005)

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
plt.plot(xp.T, Soma.iClamper.Command, linewidth=2.0)
plt.ylim([-0.5,4])
plt.ylabel('pA/um2')
plt.xlabel('time (ms)')
plt.show()

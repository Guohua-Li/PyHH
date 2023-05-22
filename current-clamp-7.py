from compartment import Compartment
from gating import NaT, KDR, gL
from clampers import Rect
from experiment import Experiment

import matplotlib.pyplot as plt


channels = {
    NaT: [0.8,  50],
    KDR: [0.18,-90],
    gL:  [0.03,-60]
}

num = 30

soma = Compartment(50, length=None, channels=channels)
P1 = [Compartment(2.0, 50.0, channels) for i in range(num)] # diameter = 2.0
P2 = [Compartment(1.5, 50.0, channels) for i in range(num)] # diameter = 1.5
P3 = [Compartment(1.0, 50.0, channels) for i in range(num)] # diameter = 1.0

for i in range(1, num):
  P1[i].attached_to(P1[i-1])
  P2[i].attached_to(P2[i-1])
  P3[i].attached_to(P3[i-1])

P1[0].attached_to(soma)
P2[0].attached_to(soma)
P3[0].attached_to(soma)

xp = Experiment([soma]+P1+P2+P3)
# relax to neuron
xp.run(5, 0.01)
xp.run(15, 0.02)


xp.Clock = 0
soma.add_iclamper()
soma.iClamper.Waveform = Rect(delay=1, width=1.5, amplitude=3.0)
xp.run(20, 0.005)

plt.figure()
plt.subplot(5,1,1)
plt.plot(xp.T, soma.Vm, linewidth=2.0)
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
plt.plot(xp.T, soma.iClamper.Command, linewidth=2.0)
plt.ylim([-0.5,4])
plt.ylabel('pA/um2')
plt.xlabel('time (ms)')
plt.show()

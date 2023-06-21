from compartment import Compartment
from experiment import Experiment
from gating import NaT, KDR, gL
from ploting import plot_Vm



channels = { NaT: [0.8,  50], KDR: [0.18,-90], gL:  [0.03,-60] }

num = 15

soma = Compartment(40, length=None, channels=channels)
P1 = [Compartment(3.0, 20.0, channels) for i in range(num)] # diameter = 2.0
P2 = [Compartment(2.0, 20.0, channels) for i in range(num)] # diameter = 1.5
P3 = [Compartment(1.0, 20.0, channels) for i in range(num)] # diameter = 1.0

for i in range(1, num):
  P1[i].attached_to(P1[i-1])
  P2[i].attached_to(P2[i-1])
  P3[i].attached_to(P3[i-1])

P1[0].attached_to(soma)
P2[0].attached_to(soma)
P3[0].attached_to(soma)

xp = Experiment([soma]+P1+P2+P3)
xp.run(15, 0.01)

soma.add_iclamper()
soma.iClamper.set_waveform('rect', delay = 2, width=1.5, amplitude=9)
xp.run(15, 0.001)

plot_Vm([soma], P1, P2, P3)


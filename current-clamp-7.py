from compartment import Compartment
from experiment import Experiment
from gating import NaT, KDR, gL
from ploting import plot_Vm



channels_medium = { NaT: [0.6,  50], KDR: [0.18,-90], gL:  [0.03,-60] }
channels_high   = { NaT: [1.8,  50], KDR: [0.18,-90], gL:  [0.03,-60] }
channels_low    = { NaT: [0.2,  50], KDR: [0.18,-90], gL:  [0.03,-60] }

num = 15
diam = 2.5
leng = 20

soma = Compartment(40, length=None)

P1 = [Compartment(diam, leng, channels_low) for i in range(num)] # For the sodium channel gMax = 0.6
P2 = [Compartment(diam, leng, channels_medium) for i in range(num)] # For the sodium channel gMax = 0.2
P3 = [Compartment(diam, leng, channels_high) for i in range(num)] # For the sodium channel gMax = 1.8

P1[0].attached_to(soma)
P2[0].attached_to(soma)
P3[0].attached_to(soma)

for i in range(1,num):
    P1[i].attached_to(P1[i-1])
    P2[i].attached_to(P2[i-1])
    P3[i].attached_to(P3[i-1])

xp = Experiment([soma]+P1+P2+P3)
xp.run(15, 0.003, verbose=False)

xp.Clock = 0
soma.add_iclamper()
soma.iClamper.set_waveform('rect', delay = 1, width=2, amplitude=8)
xp.run(15)

plot_Vm([soma], P1, P2, P3, ylim=[-90,50])

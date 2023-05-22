from compartment import Compartment
from gating import NaT, KDR, gL
from clampers import Alpha
from experiment import Experiment

import matplotlib.pyplot as plt

def my_plot(t, v1, v2, cmd):
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(t, v1)
    plt.plot(t, v2, linewidth=2.0)
    plt.ylim([-80,30])
    plt.ylabel('V (mV)')
    plt.subplot(2,1,2)
    plt.plot(t, cmd, linewidth=2.0)
    plt.ylim([-1,2])
    plt.xlabel('time (ms)')
    plt.ylabel('current (pA/um2)')
    plt.show()

channels = {
    NaT: [0.6,  50],
    KDR: [0.18,-90],
    gL:  [0.03,-60]
}

soma = Compartment(50, None, channels)
dend = Compartment(1.5, 100, channels, parent=soma)

clamper = soma.add_iclamper()
clamper.Waveform = Alpha(delay=2, tau=0.5, amplitude=0.9)

preparations = [soma,dend]

xp = Experiment(preparations)
xp.run(20)
my_plot(xp.T, soma.Vm, dend.Vm, soma.iClamper.Command)

xp.Clock = 0
xp.run(10)

my_plot(xp.T, soma.Vm, dend.Vm, soma.iClamper.Command)


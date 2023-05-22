from compartment import Compartment
from gating import NaT, KDR, gL
from clampers import Alpha # stimulation waveform
from experiment import Experiment

channels = {
    NaT: [0.6,  50],
    KDR: [0.18,-90],
    gL:  [0.03,-60]
}

soma = Compartment(50, None, channels)
dend = Compartment(1.5, 100, channels, parent=soma) # soma can not be a child.

clamper = soma.add_iclamper()
clamper.Waveform = Alpha(delay=2, tau=0.5, amplitude=0.9)

preparations = [soma, dend]

xp = Experiment(preparations)
xp.run(10)
xp.plot_Vm() # xp stores the membrane potentials of all compartments

xp.run(10)
xp.plot_Vm()

xp.Clock = 0
xp.run(10)
xp.plot_Vm()

"""
# make your own plot:
import matplotlib.pyplot as plt
plt.figure()
plt.subplot(2,1,1)
plt.plot(xp.T, soma.Vm)
plt.plot(xp.T, dend.Vm, linewidth=2.0)
plt.ylim([-80,30])
plt.ylabel('V (mV)')
plt.subplot(2,1,2)
plt.plot(xp.T, soma.iClamper.Command, linewidth=2.0)
plt.ylim([-1,2])
plt.xlabel('time (ms)')
plt.ylabel('current (pA/um2)')
plt.show()
"""

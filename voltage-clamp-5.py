from compartment import Compartment
from gating import NaT, KDR, gL
from clampers import Rect
from experiment import Experiment

channels = {
    NaT: [0.6,  50],
    KDR: [0.18,-90],
    gL:  [0.03,-60]
}

soma = Compartment(diameter = 50., length=None,  channels=channels)
dend = Compartment(diameter = 1.5, length = 100, channels=channels, parent=soma)

clamper = soma.add_vclamper(-60)
clamper.Waveform = Rect(delay=1, width=5, amplitude=40)

xp = Experiment([soma, dend])
xp.run(10)

clamper.plot()


#clamper.plot(show_Jm=1,show_Jn=1,show_Jp=1)
#clamper.plot(show_Jm=1,show_Jn=0) # show Jm only

"""
import matplotlib.pyplot as plt
plt.figure()
plt.subplot(2,1,1)
plt.plot(xp.T, soma.vClamper.Jm, linewidth=2.0)
plt.plot(xp.T, soma.vClamper.Jn, linewidth=2.0)
plt.plot(xp.T, soma.vClamper.Jc, linewidth=2.0)
plt.ylabel('Jm, Jc, Jn')
plt.ylim([-5,5])
plt.subplot(2,1,2)
plt.plot(xp.T, soma.Vm, linewidth=2.0)
plt.plot(xp.T, dend.Vm, linewidth=2.0)
plt.xlabel('time (ms)')
plt.ylabel('V (mV)')
plt.show()
Jp = soma.vClamper.calc_Jp()
plt.figure()
plt.plot(xp.T, Jp, linewidth=2.0)
plt.plot(xp.T, soma.vClamper.Jm)
plt.show()
plt.figure()
Jp = [i+j for i,j in zip(soma.vClamper.Jm, soma.vClamper.Jn)]
plt.plot(xp.T, Jp, linewidth=2.0)
plt.plot(xp.T, soma.vClamper.Jm)
plt.show()
"""

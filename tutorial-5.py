from pyhh import *
import pylab as plt

soma = Compartment(diameter = 50, length=None)
soma.add_channels([NaC, KDR, gL])
dend = Compartment(diameter = 1.5, length = 100)
dend.add_channels([NaC, KDR, gL])
soma.connect(dend)

soma.add_vclamper()
soma.vClamper.Waveform = Rect(delay=1, width=5, amplitude=40)

xp = Experiment([soma,dend])
xp.run(t=10,dt=0.005)

plt.figure()
plt.subplot(3,1,1)
plt.plot(xp.T, soma.vClamper.Jp, linewidth=2.0)
plt.ylabel('Jp (pA/um2)')
plt.ylim([-5,5])
plt.subplot(3,1,2)
plt.plot(xp.T, soma.vClamper.Jm, linewidth=2.0)
plt.plot(xp.T, soma.vClamper.Jn, linewidth=2.0)
plt.plot(xp.T, soma.vClamper.Jc, linewidth=2.0)
plt.ylabel('Jm, Jc, Jn')
plt.ylim([-5,5])

plt.subplot(3,1,3)
plt.plot(xp.T, soma.Vm, linewidth=2.0)
plt.plot(xp.T, dend.Vm, linewidth=2.0)
plt.xlabel('time (ms)')
plt.ylabel('V (mV)')
plt.show()


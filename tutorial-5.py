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



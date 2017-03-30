from pyhh import *
import pylab as plt

soma = Compartment(diameter = 50, length=None)
soma.add_channels([NaC, KDR, gL])

clp = VClamper()
clp.Waveform = Rect(delay=1, width=5, amplitude=40)
clp.connect(soma)

xp = Experiment(soma)
xp.run(t=10,dt=0.005)
clp.save('c.txt')

plt.figure()
plt.subplot(3,1,1)
plt.plot(xp.T, clp.Jp, linewidth=2.0)
plt.ylabel('Jp (pA/um2)')
plt.ylim([-5,5])
plt.subplot(3,1,2)
plt.plot(xp.T, clp.Jm, linewidth=2.0)
plt.plot(xp.T, clp.Jn, linewidth=2.0)
plt.plot(xp.T, clp.Jc, linewidth=2.0)
plt.ylabel('Jc, Jm')
plt.ylim([-5,5])
plt.subplot(3,1,3)
plt.plot(xp.T, soma.Vm, linewidth=2.0)
plt.xlabel('time (ms)')
plt.ylabel('V (mV)')
plt.ylim([-65,-25])

plt.show()


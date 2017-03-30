from pyhh import *
import pylab as plt

soma = Compartment(50,length=50)
dend = Compartment(1.5,length=100)

soma.add_channels([NaC, KDR, gL])
dend.add_channels([NaC, KDR, gL])

soma.connect(dend)

soma.add_iclamper()
soma.iClamper.Waveform = Rect(delay=2, width=1.5, amplitude=1.0)

xp = Experiment([soma,dend])


xp.run(10, 0.002)
xp.save('example.txt')

plt.figure()
plt.subplot(2,1,1)
plt.plot(xp.T, soma.Vm, linewidth=2.0)
plt.plot(xp.T, dend.Vm, linewidth=2.0)
plt.ylim([-80,30])
plt.ylabel('V (mV)')

plt.subplot(2,1,2)
plt.plot(xp.T, soma.iClamper.Command, linewidth=2.0)
plt.ylim([-1,2])
plt.xlabel('time (ms)')
plt.ylabel('current (pA/um2)')

plt.show()



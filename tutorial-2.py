from pyhh import *

soma = Compartment(50,length=None)
dend = Compartment(1.5,length=100)

soma.add_channels(NaC, KDR, gL)
dend.add_channels(NaC, KDR, gL)

soma.connect(dend) # or dend.attached_to(soma)

soma.add_iclamper()
soma.iClamper.Waveform = Rect(delay=2, width=0.5, amplitude=1.45)

xp = Experiment([soma,dend])
xp.run(10, 0.002)
xp.plot()

xp.run(10, 0.005)
xp.plot()

xp.Clock = 0
xp.run(10, 0.005)
xp.plot()

"""
# make your own plot:
import pylab as plt
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


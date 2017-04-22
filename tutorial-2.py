from pyhh import *

soma = Compartment(50, None, [NaC, KDR, gL])
dend = Compartment(1.5, 100, [NaC, KDR, gL])
dend.attached_to(soma)

soma.add_iclamper()
soma.iClamper.Waveform = Alpha(delay=2, tau=0.5, amplitude=0.9)

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


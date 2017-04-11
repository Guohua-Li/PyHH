from pyhh import *

soma = Compartment(diameter = 50, length=None)
soma.add_channels(NaC, KDR, gL)

clp = VClamper()
clp.Waveform = Rect(delay=1, width=5, amplitude=40)
clp.connect(soma)

xp = Experiment(soma)
xp.run(t=10,dt=0.005)

clp.plot()
clp.plot(show_Jm=1,show_Jn=0) # show Jm, not Jn
clp.plot(show_Jm=1,show_Jn=0,show_Jp=1) # show Jm, Jp, but not Jn

"""
# make your own plot:
import pylab as plt
plt.figure()
plt.subplot(2,1,1)
plt.plot(xp.T, clp.Jm, linewidth=3.0)
plt.plot(xp.T, clp.Jn, linewidth=2.0)
plt.plot(xp.T, clp.Jc, linewidth=1.0)
plt.ylabel('Jm, Jn, Jc')
plt.ylim([-5,5])
plt.subplot(2,1,2)
plt.plot(xp.T, soma.Vm, linewidth=2.0)
plt.xlabel('time (ms)')
plt.ylabel('V (mV)')
plt.ylim([-65,-15])
plt.show()

Jp = soma.vClamper.calc_Jp()

plt.figure()
plt.plot(xp.T, Jp, linewidth=2.0)
plt.ylabel('Jp (pA/um2)')
plt.ylim([-5,5])
plt.show()
"""

from pyhh import *
import pylab as plt

cpm = Compartment(1.5, 100)
cpm.add_channels([NaC, KDR, gL])

clp = IClamper()
clp.Waveform = Rect(delay=1, width=1, amplitude=0.7)
clp.connect(cpm)

xp = Experiment([cpm])
xp.run(10, 0.005)

plt.figure()
plt.plot(xp.T, cpm.Vm, linewidth=2.0)
plt.ylim([-80,30])
plt.show()


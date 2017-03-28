from pyhh import *

cpm = Compartment(diameter = 1.5, length = 100)
cpm.add_channels([NaC, KDR, gL])

xp = Experiment(cpm)
xp.run(t=10,dt=0.005)
cpm.save('a.txt')

import pylab as plt
plt.plot(xp.T, cpm.Vm, linewidth=2.0)
plt.ylim([-80,50])
plt.show()


"""
We provide an alternative way of ploting the membrane potentials:
The plot_Vm function accepts compartment list, so items in the same list can be plotted  together
"""


from compartment import Compartment
from experiment import Experiment
from gating import NaT, KDR, gL

from ploting import plot_Vm # an alternative way of ploting the membrane potentials


channels = { NaT: [0.6,  50], KDR: [0.18,-90], gL:  [0.03,-60] }

soma = Compartment(diameter=40, length=None, channels=channels) 
compartments = [soma]

num = 30   # number of compartments to be created
diam = 1.5  # diameter of each compartment
leng = 50   # length of each compartment

for i in range(1,num):
    p = Compartment(diam, leng, channels=channels, parent=compartments[i-1])
    compartments.append(p)

"""
for i in range(num-1):
    p = Compartment(diam, leng, [NaT, KDR, gL], parent=compartments[-1])
    compartments.append(p)
"""

clamper = compartments[0].add_iclamper('rect', delay = 1.5, width=1.0, amplitude=1.9)

xp = Experiment(compartments)
xp.run(25)

#plot_Vm(compartments, ylim=[-80,30])
plot_Vm([soma], compartments[1:], ylim=[-80,30])

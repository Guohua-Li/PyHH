"""
This tutorial demostrate the how we simulate the behivior of a neuron with three compartments.

"""

from compartment import Compartment
from gating import NaT, KDR, gL
from experiment import Experiment
from ploting import plot_Vm

channels = { NaT: [0.6,  50], KDR: [0.18,-90], gL:  [0.03,-60] }

soma = Compartment(40, None, channels)
dend1 = Compartment(3, 60, channels, parent=soma) # soma can not be a child.
dend2 = Compartment(3, 30, channels, parent=soma) # soma can not be a child.

clamper = soma.add_iclamper('alpha', delay = 2, tau=0.5, amplitude=1.9)

xp = Experiment(soma, dend1, dend2)

xp.run(15)
#xp.plot_Vm() # xp stores the membrane potentials of all compartments

plot_Vm([soma], [dend1, dend2])

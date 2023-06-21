"""
This tutorial shows how the sodium (in green) and potassium (in purple) 
conductances changes during an action potential.

"""

from compartment import Compartment
from experiment import Experiment
from gating import NaT, KDR, gL
from ploting import plot_Vm


channels = { NaT: [0.6,  50], KDR: [0.18,-90], gL:  [0.03,-60] }

soma = Compartment(40, None, channels)
dend = Compartment(1.5, 100, channels, parent=soma) # soma can not be a child.

for chnn in soma.channel_list:
    print(chnn.Tag)

nav = soma.get_channel("NaT")
nav.set_recording(True)

kdr = soma.get_channel("KDR")
kdr.set_recording(True)

clamper = soma.add_iclamper("alpha", delay=2, tau=0.5, amplitude=0.9)

xp = Experiment(soma, dend)

xp.run(15)
#xp.plot_Vm() # xp stores the membrane potentials of all compartments

#nav.plot_conductance(xp.T)

plot_Vm([soma], g_grid=True)

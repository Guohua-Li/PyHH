"""
In this tutorial, we are going to inject a rectangle current pulse 
to see the changes in membrane potentials. For this purpose, 
we need to import the IClamper and the Rect classes.

"""

from compartment import Compartment
from experiment import Experiment
from gating import NaT, KDR, gL
from clampers import IClamper
from ploting import plot_recording




chann_dist = {
    NaT: [0.6,  50],
    KDR: [0.18,-90],
    gL:  [0.03,-60]
}

cell = Compartment(diameter=40, length=None, channels = chann_dist)

clamper = IClamper(cell) # first, we define a clamper
clamper.set_waveform('rect', delay = 2.5, width=1)

xp = Experiment(cell)


VmTraces = [] # store potential from each run
StTraces = [] # store stimulus from each run
amplitudes = [-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.52, 0.54, 0.56]

for amp in amplitudes:
    clamper.set_amplitude(amp)
    xp.Clock = 0
    xp.run(t = 20)
    VmTraces.append(cell.Vm)
    StTraces.append(clamper.Command)




plot_recording(xp.T, VmTraces, StTraces)

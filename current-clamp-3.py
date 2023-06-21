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




chann_dist = { NaT: [0.6,  50], KDR: [0.18,-90], gL:  [0.03,-60] }

cell = Compartment(diameter=40, length=None, channels = chann_dist)
clamper = IClamper(cell) # first, we define a clamper
clamper.set_waveform('rect', delay = 2.5, amplitude=0.45)

xp = Experiment(cell)


VmTraces = [] # store potential from each run
StTraces = [] # store stimulus from each run
widths = [0.2, 0.4, 0.6, 0.8,1, 1.2, 1.4, 1.6]

for width in widths:
    clamper.set_width(width)
    xp.run(t = 20)
    VmTraces.append(cell.Vm)
    StTraces.append(clamper.Command)

plot_recording(xp.T, VmTraces, StTraces)

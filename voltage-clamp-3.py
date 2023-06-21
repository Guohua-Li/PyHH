from compartment import Compartment
from gating import KDR
from clampers import VClamper
from experiment import Experiment

channels = {  KDR: [0.18,-90],  }

cell = Compartment(50, None, channels)

clamper = VClamper(cell, -60)

xp = Experiment(cell)
xp.run(t=30)

clamper.plot_current()

clamper.set_waveform('rect', delay = 2, width=20, amplitude=40)
xp.Clock = 0
xp.run(t=30)
clamper.plot_current()

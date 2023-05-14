from compartment import Compartment
from channels import KDR
from clampers import VClamper, Rect
from experiment import Experiment


cell = Compartment(50, None, [KDR])

clamper = VClamper(-60)
clamper.Waveform = Rect(delay=5, width=20, amplitude=0)
clamper.connect(cell)

xp = Experiment(cell)
xp.run(t=30)

clamper.plot()

clamper.Waveform = Rect(delay=2, width=20, amplitude=40)
xp.Clock = 0
xp.run(t=30)
clamper.plot()

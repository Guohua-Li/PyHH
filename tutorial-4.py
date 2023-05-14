from compartment import Compartment
from channels import NaC, KDR, gL
from clampers import Rect, VClamper
from experiment import Experiment

cell = Compartment(diameter=50, length=None, channel_list=[NaC, KDR, gL])

clamper = VClamper(baseline = -60)
clamper.Waveform = Rect(delay=5, width=20, amplitude=0)
clamper.connect(cell)

xp = Experiment(cell)
xp.run(30)
clamper.plot()

xp.Clock = 0
clamper.set_amplitude(40) # now, delay=5, width=20, amplitude=40
xp.run(30)
clamper.plot()

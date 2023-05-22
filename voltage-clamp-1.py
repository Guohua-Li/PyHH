from compartment import Compartment
from gating import NaT, KDR, gL
from clampers import Rect, VClamper
from experiment import Experiment


channels = {
    NaT: [0.6,  50],
    KDR: [0.18,-90],
    gL:  [0.03,-60]
}


cell = Compartment(diameter=50, length=None, channels=channels)

clamper = VClamper(baseline = -60)
clamper.Waveform = Rect(delay=5, width=20, amplitude=0)
clamper.clamp(cell)

xp = Experiment(cell)
xp.run(30)
clamper.plot() # volt clampers plot Jm, Jn, or Jp

xp.Clock = 0
clamper.set_amplitude(40) # now, delay=5, width=20, amplitude=40
xp.run(30)
clamper.plot()


from compartment import Compartment
from gating import KDR
from clampers import VClamper, Rect
from experiment import Experiment

channels = {
    KDR: [0.18,-90],
}

cell = Compartment(50, None, channels)

clamper = VClamper(-60)
clamper.Waveform = Rect(delay=5, width=20, amplitude=0)
clamper.clamp(cell)

xp = Experiment(cell)
xp.run(t=30)

clamper.plot()

clamper.Waveform = Rect(delay=2, width=20, amplitude=40)
xp.Clock = 0
xp.run(t=30)
clamper.plot()
#clamper.plot(show_Jm=1,show_Jn=0) # show Jm, not Jn
#clamper.plot(show_Jm=1,show_Jn=0,show_Jp=1) # show Jm, Jp, but not Jn


from compartment import Compartment
from gating import NaT
from clampers import VClamper, Rect
from experiment import Experiment

channels = {
    NaT: [0.6, 50],
}

soma = Compartment(50, None, channels)

# another way of adding a clamper:
clamper = soma.add_vclamper(-60)
clamper.Waveform = Rect(delay=5, width=20, amplitude=0)

xp = Experiment(soma)
xp.run(30)
clamper.plot()

clamper.Waveform = Rect(delay=2, width=20, amplitude=40) # amplitude relative to the baseline
xp.Clock = 0
xp.run(t=30, dt=0.002) # we set the step to be 0.002
clamper.plot()

#clamper.plot(show_Jm=1,show_Jn=0) # show Jm, not Jn, need more explaination
#clamper.plot(show_Jm=1,show_Jn=0,show_Jp=1) # show Jm, Jp, but not Jn

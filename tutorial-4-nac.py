from compartment import *
from channels import NaC
from clampers import VClamper, Rect
from experiment import Experiment


soma = Compartment(50, None, channel_list = [NaC])

clamper = VClamper(-60) # holding at -60 mV baseline
clamper.Waveform = Rect(delay=5, width=20, amplitude=0)
clamper.connect(soma)
# instead, you can add a voltage clamper like:
#clamper = soma.add_vclamper(-60)
#soma.vClamper.Waveform = Rect(delay=5, width=20, amplitude=0)

xp = Experiment(soma)
xp.run(30) # default integration step is 0.005
clamper.plot()

clamper.Waveform = Rect(delay=2, width=20, amplitude=40) # amplitude relative to the baseline
xp.Clock = 0
xp.run(t=30, dt=0.002) # we set the step to be 0.002
clamper.plot()

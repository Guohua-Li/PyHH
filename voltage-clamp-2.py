from compartment import Compartment
from gating import NaT
from clampers import VClamper
from experiment import Experiment

channels = {
    NaT: [0.6, 50],
}

soma = Compartment(50, None, channels)


clamper = soma.add_vclamper(-60)


xp = Experiment(soma)
xp.run(30)
clamper.plot_current()

clamper.set_waveform('rect', delay = 2, width=20, amplitude=40)
xp.run(t=30, dt=0.002) # we set the step to be 0.002
clamper.plot_current()

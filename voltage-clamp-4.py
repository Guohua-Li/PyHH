from compartment import Compartment
from gating import NaT, KDR, gL
from experiment import Experiment


channels = { NaT: [0.6,  50], KDR: [0.18,-90], gL:  [0.03,-60] }

soma = Compartment(diameter = 50, length=None) # if length is set to be None, the compartment is treated as a sphere
soma.add_channels(channels)
dend = Compartment(diameter = 1.5, length = 100)
dend.add_channels(channels)
dend.attached_to(soma) # or soma.connect(dend)

clamper = soma.add_vclamper(-60)
clamper.set_waveform('rect', delay = 1, width=4, amplitude=40)

xp = Experiment([soma, dend])
xp.run(10)

soma.vClamper.plot_current()

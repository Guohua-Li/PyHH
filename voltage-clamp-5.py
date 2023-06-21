from compartment import Compartment
from gating import NaT, KDR, gL
from experiment import Experiment


channels = { NaT: [0.6,  50], KDR: [0.18,-90], gL:  [0.03,-60] }

soma = Compartment(diameter = 50., length=None,  channels=channels)
dend = Compartment(diameter = 1.5, length = 100, channels=channels, parent=soma)

clamper = soma.add_vclamper(-60)
clamper.set_waveform('rect', delay = 1, width=5, amplitude=40)

xp = Experiment([soma, dend])
xp.run(10)

clamper.plot_current(show_Jp=True)

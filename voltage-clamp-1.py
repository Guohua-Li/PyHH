from compartment import Compartment
from gating import NaT, KDR, gL
from clampers import VClamper
from experiment import Experiment
from ploting import plot_Vm


channels = {
    NaT: [0.6,  50],
    KDR: [0.18,-90],
    gL:  [0.03,-60]
}


cell = Compartment(diameter=50, length=None, channels=channels)

clamper = VClamper(cell, baseline = -60)
#clamper.clamp(cell)

nav = cell.get_channel("NaT")
nav.set_recording(True)

xp = Experiment(cell)
xp.run(30)
clamper.plot_current() # volt clampers plot Jm, Jn, or Jp

xp.Clock = 0
clamper.set_waveform('rect', delay = 5, width=20, amplitude=40)

xp.run(30)
clamper.plot_current()

plot_Vm([cell])

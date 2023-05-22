from compartment import Compartment
from gating import NaT, KDR, gL
from clampers import IClamper, Rect # IClamoer: Current clamper,  Rect: rectangular waveform for stimulation
from experiment import Experiment


cell = Compartment(diameter=40, length=None)

channels = {
    NaT: [0.6,  50],
    KDR: [0.18,-90],
    gL:  [0.03,-60]
}

cell.add_channels(channels)

clamper = IClamper()
clamper.Waveform = Rect(delay=1, width=1, amplitude=0.55)
clamper.clamp(cell)

xp = Experiment(cell)
xp.run(t = 10)
cell.plot_Vm(vlimit=[-80, 20])

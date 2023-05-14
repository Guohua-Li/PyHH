from compartment import Compartment
from channels import NaC, KDR, gL
from clampers import Rect
from experiment import Experiment


#cell = Compartment(diameter=1.5, length=100)
cell = Compartment(diameter=1.5, length=None)
cell.add_channels(NaC, KDR, gL)
cell.add_iclamper()
cell.iClamper.Waveform = Rect(delay=1, width=1, amplitude=0.55)

xp = Experiment(cell)
xp.run(t = 10) # , dt=0.005
cell.plot()

from compartment import Compartment
from channels import NaC, KDR, gL
from clampers import Rect
from experiment import Experiment


soma = Compartment(diameter = 50, length=None) # if length is set to be None, the compartment is treated as a sphere
soma.add_channels(NaC, KDR, gL)
dend = Compartment(diameter = 1.5, length = 100)
dend.add_channels(NaC, KDR, gL)
dend.attached_to(soma) # or soma.connect(dend)

soma.add_vclamper(-60)
soma.vClamper.Waveform = Rect(delay=1, width=5, amplitude=40)

xp = Experiment([soma, dend])
xp.run(10)

soma.vClamper.plot()

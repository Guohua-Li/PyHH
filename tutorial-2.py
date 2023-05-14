from compartment import Compartment
from channels import NaC, KDR, gL
from clampers import Alpha
from experiment import Experiment


soma = Compartment(50, None, [NaC, KDR, gL])
dend = Compartment(1.5, 100, [NaC, KDR, gL])
dend.attached_to(soma)

clamper = soma.add_iclamper()
clamper.Waveform = Alpha(delay=2, tau=0.5, amplitude=0.9)

xp = Experiment([soma,dend])
xp.run(10, 0.002)
xp.plot()

xp.run(10, 0.005)
xp.plot()

xp.Clock = 0
xp.run(10, 0.005)
xp.plot()

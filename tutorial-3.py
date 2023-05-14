from compartment import Compartment
from channels import NaC, KDR, gL
from clampers import Rect
from experiment import Experiment


N = 50   # number of compartments to be created
D = 1.5  # diameter of each compartment
L = 50   # length of each compartment

P = [ Compartment(D, L, [NaC, KDR, gL]) for i in range(N) ]

for i in range(1, N):
    P[i].attached_to(P[i-1])

clamper = P[0].add_iclamper()
clamper.Waveform = Rect(delay=1.5, width=1., amplitude=1.65)

xp = Experiment(P)
xp.run(30)
xp.plot()

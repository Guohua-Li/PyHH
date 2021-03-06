from pyhh import *

soma = Compartment(diameter=50, length=None, channel_list=[NaC, KDR, gL])

clp = VClamper()
clp.Waveform = Rect(delay=5, width=20, amplitude=0)
clp.connect(soma)

xp = Experiment(soma)
xp.run(t=30,dt=0.005)

clp.plot()

clp.set_amplitude(40)

xp.Clock = 0
xp.run(t=30,dt=0.005)
clp.plot()


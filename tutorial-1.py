from pyhh import *
cpm = Compartment(diameter=1.5, length=100)
cpm.add_channels(NaC, KDR, gL)
cpm.add_iclamper()
cpm.iClamper.Waveform = Rect(delay=1, width=1, amplitude=0.55)
xp = Experiment(cpm)
xp.run(t=10, dt=0.005)
cpm.plot()


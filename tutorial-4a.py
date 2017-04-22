from pyhh import *

soma = Compartment(50, None, channel_list = [NaC])

clp = VClamper()
clp.Waveform = Rect(delay=5, width=20, amplitude=0)
clp.connect(soma)

xp = Experiment(soma)
xp.run(30,0.005)

clp.plot()

clp.Waveform = Rect(delay=2, width=20, amplitude=40)
xp.Clock = 0
xp.run(t=30,dt=0.005)
clp.plot()
clp.plot(show_Jm=1,show_Jn=0) # show Jm, not Jn
clp.plot(show_Jm=1,show_Jn=0,show_Jp=1) # show Jm, Jp, but not Jn


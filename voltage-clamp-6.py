from compartment import Compartment
from gating import NaT
from clampers import VClamper
from experiment import Experiment
from ploting import plot_recording



channels = {
    NaT: [0.6, 50],
}

soma = Compartment(50, None, channels)

# another way of adding a clamper:
clamper = soma.add_vclamper(-60)
clamper.set_waveform('rect', delay=2, width=20)

xp = Experiment(soma)
xp.run(30)
#clamper.plot_current()

JiTraces = []
StTraces = []
amplitudes = [0.0, 10, 20, 30, 40, 50, 60]

for amp in amplitudes:
    clamper.set_amplitude(amp)
    xp.run(t=30, dt=0.001)
    JiTraces.append(clamper.J_ion)
    StTraces.append(clamper.Command)
    #clamper.plot_current()

#clamper.plot(show_Jm=1,show_Jn=0) # show Jm, not Jn, need more explaination
#clamper.plot(show_Jm=1,show_Jn=0,show_Jp=1) # show Jm, Jp, but not Jn

plot_recording(xp.T, JiTraces, StTraces, 'vclamp')

from compartment import Compartment
from lgic import LGIC, Ligand
from clampers import VClamper, CClamper
from experiment import Experiment
from ploting import plot_Vm



GABA = Ligand()

gabar_transit = {
                  'Cu': {'Cb': 5.0},
                  'Cb': {'Cu': 1.0, 'O': 0.8, 'D': 0.8},
                  'O':  {'Cb': 0.5},
                  'D':  {'Cb': 0.1},
                }

gabar_binding = {GABA:
                  {'Cu': 'Cb'}
                }


GABA_R = LGIC(gabar_transit, gabar_binding)

channels = {
    GABA_R: [0.01,  0], # about -70 physiolically
}

soma = Compartment(diameter = 50, length = None)
gabar = soma.add_channels(channels)

"""
another
gabar = soma.add_lgic(gabar_transit, gabar_binding, gMax=0.01, ER=0.0)
"""

gabar.set_recording(True) #gabar = soma.get_channel("LGIC")

vclamp = VClamper(baseline = -60)
vclamp.clamp(soma)

deliver = CClamper(gabar.Ligand)
xp = Experiment(soma)

xp.run(150)

deliver.set_waveform('rect', delay=20, width=80, amplitude=150)
#Waveform = Rect(delay=20, width=80, amplitude=150)
xp.Clock = 0
xp.run(150)

xp.plot_lgic()
#gabar.plot_conductance(xp.T)
plot_Vm([soma])

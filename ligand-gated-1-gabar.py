from compartment import Compartment
from lgic import LGIC, Ligand
from clampers import VClamper, CClamper, Rect
from experiment import Experiment

import matplotlib.pyplot as plt

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
gaba_r = soma.add_channels(channels)

"""
another
gaba_r = soma.add_lgic(gabar_transit, gabar_binding, gMax=0.01, ER=0.0)
"""

vclamp = VClamper(-60)
vclamp.clamp(soma)


deliver = CClamper(gaba_r.Ligand)

xp = Experiment(soma)

xp.run(150)

deliver.Waveform = Rect(delay=20, width=80, amplitude=150)
xp.run(150)
plt.figure()
plt.subplot(2,1,1)
plt.plot(xp.T, vclamp.J_ion, linewidth=2.0)
#plt.ylim([-1,0.1])

plt.subplot(2,1,2)
plt.plot(xp.T, deliver.Command, linewidth=2.0)
#plt.ylim([-0.5,1.5])

plt.show()


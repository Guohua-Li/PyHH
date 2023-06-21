from compartment import Compartment
from lgic import Ligand, LGIC
from clampers import CClamper, VClamper
from experiment import Experiment



Glu = Ligand()

# This scheme is from Krampfl et al., 2002 # in Eur J Neurosc 15:51-62.
ampar_transit = {
                 'C0': {'C1': 4.0},
                 'C1': {'C0': 2.0,   'D1': 0.15, 'C2': 2.0},
                 'C2': {'C1': 4.0,   'D2': 0.70, 'O': 20.0},
                 'D1': {'C1': 0.015, 'D2': 2.0},
                 'D2': {'C2': 0.002, 'D1': 0.875},
                 'O':  {'C2': 8.0}
                }

ampar_binding = {Glu:
                  {'C0': 'C1',
                   'C1': 'C2',
                   'D1': 'D2'
                  }
                }

AMPAR = LGIC(ampar_transit, ampar_binding)

channels = {
    AMPAR: [0.05,  0],
}

cell = Compartment(diameter = 1.5, length = 100)
ampar = cell.add_channels(channels)

vclamper = cell.add_vclamper(-60)

deliver = CClamper(ampar.Ligand)
deliver.set_waveform('rect', delay=2, width=10, amplitude=1)

xp = Experiment(cell)
xp.run(20)

xp.plot_lgic()


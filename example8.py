from pyhh import *
import pylab as plt

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

AMPAR = LGIC(ampar_transit, ampar_binding, gMax = 0.05, ER = 0)

cpm = Compartment(diameter = 1.5, length = 100)
ampar = cpm.add_channel(AMPAR)

vclp = VClamper()
vclp.Waveform = Rect(delay=0, width=150, amplitude=0)
vclp.connect(cpm)

deliver = CClamper(ampar.Ligand)
deliver.Waveform = Rect(delay=2, width=10, amplitude=1)
deliver.connect(cpm)

xp = Experiment(cpm)
xp.run(20,dt=0.005)

plt.figure()
plt.subplot(2,1,1)
plt.plot(xp.T, vclp.Jm, linewidth=2.0)
plt.ylim([-1,0.1])

plt.subplot(2,1,2)
plt.plot(xp.T, deliver.Command, linewidth=2.0)
plt.ylim([-0.5,1.5])

plt.show()


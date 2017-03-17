from PyHH import *
import pylab as plt

Glu = Ligand()

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


cpm0 = Compartment(diameter = 1.5, length = 100)
ampar0 = cpm0.add_channel(AMPAR)
cpm0.add_channel(gL)

cpm1 = Compartment(diameter = 1.5, length = 100)
ampar1 = cpm1.add_channel(AMPAR)
cpm1.add_channel(gL)

vclamp0 = VClamper()
vclamp0.Waveform = Rect(delay=0, width=150, amplitude=0)
vclamp0.connect(cpm0)

vclamp1 = VClamper()
vclamp1.Waveform = Rect(delay=0, width=150, amplitude=0)
vclamp1.connect(cpm1)

deliver0 = CClamper(ampar0.Ligand)
deliver0.Waveform = Rect(delay=2, width=10, amplitude=1)
deliver0.connect(cpm0)

deliver1 = CClamper(ampar1.Ligand)
deliver1.Waveform = Rect(delay=3, width=10, amplitude=1)
deliver1.connect(cpm1)

xp = Experiment([cpm0,cpm1])

xp.run(20,0.005)

plt.figure()
plt.subplot(2,1,1)
plt.plot(xp.T, vclamp0.Jm, linewidth=2.0)
plt.plot(xp.T, vclamp1.Jm, linewidth=2.0)
plt.ylim([-1,0.1])


plt.subplot(2,1,2)
plt.plot(xp.T, deliver0.Command, linewidth=2.0)
plt.plot(xp.T, deliver1.Command, linewidth=2.0)
plt.ylim([-0.5,1.5])

plt.show()


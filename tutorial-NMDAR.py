from pyhh import *

def func(V):
  return 5

R = TransRate(func, 'V')

def R_Mgbind(V):
  return 0.61*exp(-V/17)

def R_Reverse(V):
  return 5.4*exp(V/47)


Glu = Ligand()

k_on  = 5.0     # mM/ms
k_off = 0.0055  # /ms
alpha = 0.0916  # /ms
beta  = 0.0465  #/ms
k_r   = 0.0018  # /ms
k_d   = 0.001   # /ms

nmdar_transit = {
                 'C0': {'C1': 2*k_on},
                 'C1': {'C0': k_off, 'C2': k_on},
                 'C2': {'C1': 2*k_off, 'D0': k_r, 'O0': beta},
                 'D0': {'C2': k_d},
                 'O0': {'C2': alpha, 'M1': R_Mgbind},
                 'C3': {'C4': 2*k_on},
                 'C4': {'C3': k_off, 'C5': k_on},
                 'C5': {'C4': 2*k_off, 'D1': k_r, 'M1': beta},
                 'D1': {'C5': k_d},
                 'M1': {'C5': alpha, 'O0': R_Reverse}
                }

nmdar_binding = {Glu:
                 {'C0': 'C1',
                  'C1': 'C2'}
                }

NMDAR = LGIC(nmdar_transit, nmdar_binding, gMax=0.05, ER=60)

cpm = Compartment(1.5, 100)
nmdar = cpm.add_channels(NMDAR)

clp = VClamper()
clp.Waveform = Rect(delay=0, width=150, amplitude=0)
clp.set_baseline(-0)

clp.connect(cpm)

xp = Experiment(cpm)
xp.run(50,0.01)

dlv = CClamper(nmdar.Ligand)
dlv.Waveform = Rect(delay=5, width=100, amplitude=0.2)

xp.run(200,0.005)

import pylab as plt
plt.figure()
plt.subplot(2,1,1)
plt.plot(xp.T, clp.Jm, linewidth=2.0)
plt.ylim([-1,0.1])

plt.subplot(2,1,2)
plt.plot(xp.T, dlv.Command, linewidth=2.0)
plt.ylim([-0.5,1.5])

plt.show()


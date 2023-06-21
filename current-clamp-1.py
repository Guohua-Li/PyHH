"""
In this tutorial, we demonstrate a minimal working example of PyHH
with pre-defined sodium, potassium and chloride channels.

1. Import Compartment, Experiment, NaT, KDR and gL.

2. Define channel densities and reversal potentials.

3. Define a compartment and add channels to the compartment.

4. Define an Experiment for equation integration.

5. Run the integration.

6. Plot the result

Since there is no current stimulation, we can only see the resting potential.

If three extra lines of code are added, we can see a beautiful action potential.

"""

from compartment import Compartment
from experiment import Experiment
from gating import NaT, KDR, gL

# Please note the specific membrane conductance is given in units of nS/μm2.
chann_dist = {
    NaT: [0.6,  50], # gMax(Na) = 0.6 nS/μm2, ER (Na) = 50 mV
    KDR: [0.18,-90], # gMax(K) = 0.18 nS/μm2, ER (K) = -90 mV
    gL:  [0.03,-60], # gL = 0.03 nS/μm2, EL = -60 mV
}


cell = Compartment(diameter=40, length=None) # if length = None, the compartment is treated as a sphere
cell.add_channels(channels=chann_dist)

xp = Experiment(cell)
xp.run(t = 10) # t: time to integrate over (ms), the default integration time step (dt) is 0.005 ms

# For each run, the compartment stores the membrane potential in an array.
# In this example, there are 2000 (=10/0.005) potential samples
cell.plot_Vm(ylim=[-70, 10]) #ylim: y-axis limit

clamper = cell.add_iclamper('rect', delay = 2, width=0.5, amplitude=1.0)
xp.run(t = 16)
cell.plot_Vm(ylim=[-90, 40])

# we can change the waveform of the stimulus
clamper.set_waveform('alpha', delay = 2, tau=0.5, amplitude=0.6)
xp.run(t = 16)
cell.plot_Vm(ylim=[-90, 40])

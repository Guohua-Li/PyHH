from compartment import Compartment
from gating import NaT, KDR, gL # pre-defined ion channels
from experiment import Experiment # for equation integration

channels = {
    NaT: [0.6,  50], # gMax = 0.6 (nS/um2), Er = 50 (mV)
    KDR: [0.18,-90],
    gL:  [0.03,-60]
}

cell = Compartment(diameter=40, length=None) # if length = None, the compartment is treated as a sphere
cell.add_channels(channels)

xp = Experiment(cell)
xp.run(t = 10) # time to integrate over, the default integration time step (dt) is 0.005 ms

# For each run, the compartment stores the membrane potential in an array
# In this example, there are 2000 (=10/0.005) potential samples
cell.plot_Vm(vlimit=[-70, 10])

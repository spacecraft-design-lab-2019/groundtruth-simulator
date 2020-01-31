
import numpy as np
import matplotlib.pyplot as plt
import time
import sys, os
import sim_config as config
from simulator import Simulator

#-----------------------Clear Figures----------------------------
# clear figures
plt.close('all')


#-----------------Configuration / Parameters--------------------
tspan = np.array([0, 900*60])    # [sec]


#----------------Initialize / Setup Workspace------------------
# setup sim
config.tstep = 1
sim = Simulator(config)

# preallocate memory
T = np.arange(0, tspan[1]+config.tstep, config.tstep)
state_history = np.zeros((np.shape(T)[0], np.shape(sim.state)[0]))
eclipse_history = np.zeros((np.shape(T)[0], 1))
state_history[0, :] = sim.state


#---------------------Propagate---------------------------
t = time.time()

for i, elapsed_t in enumerate(T[0:-1]):
	sim.step(config.tstep)
	state_history[i+1, :] = sim.state
	eclipse_history[i+1, :] = sim.debug_output[2] * 1

elapsed = time.time() - t
print(elapsed)


#------------------------Plot-----------------------------
plt.figure()
plt.plot(T/3600, state_history[:,13])
plt.xlabel('time [hr]')
plt.ylabel('Temperature [K]')
plt.grid()
plt.savefig('./thermal_plots/temperature.png')

plt.figure()
plt.plot(T/3600, eclipse_history)
plt.xlabel('time [hr]')
plt.ylabel('In Eclipse')
plt.grid()
plt.savefig('./thermal_plots/eclipse.png')
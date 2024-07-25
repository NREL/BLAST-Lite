from blast.models.nmc111_gr_Kokam75Ah_2017 import Nmc111_Gr_Kokam75Ah_Battery
import matplotlib.pyplot as plt
import numpy as np




# Instantiate the cell
cell = Nmc111_Gr_Kokam75Ah_Battery()

# Run the simulation

cell.simulate_battery_life(input)

plt.plot(cell.stressors['t_days']/365, cell.outputs['q'])
plt.xlabel('Time (years)')
plt.ylabel('Relative discharge capacity')
plt.ylim((0.7, 1.02))
plt.show()

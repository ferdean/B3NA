import numpy as np
import matplotlib.pyplot as plt

from lib import *


# %% Shape functions validation
 
colors = ('#fcfdbf', '#fed395', '#fea973', '#fa7d5e', '#e95462', '#c83e73',
          '#a3307e', '#7e2482', '#59157e', '#331067', '#120d31', '#000004')

x_plot =  np.linspace(0, 1, 1000)
grid = np.linspace(0, 1, 12)

fig, ax = plt.subplots(2, 1, dpi = 330)

ax[0].scatter(grid, np.zeros(grid.shape), marker= 'x', color = 'k')
ax[1].scatter(grid, np.zeros(grid.shape), marker= 'x', color = 'k')

for i in range(12):
    
    phi = get_phi(grid, i)
    
    ax[0].plot(x_plot, phi(x_plot)[0], color= colors[i])
    ax[1].plot(x_plot, phi(x_plot)[1], color= colors[i])
    

ax[0].set_title('Odd shape functions (position)')
ax[1].set_title('Even shape functions (derivative)')

plt.tight_layout()
plt.show()
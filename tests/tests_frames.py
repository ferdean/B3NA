import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from lib import *

exec(open('../src/frame_classes.py').read())
try:
    exec(open('../src/lib.py').read())
except:
    import lib
  
# %% Bridge

name      = 'bridge_gud.txt'

directory = '../frames/' + name

x         = Structure(directory)

x.assemble_matrices()
x.solve_system()

nB     = len(x.beams)  # Number of beams
nN     = len(x.nodes)  # Number of nodes
nNb    = 2                # Number of nodes per beam
nDOFn  = 3                # Number of DOF per node
nDOFb  = nNb * nDOFn      # Number of DOF per beam

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size' : 9})

fig, ax = plt.subplots(3, 1, figsize=(5, 5), dpi = 150, sharex = False)

scaler = 5E4    

for idxBeam in range(nB):

    DOF = range(idxBeam * nDOFb, (idxBeam + 1) * nDOFb)
    
    sol_L = x.dof[DOF[0:2]]
    sol_T = x.dof[DOF[2:]]
    
    grid = np.array([0, x.beams[idxBeam].length]) # TO BE ERASED
    
    v = interp1d(grid, sol_L * scaler) 
    w = get_sol(grid, sol_T  * scaler)  
    
    x0    = x.beams[idxBeam].offset
    
    if abs(x.beams[idxBeam].direction[0]) > 1E-6:
        theta = np.arctan(x.beams[idxBeam].direction[1]/x.beams[idxBeam].direction[0]) 

    else:
        theta = np.pi/2

    if x.beams[idxBeam].direction[0] < 0:
        theta = theta + np.pi
        
    if abs((abs(x.beams[idxBeam].direction[1]) - 1)) < 1E-6 and x.beams[idxBeam].direction[1] < 0:
        theta = theta + np.pi
    
    x.beams[idxBeam].theta = theta    
    
    rotBeam  = rotateBeam((v, w), grid[-1], x0, theta)
    original = rotateBeam((lambda x: x * 0, lambda x: x * 0), grid[-1], x0, theta)
             
    x_plot = np.linspace(0, x.beams[idxBeam].length, 20)
    
    ax[2].plot(rotBeam(x_plot)[0, :],  rotBeam(x_plot)[1, :],  color = 'k')
    ax[2].plot(original(x_plot)[0, :], original(x_plot)[1, :], color = 'gray', linestyle = '--')
    ax[1].plot(original(x_plot)[0, :], original(x_plot)[1, :], color = 'k')

ax[1].arrow(4, 2, 0, -1.5, length_includes_head = True, head_width = 0.1, head_length = 0.4, fc='r', ec='r', label = 'point force', zorder=10)
ax[1].arrow(5, 2, 0, -1.5, length_includes_head = True, head_width = 0.1, head_length = 0.4, fc='r', ec='r', label = 'point force', zorder=10)
ax[1].arrow(6, 2, 0, -1.5, length_includes_head = True, head_width = 0.1, head_length = 0.4, fc='r', ec='r', label = 'point force', zorder=10)
ax[1].arrow(3, 2, 0, -1.5, length_includes_head = True, head_width = 0.1, head_length = 0.4, fc='r', ec='r', label = 'point force', zorder=10)
ax[1].arrow(7, 2, 0, -1.5, length_includes_head = True, head_width = 0.1, head_length = 0.4, fc='r', ec='r', label = 'point force', zorder=10)
ax[1].scatter([1, 9], [2, 2], marker = 'o', color = 'k', s = 10)
ax[2].scatter([1, 9], [2, 2], marker = 'o', color = 'k', s = 10)
    
ax[2].set_xlabel('x direction (m)  - deformation scaler: %.2E'%scaler)
ax[1].set_ylabel('y direction (m)')
ax[2].set_ylabel('y direction (m)')

ax[1].tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
               top= False, right= False, left=True, width = 1)
ax[2].tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
               top= False, right= False, left=True, width = 1)

ax[1].set_ylim([0, 9])
ax[2].set_ylim([0, 9])

img = mpimg.imread('bridge.jpg')
ax[0].imshow(img)

ax[0].tick_params(bottom= False, top= False, right= False, left=False, width = 1)
ax[0].axes.xaxis.set_visible(False)
ax[0].axes.yaxis.set_visible(False)

ax[0].set_title('bridge over Rush Creek (Oklahoma)')

fig.savefig('bridge_frame.png')
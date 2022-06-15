import numpy as np
import matplotlib.pyplot as plt

from lib import *

# %% Implementation

### Grid and original vector
L    = 1       
nN   = 2 
grid = np.linspace(0, L, nN)

### Deformation
def v(x):
    return 0 * x

def v_2(x):
    return 0.1 + 0*x

def w(x):
    return 0 * x

### Initial position
x0 = np.array([0, 0])


### Function
def rotateBeam(locDef, L, x0, theta):    
    """
    Returns the rotated solution function

    Parameters
    ----------
    locDef: {tuple}
        Contains the deformation functions in both directions.
            * v(x): longitudinal deformation
            * w(x): transversal deformation
    L: {scalar}
        Beam length.
    x0: {array}
        (x,y) position of the initial point.
    theta: {scalar}
        Rotation angle.

    Returns
    -------
    rotbeam: {function}
        Parametric representation of the solution.
    """
    v, w = locDef
    R    = np.array([[np.cos(theta), - np.sin(theta)],
                    [np.sin(theta),   np.cos(theta)]])

    def rotbeam(x):
        beta    = np.array([v(x) + x, w(x)])      
        globDef = x0.reshape((-1,1)) + R @ beta
        
        return globDef
    
    return rotbeam


# %% Test 1
x_plot = np.linspace(0, L, 20)

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size' : 9})

originalBeam = rotateBeam((v, w), L, x0, 0)
originalDefo = beam

x0 = np.array([0.2, 0.2])
rotBeam = rotateBeam((v, w), L, x0, np.pi/360)
rotDefo = rotateBeam((v_2, beam), L, x0, np.pi/360)

fig, ax = plt.subplots(2, 1, figsize=(5, 3), dpi = 150)

ax[0].plot(originalBeam(x_plot)[0, :], originalBeam(x_plot)[1, :], color = 'k')
ax[0].plot(x_plot, originalDefo(x_plot), color = 'gray', linestyle = '--')

ax[0].tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
        top=True, right= False, left=True, width = 1)

ax[1].plot(rotBeam(x_plot)[0, :], rotBeam(x_plot)[1, :], color = 'k')
ax[1].plot(rotDefo(x_plot)[0, :], rotDefo(x_plot)[1, :], color = 'gray', linestyle = '--')

ax[1].tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
        top=True, right= False, left=True, width = 1)


# %% Test 2

x_plot = np.linspace(0, L, 20)

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size' : 9})


x01    = np.array([0, 0])
beam_1 = rotateBeam((v, w), L, x01, np.pi/6)

x02    = np.array([beam_1(x_plot[-1])[0, 0], beam_1(x_plot[-1])[0, 1]])
beam_2 = rotateBeam((v, w), L, x02, -np.pi/12)

fig, ax = plt.subplots(figsize=(5, 3), dpi = 150)

ax.plot(beam_1(x_plot)[0, :], beam_1(x_plot)[1, :], color = 'k')
ax.plot(beam_2(x_plot)[0, :], beam_2(x_plot)[1, :], color = 'k')
ax.axvline(x=0, color="black", linestyle="-", linewidth = 5) 

ax.set_aspect('equal')

ax.tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
        top=False, right= False, left=True, width = 1)

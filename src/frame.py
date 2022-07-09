# %%
# ++++++++++++++++++++++++++++++++++++++++++++++
# +                                            +
# +     ██████╗░██████╗░███╗░░██╗░█████╗░      +
# +     ██╔══██╗╚════██╗████╗░██║██╔══██╗      +
# +     ██████╦╝░█████╔╝██╔██╗██║███████║      +
# +     ██╔══██╗░╚═══██╗██║╚████║██╔══██║      +
# +     ██████╦╝██████╔╝██║░╚███║██║░░██║      +
# +     ╚═════╝░╚═════╝░╚═╝░░╚══╝╚═╝░░╚═╝      +     
# +                                            +
# ++++++++++++++++++++++++++++++++++++++++++++++
# +                                            +
# + Authors: F. de Andrés, F. Sindy, V. Yelnyk +
# +                                            +
# ++++++++++++++++++++++++++++++++++++++++++++++

exec(open('frame_classes.py').read())
try:
    exec(open('lib.py').read())
except:
    import lib
    
# %% Part 1.- create a mesh

# name      = 'test'
# directory = '../frames/' + name

# mesher    = Mesher(directory, plotflag = False)

# %% Part 2.- Frame simulation

# name      = 'bridge_gud.txt'
name      = 'crane.txt'
directory = '../frames/' + name

x         = Structure(directory)

x.E       = 10     
x.I       = 1000    
x.A       = 1
x.mu      = 1

_, _, _   = x.assemble_matrices()

x.solve_system()
x.plot_frame(scaler = 1e-2)

# %% Part 3.- Eigenmodes

# x.eigen_freq_modes(4, 0)
# x.plot_frame(scaler = 1e4)

# %% Part 4.- Dynamics

scaler = 5e-3
lim    = np.array([1.5, 6.5])

sol, t = x.solve_dynamic(0.01, 0, 60)
x.animate_frame(xlim = lim, ylim = lim)

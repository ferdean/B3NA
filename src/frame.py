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

exec(open('src/frame_classes.py').read())
try:
    exec(open('src/lib.py').read())
except:
    import lib
    
# %% Part 1.- create a mesh

# name      = 'test'
# directory = '../frames/' + name

# mesher    = Mesher(directory, plotflag = False)

# %% Part 2.- Frame simulation

# name      = 'bridge_gud.txt'
name      = 'bridge.txt'
directory = 'frames/' + name

x         = Structure(directory)

x.E       = 10     
x.I       = 1000    
x.A       = 1
x.mu      = 1

_, _, _   = x.assemble_matrices()

x.solve_system()
x.plot_frame(scaler = 1e-2)

# %%part 3 Dynamics
#scaler = 1e-2

#lim = np.array([1.5, 6.5])

#sol, t = x.solve_dynamic(0.01, 0, 60)
#x.animate_frame(xlim = lim, ylim = lim)

# %% Part 4.- Eigenmodes
x.eigen_freq_modes(2)
x.plot_frame(scaler = 1e3)

# %%part 5 Dynamics Eigenmodes
scaler = 1e3
sol = x.eigen_freq_modes(0, dynamic = True,t_0 = 0, t_f = 10, Nt = 200, modes = np.array([1]))
x.animate_frame(xlim = (0,10), ylim = (0,9))

x.ani.save('test2.gif', writer='imagemagick', fps= 30)
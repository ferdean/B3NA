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

test,_,_ = x.assemble_matrices()

x.solve_system()

x.plot_frame(scaler = 5E2)

# %%part 3 Dynamics
name      = 'crane.txt'

directory = '../frames/' + name

x         = Structure(directory)

x.assemble_matrices()
sol,t = x.solve_dynamic(0.1,0,50)

x.animate_frame()
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

# name      = 'frame_3'
# directory = '../frames/' + name

# mesher    = Mesher(directory)

# %% Part 2.- Frame simulation

name      = 'crane.txt'
directory = '../frames/' + name

x         = Structure(directory)

x.assemble_matrices()
x.solve_system()

x.plot_frame(scaler = 5E3)
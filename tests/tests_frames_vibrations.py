exec(open('src/frame_classes.py').read())
try:
    exec(open('src/lib.py').read())
except:
    import lib

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

# %%Eigenmodes
for i in range(3): 
    x.eigen_freq_modes(i+1)
    name_ = "Eigenmode_"+str(i+1)+"_Bridge.png"
#    x.plot_frame(scaler = 1e3,save = True, name = name_)

# %%Dynamics Eigenmodes
for i in range(3):
    scaler = 1e3
    sol = x.eigen_freq_modes(0, dynamic = True,t_0 = 0, t_f = 10, Nt = 200, modes = np.array([1]))
    x.animate_frame(xlim = (0,10), ylim = (0,9))
#    x.ani.save('test2.gif', writer='imagemagick', fps= 30)
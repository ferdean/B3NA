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

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

try:
    import Tkinter as tk
except:
    import tkinter as tk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

try:
    exec(open('lib.py').read())
except:
    import lib

### Initialization
def main():
    root = tk.Tk()
    GUI  = mesher(root)
    GUI.root.mainloop()
    
    return None

### Meshing window
class mesher:
    def __init__(self, root):
        self.root = root
        self.root.title('B3NA mesher')
        self.root.geometry('1200x500')
        self.root.configure(background='white')  
                
        self.root.columnconfigure(0, weight= 10)
        self.root.columnconfigure(1, weight= 1)
        self.root.columnconfigure(2, weight= 1)
        self.root.columnconfigure(3, weight= 1)
        self.root.columnconfigure(4, weight= 1)
        
        self.mesh = np.array([0, 1])
        
        self.L   = 1.0      # [m]
        self.nN  = 30       # [-]
      
        titleFont    = ("CMU Sans Serif", 14, "bold")
        textFont     = ("CMU Sans Serif", 10)
                
        ### Meshing
        tk.Label(self.root, text = 'Mesh parameters', font = titleFont, bg = 'white', bd = 2).grid(row = 0, columnspan = 4, column = 1)
        tk.Label(self.root, text = 'Length: ',     font = textFont,  bg = 'white', justify = 'right').grid(row = 1, column = 1)
        self.L_input = tk.Entry(self.root, width = 10)        
        self.L_input.grid(row = 1, columnspan = 2, column = 2)
        tk.Label(self.root, text = '[m]', font = textFont, bg = 'white').grid(row = 1, column = 4)

        tk.Label(self.root, text = 'Number of nodes: ',     font = textFont,  bg = 'white', justify = 'right').grid(row = 2, column = 1)
        self.nN_input = tk.Entry(self.root, width = 10)        
        self.nN_input.grid(row = 2, columnspan = 2, column = 2)
        tk.Label(self.root, text = '[-]', font = textFont, bg = 'white').grid(row = 2, column = 4)
        
        ### Action buttons   
        
        regularButton  = tk.Button(self.root, text = 'Regular mesh', font = textFont, command = self.regularMesh)
        regularButton.grid(row = 3, columnspan = 2, column = 1)
   
        chebyButton = tk.Button(self.root, text = 'Chebyshev mesh', font = textFont, command = self.chebyshevMesh)
        chebyButton.grid(row = 3, columnspan = 2, column = 3)
        
        btn = tk.Button(self.root, text = 'Accept', bg = 'lightgreen', font = textFont, width = 20, command = self.openSolver)
        btn.grid(row = 5, columnspan = 4, column = 1)               
        
        self.plot()

    def updateParams(self, event = None):
        
        self.nN  = int(self.nN_input.get())
        self.L   = float(self.L_input.get())
        
        return None
    
    def regularMesh(self):
        
        self.updateParams()
        self.mesh = np.linspace(0, self.L, self.nN)
        
        self.plot()
        
        return None
    
    def chebyshevMesh(self):
        
        self.updateParams()
        
        t = (2 * np.arange(self.nN -2 ) + 1 ) / ( 2 * ( self.nN -2 ) ) * np.pi       
        self.mesh = np.hstack((0, 0.5 * (self.L) * (1 + np.cos(t)), self.L))       
       
        self.plot()
        
        return None
        
    def plot(self):
        
        figure = plotMesh(self.mesh)
        
        chart = FigureCanvasTkAgg(figure, self.root)
        chart.get_tk_widget().grid(rowspan = 10, row = 0, column = 0) 
    
        return None
        
    def openSolver(self):
        solverWindow = tk.Toplevel(self.root)
        
        GUI = window(solverWindow, self)
        GUI.root.withdraw()

        solverWindow.deiconify()

        return None



### Main class
class window:
    def __init__(self, root, mesher):
        
        self.root = root
        self.root.title('B3NA - Benouilli beam bending numerical analysis')
        self.root.geometry('1200x500')
        self.root.configure(background='white')
        
        self.root.columnconfigure(0, weight= 10)
        self.root.columnconfigure(1, weight= 1)
        self.root.columnconfigure(2, weight= 1)
        self.root.columnconfigure(3, weight= 1)
        self.root.columnconfigure(4, weight= 1)
        
        self.E   = 210      # [N/mm2]
        self.I   = 3.3e7    # [mm4]
        self.mu  = 20       # [kg/m]
        
        self.T  = 10
        self.h  = 1E-3
                
        self.mesh = mesher.mesh
        self.L    = mesher.L
        self.nN   = mesher.nN
        
        self.staticSol = np.zeros((self.nN * 2 + 2,))
        self.u_1_0 = np.zeros(self.staticSol.shape)
        self.u_2_0 = np.zeros(self.staticSol.shape)
        
        titleFont    = ("CMU Sans Serif", 14, "bold")
        textFont     = ("CMU Sans Serif", 10)
        
        ### Material properties
        tk.Label(self.root, text = 'Material properties', font = titleFont, bg = 'white', bd = 2).grid(row = 0, columnspan = 4, column = 1)
        tk.Label(self.root, text = 'Young modulus: ',     font = textFont,  bg = 'white', justify = 'right').grid(row = 1, column = 1)
        self.E_input = tk.Entry(self.root, width = 10)
        self.E_input.grid(row = 1, columnspan = 2, column = 2)
        tk.Label(self.root, text = '[N/mm^2]', font = textFont, bg = 'white').grid(row = 1, column = 4)
        
        tk.Label(self.root, text = 'Inertia: ', font = textFont, bg = 'white', justify = 'right').grid(row = 2, column = 1)
        self.I_input = tk.Entry(self.root, width = 10)
        self.I_input.grid(row = 2, columnspan = 2, column = 2)
        tk.Label(self.root, text = '[mm^4]', font = textFont, bg = 'white').grid(row = 2, column = 4)
                
        tk.Label(self.root, text = 'Density: ', font = textFont,  bg = 'white', justify = 'right').grid(row = 3, column = 1)
        self.mu_input = tk.Entry(self.root, width = 10)
        self.mu_input.grid(row = 3, columnspan = 2, column = 2)
        tk.Label(self.root, text = '[kg/m]', font = textFont, bg = 'white').grid(row = 3, column = 4)
        
        ### Boundary Conditions
        tk.Label(self.root, text = 'Boundary conditions', font = titleFont, bg = 'white').grid(row = 4, columnspan = 4, column = 1)
        tk.Label(self.root, text = 'a: ', font = textFont, bg = 'white', justify = 'right').grid(row = 5, column = 1)
        tk.Label(self.root, text = 'b: ', font = textFont, bg = 'white', justify = 'right').grid(row = 5, column = 3)
        tk.Label(self.root, text = 'QL: ', font = textFont, bg = 'white', justify = 'right').grid(row = 6, column = 1)
        tk.Label(self.root, text = 'ML: ', font = textFont, bg = 'white', justify = 'right').grid(row = 6, column = 3)
        
        self.a_input = tk.Entry(self.root, width = 5)
        self.a_input.grid(row = 5, column = 2)
        self.b_input = tk.Entry(self.root, width = 5)
        self.b_input.grid(row = 5, column = 4)
        self.QL_input = tk.Entry(self.root, width = 5)
        self.QL_input.grid(row = 6, column = 2)
        self.ML_input = tk.Entry(self.root, width = 5)
        self.ML_input.grid(row = 6, column = 4)
        
        ### Simulation properties
        tk.Label(self.root, text = 'Simulation properties', font = titleFont, bg = 'white').grid(row = 7, columnspan = 4, column = 1)
        tk.Label(self.root, text = 'Tmax (s): ', font = textFont, bg = 'white', justify = 'right').grid(row = 8, column = 1)
        tk.Label(self.root, text = 'Stepsize (-): ', font = textFont, bg = 'white', justify = 'right').grid(row = 8, column = 3)
        
        self.T_input = tk.Entry(self.root, width = 5)
        self.T_input.grid(row = 8, column = 2)
        self.h_input = tk.Entry(self.root, width = 5)
        self.h_input.grid(row = 8, column = 4)
        
        
        ### Action buttons        
        tk.Label(self.root, text = 'Solver', font = titleFont, bg = 'white').grid(row = 9, columnspan = 4, column = 1)
              
        staticButton  = tk.Button(self.root, text = 'Compute static def.', font = textFont, command = self.static)
        staticButton.grid(row = 10, columnspan = 2, column = 1)
   
        dynamicButton = tk.Button(self.root, text = 'Simulate', font = textFont, command = self.dynamic)
        dynamicButton.grid(row = 10, columnspan = 2, column = 3)
        
        self.root.bind("<Return>", self.updateParams)
        self.plot()
        
        ### State label     
        self.stateLabel = tk.StringVar()
        self.stateLabel.set('')
        
        self.stateLabelColor = 'red'
        
        self.stateLabelHandle = tk.Label(self.root, textvariable = self.stateLabel, font = textFont, bg = 'white', fg = self.stateLabelColor).grid(row = 12, columnspan = 4, column = 1)

   
    def updateParams(self, event = None):
        
        self.E  = float(self.E_input.get())
        self.I  = float(self.I_input.get())
        self.mu = float(self.mu_input.get())
        
        self.a  = float(self.a_input.get())
        self.b  = float(self.b_input.get())
        self.QL = float(self.QL_input.get())
        self.ML = float(self.ML_input.get())
        
        self.BC  = (self.a, self.b, self.QL, self.ML)
                
        self.plot()
        
        return None
    
    def updateParamsDynamic(self, event = None):
        
        self.T  = float(self.T_input.get())
        self.h  = float(self.h_input.get())
        
        return None

    def static(self):
        
        self.updateParams()
        self.stateLabel.set('Computing...')
        
        def q(x):
            return - 1000 * x

        S, M  = getMatrices(self.mesh, self.E, self.I, self.mu, quadrature = True)
        
        RHS   = getRHS(self.mesh, q)       
        
        nN = self.nN
        
        e0 = np.zeros(nN*2);    e0[0]  = 1.0
        eL = np.zeros(nN*2);    eL[-1] = 1.0
        d0 = np.zeros(nN*2);    d0[1]  = 1.0
        dL = np.zeros(nN*2);    dL[-2] = 1.0
        
        # Apply BCs
        self.Me, self.Se, self.RHSe = fixBeam(M, S, RHS, (e0, eL), (d0, dL), self.BC)
        
        # Solve
        self.staticSol      = sparse.linalg.spsolve(self.Se, self.RHSe)
        
        self.stateLabel.set('Static solution available')
        
        self.plot()
        
        return None

    
    def dynamic(self):
        
        self.updateParamsDynamic()
        self.stateLabel.set('Computing...')
        
        self.t0 = 0  
        
        self.RHSe = np.zeros(self.staticSol.shape) 
               
        self.initialConds = (self.staticSol, self.u_1_0, self.u_2_0)
        
        self.dynamicSol, self.time     = newmarkMethod(self.Me, self.Se, self.RHSe, 
                                                       self.initialConds, self.h, self.t0, self.T, 
                                                       verbose = False)
        self.stateLabel.set('Done')
        self.stateLabelHandle = tk.Label(self.root, textvariable = self.stateLabel, font = ("CMU Sans Serif", 10), bg = 'lightgreen', fg = 'black').grid(row = 12, columnspan = 4, column = 1)
      
        self.animate()

        return None
    
    def plot(self):
        
        figure, self.ymax = plotBeam(self.mesh, self.staticSol[:-2], 50, -1)
        
        chart = FigureCanvasTkAgg(figure, self.root)
        chart.get_tk_widget().grid(rowspan = 12, row = 0, column = 0)
        
    def animate(self):
        
        self.fig = plt.Figure(figsize = (5, 3), dpi = 150)
        
        canvas = FigureCanvasTkAgg(self.fig, master= self.root)
        canvas.get_tk_widget().grid(rowspan = 12, row = 0, column = 0)
        
        self.ax = self.fig.add_subplot(111)
                
        x_plot = np.linspace(0, 1, 50)

        beam  = get_sol(self.mesh, self.dynamicSol[0:-2, 0])
        line, = self.ax.plot(x_plot, beam(x_plot)*1e3, color= '#808080')

        self.ax.plot([self.mesh.min(), self.mesh.max()], [beam(x_plot)[0]*1e3, beam(x_plot)[0]*1e3], color= '#959595', linestyle= '--')   
        self.ax.axvline(x = 0, color="black", linestyle="-", linewidth = 5)
        
        self.ax.set_ylabel('deformation (mm)')
        self.ax.set_xlabel('x-dimension (-)')
        
        self.ax.tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
            top=True, right= False, left=True, width = 1)
            
        self.ax.set_ylim((-self.ymax, self.ymax))
                
        text = self.ax.text(0.475 * self.L, 1.1 * self.ymax, 't = 0.00 s')
        
        def animationFrame(i):
            beam = get_sol(self.mesh, self.dynamicSol[0:-2, i+1])    
            line.set_ydata(beam(x_plot) * 1e3)  # update the data
            text.set_text('t = %.2f s'%(self.time[i]))
            return line, text,
        
        self.ani = animation.FuncAnimation(self.fig, animationFrame, np.arange(0, self.dynamicSol.shape[1]), interval = 200, blit=False)
                
        
### TODO: Click on mesh to set up forces interactively. Some references:
#           * https://stackoverflow.com/questions/27565939/getting-the-location-of-a-mouse-click-in-matplotlib-using-tkinter
#           * https://stackoverflow.com/questions/25521120/store-mouse-click-event-coordinates-with-matplotlib
#           * https://matplotlib.org/3.1.0/gallery/user_interfaces/embedding_in_tk_sgskip.html

### TODO: Allow non-constant material properties

### TODO: Posibility to specify force
        
main()








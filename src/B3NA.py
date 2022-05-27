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
import tkinter as tk
from lib import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


### Testing (to be erased)

L  = 1
nN = 30

grid = np.linspace(0, L, nN)




### Initialization
def main():
    root = tk.Tk()
    GUI  = window(root)
    GUI.root.mainloop()
    
    return None

### Main class
class window:
    def __init__(self, root):
        self.root = root
        self.root.title('B3NA - Benouilli beam bending numerical analysis')
        self.root.geometry('1200x500')
        self.root.configure(background='white')
        
        self.root.columnconfigure(0, weight= 12)
        self.root.columnconfigure(1, weight= 1)
        self.root.columnconfigure(2, weight= 1)
        self.root.columnconfigure(3, weight= 1)
        
        self.E   = 210      # [N/mm2]
        self.I   = 3.3e7    # [mm4]
        self.L   = 1        # [m]
        self.nN  = 30       # [-]
        self.mu  = 20       # [kg/m]
        
        self.sol = np.zeros((self.nN * 2 + 2,))
        
        u_1_0 = np.zeros(sol.shape)
        u_2_0 = np.zeros(sol.shape)
        
        self.initialConds = (self.sol, u_1_0, u_2_0)
                
        self.grid = np.linspace(0, self.L, self.nN)
     
        titleFont    = ("CMU Sans Serif", 14, "bold")
        textFont     = ("CMU Sans Serif", 10)
        
        ### Material properties
        tk.Label(self.root, text = 'Material properties', font = titleFont, bg = 'white', bd = 2).grid(row = 0, columnspan = 3, column = 1)
        tk.Label(self.root, text = 'Young modulus: ',     font = textFont,  bg = 'white', justify = 'left').grid(row = 1, column = 1)
        self.E_input = tk.Entry(self.root, width = 10)
        self.E_input.grid(row = 1, column = 2)
        tk.Label(self.root, text = '[N/mm2]', font = textFont, bg = 'white').grid(row = 1, column = 3)
        
        tk.Label(self.root, text = 'Inertia: ', font = textFont, bg = 'white', justify = 'left').grid(row = 2, column = 1)
        self.I_input = tk.Entry(self.root, width = 10)
        self.I_input.grid(row = 2, column = 2)
        tk.Label(self.root, text = 'mm4', font = textFont, bg = 'white').grid(row = 2, column = 3)
        
        tk.Label(self.root, text = 'Length: ',font = textFont,  bg = 'white', justify = 'left').grid(row = 3, column = 1)
        self.L_input = tk.Entry(self.root, width = 10)
        self.L_input.grid(row = 3, column = 2)
        tk.Label(self.root, text = 'm', font = textFont, bg = 'white').grid(row = 1, column = 3)
        
        tk.Label(self.root, text = 'Density: ', font = textFont,  bg = 'white', justify = 'left').grid(row = 4, column = 1)
        self.mu_input = tk.Entry(self.root, width = 10)
        self.mu_input.grid(row = 4, column = 2)
        tk.Label(self.root, text = 'kg/m', font = textFont, bg = 'white').grid(row = 4, column = 3)
               
        ### Action buttons        
        tk.Label(self.root, text = 'Solver', font = titleFont, bg = 'white').grid(row = 5, columnspan = 3, column = 1)
              
        staticButton  = tk.Button(self.root, text = 'Compute static def.', font = textFont, command = self.static)
        staticButton.grid(row = 6, columnspan = 2, column = 1)
   
        dynamicButton = tk.Button(self.root, text = 'Simulate', font = textFont, command = self.dynamic)
        dynamicButton.grid(row = 6, columnspan = 2, column = 2)
        
        self.root.bind("<Return>", self.updateParams)
        self.plot()
        
        ### State label     
        self.stateLabel = tk.StringVar()
        self.stateLabel.set('')
        
        self.stateLabelColor = 'red'
        
        self.stateLabelHandle = tk.Label(self.root, textvariable = self.stateLabel, font = textFont, bg = 'white', fg = self.stateLabelColor).grid(row = 7, columnspan = 3, column = 1)

   
    def updateParams(self, event = None):
        
        self.E  = float(self.E_input.get())
        self.I  = float(self.I_input.get())
        self.L  = float(self.L_input.get())
        self.mu = float(self.mu_input.get())
        
        self.plot()
        
        return None

    def static(self):
        
        self.updateParams()
        self.stateLabel.set('Computing...')
        
        def q(x):
            return 1000 * x

        S, M  = getMatrices(self.grid, self.E, self.I, self.mu, quadrature = True)
        
        RHS   = getRHS(grid, q)       
        
        BC   = (0, 0, 0, 0)
        
        e0 = np.zeros(nN*2);    e0[0]  = 1.0
        eL = np.zeros(nN*2);    eL[-1] = 1.0
        d0 = np.zeros(nN*2);    d0[1]  = 1.0
        dL = np.zeros(nN*2);    dL[-2] = 1.0
        
        # Apply BCs
        self.Me, self.Se, self.RHSe = fixBeam(M, S, RHS, (e0, eL), (d0, dL), BC)
        
        # Solve
        self.sol      = sparse.linalg.spsolve(Se, RHSe)
        
        self.stateLabel.set('Static solution available')
        
        self.plot()
        
        return None

    
    def dynamic(self):
        
        # La GUI necesita una zona para especificar las condiciones de la simulación
        # (t0, T y h)
        self.updateParams()
        self.stateLabel.set('Computing...')
        
        self.t0 = 0
        self.T  = 10
        self.h  = 1E-3
        
        nS    = int((self.T - self.t0)//self.h)
        
        for idx in range(nS):
            self.sol, self.u_1, self.u_2 = Newmarkmethod_step(self.sol, self.u_1, self.u_2, self.h, self.Me, self.Se, self.RHSe, beta = 1/4, gamma = 1/2)           
            self.plot()
            
            print(idx)

        return None
    
    def plot(self):
        
        # self.sol = 
        
        figure   = plotBeam(self.grid, self.sol[:-2], 50, (-1.4e-5, 1.4e-5))
        
        chart = FigureCanvasTkAgg(figure, self.root)
        chart.get_tk_widget().grid(rowspan = 6, row = 0, column = 0)
    
main()








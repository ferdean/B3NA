import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from lib import *
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as Tk

np.set_printoptions(precision = 2)

class Mesher:
    def __init__(self, directory, whole = True, zeroindex = True, plotflag = False):
        
        self.indexing  = 0 if zeroindex else 1
        self.directory = directory
        self.plotflag  = plotflag
        
        plt.rcParams['text.usetex'] = True
        plt.rcParams.update({'font.size' : 9})
        
        fig, ax = plt.subplots(figsize=(5, 3), dpi = 150)
        
        ax.title.set_text('Add nodes by clicking')
        
        ax.plot(0, 0)
        
        ax.set_xlim([0, 20])
        ax.set_ylim([0, 20])
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        
        ax.tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
                      top=True, right= False, left=True, width = 1)       
        
        ax.grid()
        
        self.coords = []
    
        def onclick(event):
            
            ix, iy = event.xdata, event.ydata
        
            if ix is None or iy is None:
                return
            
            if whole:
                ix = round(ix); iy = round(iy)
                
            if (ix, iy) in self.coords:
                return
            
            self.coords.append((ix, iy))
        
            ax.plot(ix, iy, marker = 'o', color = 'k')
            ax.set_xlim([0, 20])
            ax.set_ylim([0, 20])
            fig.canvas.draw()
            
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        
        plt.show()
        
        if len(self.coords) == 0:
            exit()

        fig_edge, ax_edge = plt.subplots(figsize=(5, 3), dpi = 150)
        
        ax_edge.title.set_text('Draw edges')
        
        x = [x[0] for x in self.coords]
        y = [y[1] for y in self.coords]
        
        ax_edge.plot(x, y, "o", color="black")
        
        ax_edge.set_xlim([0, 20])
        ax_edge.set_ylim([0, 20])
        ax_edge.set_xlabel('x (m)')
        ax_edge.set_ylabel('y (m)')

        ax_edge.tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
                            top=True, right= False, left=True, width = 1)  
        
        for i,(x,y) in enumerate(self.coords):
            ax_edge.annotate(i + self.indexing, (x + 0.1, y + 0.1))
        
        self.edges        = []
        line_coords       = []
        line_indices      = []
        node_indices      = []

        def onclick(event):
            ix, iy = event.xdata, event.ydata
            try:
                node_index = self.find_closest((ix,iy))
                node_indices.append(node_index + self.indexing)
                
            except TypeError:
                return
            
            line_coords.append(self.coords[node_index])
            
            if len(line_coords) == 2:
                if sorted([x for x in node_indices]) not in [sorted(x) for x in self.edges]:
                    ax_edge.plot([line_coords[0][0], line_coords[1][0]], [line_coords[0][1],line_coords[1][1]], color="black")
                    ax_edge.set_xlim([0, 20])
                    ax_edge.set_ylim([0, 20])
                    
                    self.edges.append((node_indices[0],node_indices[1]))
                
                line_coords.pop()
                line_coords.pop()
                node_indices.pop()
                node_indices.pop()
        
            fig_edge.canvas.draw()

        cid = fig_edge.canvas.mpl_connect('button_press_event', onclick)
        
        plt.show()

        
        self.classify()
    

    def find_closest(self, point):

        min_length = None
        index      = None

        for i, coord in enumerate(self.coords):
            
            length = np.linalg.norm(np.array(coord)-np.array(point))
            
            if min_length is None or length < min_length:
                min_length = length
                index = i
                
        return index
    
    def classify(self):  
        print('Determine type of each node:\n 1 -> FIXED\n 2 -> FREE\n 3 -> FORCE\n 4 -> MOVABLE')
        
        allowed_vals  = ['1','2','3','4']
        typedict      = {'1':'FIXED','2':'FREE','3':'FORCE','4':'MOVABLE'}
        self.typesstr = []
        self.types    = []
        
        print(len(self.coords))
        
        for i in range(self.indexing, len(self.coords) + self.indexing):
            nr = ''
            
            while nr not in allowed_vals:
                nr = input(f'Node {i}: ')
                
            if nr == '3':
                Fx = input('F_x = ')
                Fy = input('F_y = ')
                M  = input('M = ')
                
                self.typesstr.append(f'{i} FORCE [{Fx} {Fy}] [{M}]\n')
                
            elif nr == '4':
                x = input('x-val = ')
                y = input('y-val = ')
                
                self.typesstr.append(f'{i} MOVABLE [{x} {y}]\n')
                
            else:
                self.typesstr.append(f'{i} {typedict[nr]}\n')
                
            self.types.append(typedict[str(nr)])
            
        self.savefile()
        
        return None
        
    def savefile(self):
        
        with open(self.directory + '.txt', 'w') as f:
            f.write('NODES\n')
            for coord in self.coords:
                f.write(f'{coord[0]} {coord[1]}\n')
            f.write('BEAMS\n')
            for edge in self.edges:
                f.write(f'{edge[0]} {edge[1]}\n')               
            f.write('TYPE\n')
            for type_ in self.typesstr:
                f.write(type_)

        if self.plotflag:
            self.plot_mesh()
            
    def plot_mesh(self):
        
        plt.rcParams['text.usetex'] = True
        plt.rcParams.update({'font.size' : 9})
        
        fig, ax = plt.subplots(figsize=(5, 3), dpi = 150)
        
        ax.set(xlim=(0,10),ylim=(0,10))
        ax.grid()
        
        colors = {'FIXED':'black','FREE':'blue','FORCE':'red','MOVABLE':'green'}
       
        fig.text(0.2,0.9,'FIXED',ha='center',va='bottom',size='large',color=colors['FIXED'])
        fig.text(0.4,0.9,'FREE',ha='center',va='bottom',size='large',color=colors['FREE'])
        fig.text(0.6,0.9,'FORCE',ha='center',va='bottom',size='large',color=colors['FORCE'])
        fig.text(0.8,0.9,'MOVABLE',ha='center',va='bottom',size='large',color=colors['MOVABLE'])
        
        for edge in edges:
            x1,y1 = self.coords[edge[0] - self.indexing]
            x2,y2 = self.coords[edge[1] - self.indexing]
            
            ax.plot([x1,x2], [y1,y2], color='black')
            
        for (i,coord,type) in zip(range(len(self.coords)), self.coords, self.types):
            ax.plot(coord[0],coord[1],color = colors[type], marker='o')
            
            ax.annotate(i + self.indexing, (coord[0] + 0.1, coord[1] + 0.1), color='black')
            
        plt.savefig(self.directory + '.png')


class Structure:
    def __init__(self, filename):
        
        self.nodes     = []
        self.beams     = []
        
        self.C_matrix  = None   # Constraints
        self.S_matrix  = None   # Free stiffness
        self.Se_matrix = None   # Constrained stiffness
        self.M_matrix  = None   # TODO
        self.RHS       = None   # Right-hand side
        
        self.dof       = 0         # Solved degrees of freedom
        self.E         = 210      # [N/mm2]
        self.I         = 3.3E7    # [mm4]
        self.A         = 22E4     # [mm2]
        self.mu        = 0.1      # [N s^2/mm^4]
        
        mode          = 0
        
        f = open(filename, "r")  # Read the file containing the mesh

        for line in f:

            if line[0] == "#":
                continue
            
            if line.strip() == "NODES:":
                mode = 1
                node_index = 0
                continue

            if line.strip() == "BEAMS:":
                mode = 2
                beam_index = 0
                continue

            if mode == 1:
                content = line.strip().split()
                coord = np.array(
                    [float(content[0]), float(content[1])])  # x,y coordinates of a node in a global coordinate system
                status = content[2]  # get status of the node
                node = Node(node_index, coord, status)  # create node object
                if status == "FORCE":  # if there is a force applied to the node, store it in the node object
                    node.force = np.array([float(content[3]), float(content[4]), float(content[5])])
                self.nodes.append(node)  # store this node in a list with all other nodes in the structure
                node_index += 1
                
            if mode == 2:
                content = line.strip().split()
                nodes = (self.nodes[int(content[0])], self.nodes[int(content[1])])  # 2 nodes that define a beam
                beam = Beam(beam_index, nodes)  # create beam object
                nodes[0].beams.append(beam)  # add a reference to this beam for every node that is a part of this beam
                nodes[1].beams.append(beam)
                # print(beam.index)
                self.beams.append(beam)  # store this beam in a list with all other beams in the structure
                beam_index += 1
                
        f.close()
                  
   
    def plot_frame(self, scaler = 1e6):
        
        nB     = len(self.beams)  # Number of beams
        nN     = len(self.nodes)  # Number of nodes
        nNb    = 2                # Number of nodes per beam
        nDOFn  = 3                # Number of DOF per node
        nDOFb  = nNb * nDOFn      # Number of DOF per beam
        
        plt.rcParams['text.usetex'] = True
        plt.rcParams.update({'font.size' : 9})
        fig, ax = plt.subplots(figsize=(5, 3), dpi = 150)
            
        for idxBeam in range(nB):
        
            DOF = range(idxBeam * nDOFb, (idxBeam + 1) * nDOFb)
            
            sol_L = self.dof[DOF[0:2]]
            sol_T = self.dof[DOF[2:]]
            
            grid = np.array([0, self.beams[idxBeam].length]) # TO BE ERASED
            
            v = interp1d(grid, sol_L * scaler) 
            w = get_sol(grid, sol_T  * scaler)  
            
            x0    = self.beams[idxBeam].offset
            
            if abs(self.beams[idxBeam].direction[0]) > 1E-6:
                theta = np.arctan(self.beams[idxBeam].direction[1]/self.beams[idxBeam].direction[0]) 

            else:
                theta = np.pi/2

            if self.beams[idxBeam].direction[0] < 0:
                theta = theta + np.pi
                
            if abs((abs(self.beams[idxBeam].direction[1]) - 1)) < 1E-6 and self.beams[idxBeam].direction[1] < 0:
                theta = theta + np.pi
            
            self.beams[idxBeam].theta = theta    
            
            rotBeam  = rotateBeam((v, w), grid[-1], x0, theta)
            original = rotateBeam((lambda x: x * 0, lambda x: x * 0), grid[-1], x0, theta)
                     
            x_plot = np.linspace(0, self.beams[idxBeam].length, 20)
            
            ax.plot(rotBeam(x_plot)[0, :],  rotBeam(x_plot)[1, :],  color = 'k')
            ax.plot(original(x_plot)[0, :], original(x_plot)[1, :], color = 'gray', linestyle = '--')
            
        ax.set_xlabel('x direction (m)')
        ax.set_ylabel('y direction (m)')
        
        ax.tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
                       top=True, right= False, left=True, width = 1)
        
        ax.set_title('(deformation scaler: %.2E)'%scaler, fontsize = 8)

        plt.show()
    

    def is_origin(self, node, beam):
        
        if beam.nodes[0] == node:
            return 1
        
        elif beam.nodes[1] == node:
            return 2
        
        else:
            print("node doesnt belong to this beam")
            return 0

    def assemble_matrices(self):
        n_beams = len(self.beams)
        S = np.zeros((n_beams * 6, n_beams * 6))  # global stiffness matrix
        M_tilde = np.zeros((n_beams * 6, n_beams*6))

        M1_loc = np.array([[1/3, 1/6],[1/6, 1/3]])
        M2_loc = np.array([[ 0.37142857,  0.05238095,  0.12857143, -0.03095238],
                           [ 0.05238095,  0.00952381,  0.03095238, -0.00714286],
                           [ 0.12857143,  0.03095238,  0.37142857, -0.05238095],
                           [-0.03095238, -0.00714286, -0.05238095,  0.00952381]])

        S1_loc = np.array([[1.0, -1.0], [-1.0, 1.0]])
        S2_loc = np.array([[12., 6., -12., 6.],
                           [6., 4., -6., 2.],
                           [-12., -6., 12., -6.],
                           [6., 2., -6., 4.]])

        mode1 = np.array([[1, 0, 1, 0],
                          [0, 0, 0, 0],
                          [1, 0, 1, 0],
                          [0, 0, 0, 0]])

        mode2 = np.array([[0, 1, 0, 1],
                          [1, 0, 1, 0],
                          [0, 1, 0, 1],
                          [1, 0, 1, 0]])

        mode3 = np.array([[0, 0, 0, 0],
                          [0, 1, 0, 1],
                          [0, 0, 0, 0],
                          [0, 1, 0, 1]])

        # loop through beams and construct matrix S:
        for beam in self.beams:
            global_index = beam.index * 6
            h = beam.length
            h_array = (mode1 * h ** -3) + (mode2 * h ** -2) + (mode3 * h ** -1)
            S[global_index:global_index + 2, global_index:global_index + 2] = self.E * self.A * S1_loc / h
            S[global_index + 2:global_index + 6, global_index + 2:global_index + 6] = self.E * self.I * S2_loc * h_array

            M_tilde[global_index:global_index + 2, global_index:global_index + 2] = self.mu * M1_loc * h
            M_tilde[global_index + 2:global_index + 6, global_index + 2:global_index + 6] = self.mu *  M2_loc * h_array * h**4

        # loop through nodes and get constraints and forces
        RHS = np.zeros(n_beams * 6)
        C = np.zeros((n_beams * 6, 0))

        for node in self.nodes:
            if node.status == "FREE" or node.status == "FORCE":

                anchor_beam = node.beams[0]  # beam to which all other beams in this node will be pairwise attached
                cos_phi1 = anchor_beam.direction[0]
                sin_phi1 = anchor_beam.direction[1]
                local_index_1 = self.is_origin(node, anchor_beam)
                global_index_1 = anchor_beam.index * 6  # where dof of this vector are located in the global dof vector (every beam has 6 dof)
                index_1_v = global_index_1 + local_index_1 - 1  # location of v(0) or v(L) dof (depending if the node is origin of this beam or not)  in the global dof vector
                index_1_w = global_index_1 + local_index_1 * 2  # location of w(0) or w(L) dof (depending if the node is origin of this beam or not)  in the global dof vector
                index_1_w_prime = global_index_1 + local_index_1 * 2 + 1  # location of w'(0) or w'(L) dof ................ in the global dof vector

                if node.status == "FORCE":
                    fx = node.force[0]
                    fy = node.force[1]
                    M = node.force[2]
                    RHS[index_1_v] += fx * cos_phi1 + fy * sin_phi1  # I don't know if this is correct !!!!!!
                    RHS[index_1_w] += -fx * sin_phi1 + fy * cos_phi1  # I don't know if this is correct !!!!!!
                    RHS[index_1_w_prime] += M  # reserved for the moment

                for beam in node.beams[1:]:
                    cos_phi2 = beam.direction[0]
                    sin_phi2 = beam.direction[1]
                    local_index_2 = self.is_origin(node, beam)
                    global_index_2 = beam.index * 6  # where dof of this vector are located in the global dof vector (every beam has 6 dof)
                    index_2_v = global_index_2 + local_index_2 - 1  # location of v(0) or v(L) dof (depending if the node is origin of this beam or not)  in the global dof vector
                    index_2_w = global_index_2 + local_index_2 * 2  # location of w(0) or w(L) dof (depending if the node is origin of this beam or not)  in the global dof vector
                    index_2_w_prime = global_index_2 + local_index_2 * 2 + 1  # location of w'(0) or w'(L) dof ................ in the global dof vector

                    constrain_vectror = np.zeros((6 * n_beams,
                                                  3))  # initialize constraint matrix (3 equations for every pair of beams: x,y compliance + stiff angle )
                    # first equation (x constraint)
                    constrain_vectror[index_1_v, 0] += cos_phi1
                    constrain_vectror[index_1_w, 0] += -sin_phi1
                    constrain_vectror[index_2_v, 0] += -cos_phi2
                    constrain_vectror[index_2_w, 0] += sin_phi2
                    # second equation(y constraint)
                    constrain_vectror[index_1_v, 1] += sin_phi1
                    constrain_vectror[index_1_w, 1] += cos_phi1
                    constrain_vectror[index_2_v, 1] += -sin_phi2
                    constrain_vectror[index_2_w, 1] += -cos_phi2
                    # third equation (stiff ange condition)
                    constrain_vectror[index_1_w_prime, 2] += 1
                    constrain_vectror[index_2_w_prime, 2] += -1
                    # print("FREE node ",constrain_vectror)
                    C = np.hstack([C, constrain_vectror])


            elif node.status == "FIXED":
                for beam in node.beams:
                    local_index = self.is_origin(node, beam)
                    global_index = beam.index * 6  # where dof of this vector are located in the global dof vector (every beam has 6 dof)
                    index_v = global_index + local_index - 1  # location of v(0) or v(L) dof in the global dof vector
                    index_w = global_index + local_index * 2  # location of w(0) or w(L) dof in the global dof vector
                    index_w_prime = global_index + local_index * 2 + 1  # location of w'(0) or w'(L) dof in the global dof vector
                    constrain_vectror = np.zeros((6 * n_beams,
                                                  3))  # initialize constraint matrix (3 equations for every beam: v, w fixed + w' fixed )
                    constrain_vectror[index_v, 0] += 1
                    constrain_vectror[index_w, 1] += 1
                    constrain_vectror[index_w_prime, 2] += 1
                    C = np.hstack([C, constrain_vectror])

        Se = np.vstack([np.hstack([S, C]), np.hstack([C.T, np.zeros((C.shape[1], C.shape[1]))])])

        Me = np.vstack([np.hstack([M_tilde, 0*C]), np.hstack([0*C.T, np.zeros((C.shape[1], C.shape[1]))])])
        self.Se_matrix = Se
        self.S_matrix = S
        self.C_matrix = C
        self.RHS = RHS
        self.Me_matrix = Me
        self.M_matrix  = M_tilde
        return M_tilde, Se, RHS

    def solve_system(self):  # when the matrices and RHS are constructed, steady solution can be obtained
        RHS2 = np.hstack([self.RHS, np.zeros(self.Se_matrix.shape[0] - len(self.RHS))])
        dof = np.linalg.solve(self.Se_matrix, RHS2)
        self.dof = dof
        return dof

    def solve_dynamic(self,h,t0,T):
        RHS2 = np.hstack([self.RHS, np.zeros(self.Se_matrix.shape[0] - len(self.RHS))])*0
        init = self.solve_system()
        init_tup = (init,0*init,0*init)
        sol, t = newmarkMethod(self.Me_matrix, self.Se_matrix, RHS2 , init_tup, h, t0, T)
        self.sol_dyn = sol
        self.t = t
        return sol, t

    def animate_frame(self, xlim = np.array([0, 10]), ylim = np.array([0, 10])):

            nB     = len(self.beams)  # Number of beams
            nN     = len(self.nodes)  # Number of nodes
            nNb    = 2                # Number of nodes per beam
            nDOFn  = 3                # Number of DOF per node
            nDOFb  = nNb * nDOFn      # Number of DOF per beam
            
            plt.rcParams['text.usetex'] = True
            plt.rcParams.update({'font.size' : 9})
            
            ### Plot definition
            self.fig = plt.Figure(figsize = (5, 3), dpi = 150)
            
            
            ### Window definition
            self.root = Tk.Tk()
            self.root.title('Animation test')
        
            label = Tk.Label(self.root, text = "Animation test").grid(column = 0, row = 0)
        
            canvas = FigureCanvasTkAgg(self.fig, master = self.root)
            canvas.get_tk_widget().grid(column=0,row=1)

            canvas = FigureCanvasTkAgg(self.fig, master = self.root)
            canvas.get_tk_widget().grid(column=0,row=1)
            
            ### Animation definition
            self.ax = self.fig.add_subplot(111)            
            
            def animation_frame(i):
                
                self.ax.clear()
                self.ax.set_ylim(ylim)
                self.ax.set_xlim(xlim)
                 
                for idxBeam in range(nB):
                
                    DOF = range(idxBeam * nDOFb, (idxBeam + 1) * nDOFb)
                    
                    sol_L = (self.sol_dyn[:,i])[DOF[0:2]]
                    sol_T = (self.sol_dyn[:,i])[DOF[2:]]
                    
                    grid = np.array([0, self.beams[idxBeam].length]) 
                    
                    v = interp1d(grid, sol_L * scaler) 
                    w = get_sol(grid, sol_T  * scaler)  
                    
                    x0    = self.beams[idxBeam].offset
                    
                    if abs(self.beams[idxBeam].direction[0]) > 1E-6:
                        theta = np.arctan(self.beams[idxBeam].direction[1]/self.beams[idxBeam].direction[0]) 
        
                    else:
                        theta = np.pi/2
        
                    if self.beams[idxBeam].direction[0] < 0:
                        theta = theta + np.pi
                        
                    if abs((abs(self.beams[idxBeam].direction[1]) - 1)) < 1E-6 and self.beams[idxBeam].direction[1] < 0:
                        theta = theta + np.pi
                    
                    self.beams[idxBeam].theta = theta    
                    
                    rotBeam  = rotateBeam((v, w), grid[-1], x0, theta)
                    original = rotateBeam((lambda x: x * 0, lambda x: x * 0), grid[-1], x0, theta)
                             
                    x_plot = np.linspace(0, self.beams[idxBeam].length, 20)
                    
                    self.ax.plot(rotBeam(x_plot)[0, :],  rotBeam(x_plot)[1, :],  color = 'k')
                    self.ax.plot(original(x_plot)[0, :], original(x_plot)[1, :], color = 'gray', linestyle = '--')

            self.ani = animation.FuncAnimation(self.fig, animation_frame, np.arange(0, 200), interval = 10, blit= False)
            self.root.mainloop()

    def eigen_freq_modes(self, Num, index, dynamic = False, t_0 = None, t_f = None, Nt = None, modes = None):

        if dynamic: 
            sol = eigenvalue_method_dynamic(t_0,t_f,Nt,self.Me_matrix,self.Se_matrix,modes,Num)
            self.sol_dyn = sol
            return sol
        else: 
            eigenvalues, eigenmodes = eigenvalue_method(self.Me_matrix,Num,self.Se_matrix)
            self.eigenvalues = eigenvalues
            self.eigenmodes = eigenmodes
            self.dof = eigenmodes[:,index-1]
            return eigenvalues, eigenmodes

class Node:
    def __init__(self, index, coord, status):
        self.index = index
        self.coordinates = coord  # x,y coordinates in numpy array
        self.status = status  # type of joint defined by a string(or character)
        self.beams = []  # reference to all "Beam objects" that meet at this node
        self.force = np.array([0.0, 0.0, 0.0])  # force that act on this node (if any)


class Beam:
    def __init__(self, index, nodes):
        self.index = index  # number of this beam
        self.nodes = nodes  # reference to a "Node objects" that are at the ends of this beam
        self.offset = nodes[0].coordinates  # coordinates of the origin of this beam in global ref frame
        self.direction = (nodes[1].coordinates - nodes[0].coordinates)
        self.length = np.sqrt(np.sum(self.direction ** 2))  # scalar length of the beam
        self.direction = self.direction / self.length  # unit vector in the direction of the beam in global frame








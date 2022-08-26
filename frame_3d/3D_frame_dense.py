import numpy as np
import matplotlib.pyplot as plt
import scipy as sp 
from scipy import sparse
from scipy.sparse import linalg
#np.set_printoptions(precision=  4)


class Structure:
    def __init__(self,filename):
        # copy all information about structure from a text file to an object 
        
        self.nodes = []
        self.beams = []
        self.C_matrix = None # constraints 
        self.S_matrix = None  # without constraints 
        self.Se_matrix = None # with constraints 
        self.M_matrix = None # to be done 
        self.RHS = None # 
        self.dof = 0 # solved degrees of freedom 
        self.E = 1.0
        self.I = 1.0
        self.A = 1.0
        self.J = 1.0
        self.G = 1.0
        self.mu = 1.0
        mode = 0
        f = open(filename, "r") # read off nodes and beams from the file and store them in lists 

        for line in f:
            
            if line[0] == "#":
                continue

            if line.strip() == "NODES":
                mode = 1
                node_index = 0
                continue

            if line.strip() == "EDGES":
                mode = 2
                beam_index = 0
                continue

            if line.strip() == "TYPE":
                mode = 3
                continue

            if line.strip() == "PARAMETERS":
                mode = 4
                continue

            if mode == 1:
                content = line.strip().split()
                coord = np.array([float(content[0]), float(content[1]), float(content[2])]) # x,y coordinates of a node in a global coordinate system
                #status = content[2] # get status of the node 
                node = Node(coord) # create node object
                node.index = node_index
                #if status == "FORCE": # if there is a force applied to the node, store it in the node object
                #    node.force = np.array([float(content[3]), float(content[4]),float(content[5])])
                self.nodes.append(node) # store this node in a list with all other nodes in the structure 
                node_index +=1
                print(node.coordinates)

            if mode == 2:
                content = line.strip().split()
                print(content)
                nodes = (self.nodes[int(content[0])-1], self.nodes[int(content[1])-1]) # 2 nodes that define a beam

                beam = Beam(nodes) # create beam object
                beam.index = beam_index # write index of the beam
                nodes[0].beams.append(beam) # add a reference to this beam for every node that is a part of this beam 
                nodes[1].beams.append(beam)
                print(beam.index)
                self.beams.append(beam) # store this beam in a list with all other beams in the structure 
                beam_index +=1

            if mode == 3:
                content = line.strip().split()
                node_index = int(content[0])-1
                node_status = content[1]
                self.nodes[node_index].status = node_status
                print(content)
                if len(content) > 2:
                    a = np.array([[content[2],content[3],content[4]],[content[5],content[6],content[7]]])
                    print(a)
                    self.nodes[node_index].force = a.astype(np.float)

            if mode == 4:
                content = line.strip().split()
                identifier = content[0]
                value = float(content[1])
                if identifier == "E":
                    self.E = value
                elif identifier == "I":
                    self.I = value
                elif identifier == "A":
                    self.A = value
                elif identifier == "mu":
                    self.mu = value
                elif identifier == "G":
                    self.G = value
                elif identifier == "J":
                    self.J = value
        f.close()
        for node in self.nodes:
            print(node.index, "  ",node.coordinates, " ", node.status," ",node.force," ")
        for beam in self.beams:
            print(beam.index, "  ",beam.offset, " ", )


    def plot(self):
        # draws a frame with all nodes, beams and force
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        


        for node in self.nodes :
            if node.status == "FREE" :
                color = 'g'
                marker= 'o'
            elif node.status == "FORCE" :
                color = 'y'
                marker= 'o'
            elif node.status == "FIXED" :
                color = 'r'
                marker= 's'

            x,y,z = node.coordinates
            index = node.index
            ax.scatter([x],[y],[z], marker=marker, color = color)
            
        
        for beam in self.beams :
            x1,y1,z1 = beam.nodes[0].coordinates # [x,y] of first node 
            x2,y2,z2 = beam.nodes[1].coordinates # [x,y] of second node
            v = beam.direction * 0.1 *beam.length # scaled beam direction
            index = beam.index # number of the beam
            ax.plot([x1,x2],[y1,y2],[z1,z2],color = 'orange')
        
        ax.legend()
        plt.show()


    def get_line(self,beam,n):
        dof = self.dof[beam.index*12:12*beam.index+12]
        offset = beam.offset
        # generates a 3d curve in global ref frame 
        u_0,u_L, v_0, v_L, w_0, wp_0, w_L, wp_L, l_0, lp_0, l_L, lp_L = dof
        #print("dof: ",dof)
        L = beam.length
        t = np.linspace(0.0,L,n)
        t_n = np.linspace(0.0,1.0,n) # normalized
        v = v_0 + t * (v_L-v_0)/L
        
        # basis functions at sampled points  
        phi1 = 1- 3 * t_n**2 + 2* t_n**3
        phi2 = t_n* (t_n-1)**2
        phi3 = 3*t_n**2 - 2*t_n**3
        phi4 = t_n**2 * (t_n-1)
        x = t + v
        y = l_0*phi1 + lp_0*phi2 + l_L*phi3 +  lp_L*phi4
        z = w_0*phi1 + wp_0*phi2 + w_L*phi3 +  wp_L*phi4
        angles = self.get_Euler_angles(beam)
        R , _ = self.get_transformations(angles)
        A = R @  np.vstack([x,y,z])
        A = A + np.reshape(offset,(1,3)).T
        return  A

    def plot_deformed(self):
        # fing the anchor point
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        for beam in self.beams:
            points = self.get_line(beam,10)
            x,y,z = (points[0],points[1],points[2])
            a = [beam.nodes[0].coordinates[0],beam.nodes[1].coordinates[0]]
            b = [beam.nodes[0].coordinates[1],beam.nodes[1].coordinates[1]]
            c = [beam.nodes[0].coordinates[2],beam.nodes[1].coordinates[2]]
            ax.plot(a,b,c,color="orange",)
            ax.plot(x,y,z,linestyle = "--")
            ax.plot(x[:2],y[:2],z[:2])
        plt.show()

    def is_origin(self,node,beam):
        if beam.nodes[0] == node:
            return 1
        elif beam.nodes[1] == node:
            return 2
        else :
            print("node doesnt belong to this beam")
            return 0

    def get_Euler_angles(self,beam):
        u = beam.direction
        roll = 0.0
        if np.isclose(u[0], 0.0, rtol=1e-06) and  np.isclose(u[1], 0.0, rtol=1e-06):
            yaw = 0.0
            if u[2] > 0.0:
                pitch = -np.pi /2
            else:
                pitch = np.pi /2
        else:
            yaw = np.arctan2(u[1],u[0])
            pitch = np.arcsin(-u[2])
        return np.array([yaw,pitch,roll])

    def get_transformations(self,angles):
        #computes transformation matrix from local to global  coordinate frame and vise-versa based on the direction of the beam 
        yaw,pitch,roll = angles[0],angles[1],angles[2]
  
        R_yaw = np.array([[np.cos(yaw),-np.sin(yaw), 0.0],
                          [np.sin(yaw), np.cos(yaw), 0.0],
                          [0.0        ,     0.0    , 1.0]])

        R_pitch = np.array([[np.cos(pitch) , 0.0, np.sin(pitch)],
                            [    0.0       , 1.0,        0.0   ],
                            [-np.sin(pitch), 0.0, np.cos(pitch)]])  

        R_roll =  np.array([[ 1.0 ,       0.0   ,      0.0     ],
                            [ 0.0 , np.cos(roll), -np.sin(roll)],
                            [ 0.0 , np.sin(roll),  np.cos(roll)]])  

                        
        R_g_l  = R_yaw @ R_pitch @  R_roll# from local to global
        R_l_g = np.linalg.inv(R_roll) @ np.linalg.inv(R_pitch) @ np.linalg.inv(R_yaw)
        return (R_g_l, R_l_g) # local-> global, global->local



    def get_perturbation_derivatives(self,beam1,beam2,vec):
        derivatives = []
        h = 10**-8 # step size for central difference 
        for i in range(6):
            a0 = np.zeros(6)
            a1 = np.zeros(6)
            a1[i] = h
            a2 = np.zeros(6)
            a2[i] = -h
            #print("a1:  ",a1)
            d = (self.angle_constraint_function(beam1,beam2,vec,a1) - self.angle_constraint_function(beam1,beam2,vec,a0))/(h)
            derivatives.append(d)
        return np.array(derivatives)

    def assemble_matrices(self):
             
        n_beams = len(self.beams) 
        S = np.zeros((n_beams*12, n_beams*12), dtype= np.float64)
        A_loc = np.array([[ 12. ,  6. ,-12. ,  6.],
                          [  6. ,  4. , -6. ,  2.],
                          [-12. , -6. , 12. , -6.],
                          [  6. ,  2. , -6. ,  4.]])
        mode1 = np.array([ [1,0,1,0],
                            [0,0,0,0],
                            [1,0,1,0],
                            [0,0,0,0] ])  

        mode2 = np.array([ [0,1,0,1],
                            [1,0,1,0],
                            [0,1,0,1],
                            [1,0,1,0] ]) 
    
        mode3 = np.array([ [0,0,0,0],
                            [0,1,0,1],
                            [0,0,0,0],
                            [0,1,0,1] ]) 

        #loop through beams and construct matrix S:
        for beam in self.beams:
            global_index = beam.index*12
            h = beam.length
            h_array = (mode1 * h**-3) +(mode2 * h**-2) +(mode3 * h**-1)
            S_torsion = np.array([[1.0,-1.0],[-1.0,1.0]])  * self.G * self.J /h
            S_v = np.array([[1.0,-1.0],[-1.0,1.0]])  * self.E * self.A /h
            S_w = self.E * self.I * A_loc * h_array
            S_l = self.E * self.I * A_loc * h_array
            S_local = sp.linalg.block_diag(S_torsion,S_v,S_w,S_l) # 12x12 local matrix corresponding to one element
            
            S[global_index:global_index+12,global_index:global_index+12] = S_local
            
        #loop through nodes and get constraints and forces
        
        RHS = np.zeros(n_beams*12)
        C = np.zeros((0,n_beams*12))

        for node in self.nodes:
            if node.status == "FREE" or node.status == "FORCE":
                
                anchor_beam = node.beams[0] # beam to which all other beams in this node will be pairwise attached
                print("DIRECTION:", anchor_beam.direction)
                Euler_angles1 = self.get_Euler_angles(anchor_beam)
                R_g_l1 , R_l_g1 = self.get_transformations(Euler_angles1) # direct and inverse coordinate transformation matrices
                local_index_1 = self.is_origin(node,anchor_beam) # 1 if origin, 2 if the other end 
                global_index_1 = anchor_beam.index*12 # where dof of this vector are located in the global dof vector (every beam has 6 dof)
                index_1_u = global_index_1 + local_index_1-1 # location of u(0) or u(L) dof(torsion)
                index_1_v = global_index_1 + local_index_1+1 # location of v(0) or v(L) dof (depending if the node is origin of this beam or not)  in the global dof vector 
                index_1_w = global_index_1 + local_index_1*2+2 # location of w(0) or w(L) dof 
                index_1_w_prime =  global_index_1 + local_index_1*2 +3# location of w'(0) or w'(L) dof 
                index_1_l = global_index_1 + local_index_1*2+6 # location of l(0) or l(L) dof 
                index_1_l_prime =  global_index_1 + local_index_1*2 +7# location of l'(0) or l'(L) dof
                
                if node.status == "FORCE":
                    f = R_l_g1 @ node.force[0] # compute forces and moments in local coordinate frame 
                    M = R_l_g1 @ node.force[1]
                    RHS[index_1_v] += f[0]# I don't know if this is correct !!!!!! 
                    RHS[index_1_l] += f[1]# 
                    RHS[index_1_w] += f[2]# I don't know if this is correct !!!!!! 

                    RHS[index_1_u] += M[0]#  moment around x axis 
                    RHS[index_1_w_prime] += M[1]# moment around y axis 
                    RHS[index_1_l_prime] += M[2]# moment around z axis 
                
                for beam in node.beams[1:]:
                    
                    Euler_angles2 = self.get_Euler_angles(beam)
                    R_g_l2 , R_l_g2 = self.get_transformations(Euler_angles2) # direct and inverse coordinate transformation matrices
                    local_index_2 = self.is_origin(node,beam) # 1 if origin, 2 if the other end 
                    global_index_2 = beam.index*12 # where dof of this vector are located in the global dof vector (every beam has 6 dof)
                    index_2_u = global_index_2 + local_index_2-1 # location of u(0) or u(L) dof(torsion)
                    index_2_v = global_index_2 + local_index_2+1 # location of v(0) or v(L) dof (depending if the node is origin of this beam or not)  in the global dof vector 
                    index_2_w = global_index_2 + local_index_2*2+2 # location of w(0) or w(L) dof 
                    index_2_w_prime =  global_index_2 + local_index_2*2 +3# location of w'(0) or w'(L) dof 
                    index_2_l = global_index_2 + local_index_2*2+6 # location of l(0) or l(L) dof 
                    index_2_l_prime =  global_index_2 + local_index_2*2 +7# location of l'(0) or l'(L) dof
                    
                    constrain_vectror = np.zeros((3,12*n_beams), dtype= np.float64) # initialize constraint matrix (3 equations for every pair of beams: x,y compliance + stiff angle )
                    # position constraint
                    constrain_vectror[:,index_1_v] = R_g_l1[:,0]
                    constrain_vectror[:,index_1_l] = R_g_l1[:,1]
                    constrain_vectror[:,index_1_w] = R_g_l1[:,2]

                    constrain_vectror[:,index_2_v] = -R_g_l2[:,0]
                    constrain_vectror[:,index_2_l] = -R_g_l2[:,1]
                    constrain_vectror[:,index_2_w] = -R_g_l2[:,2]

                    # angle constraint 
                    angle_constraint = np.zeros((3,12*n_beams), dtype= np.float64)
                   
                    angle_constraint[:,index_1_u]        =  R_g_l1[:,0]
                    angle_constraint[:,index_1_w_prime]  = -R_g_l1[:,1]
                    angle_constraint[:,index_1_l_prime]  =  R_g_l1[:,2]

                    angle_constraint[:,index_2_u]        = -R_g_l2[:,0]
                    angle_constraint[:,index_2_w_prime]  =  R_g_l2[:,1]
                    angle_constraint[:,index_2_l_prime]  = -R_g_l2[:,2]
                     
                    C = np.vstack([C,constrain_vectror,angle_constraint])#
                    

            elif node.status == "FIXED":
                for beam in node.beams:
                    c_vec = np.zeros((6,12*n_beams), dtype= np.float64)
                    local_index = self.is_origin(node,beam) # 1 if origin, 2 if the other end 
                    global_index = beam.index*12 # where dof of this vector are located in the global dof vector (every beam has 6 dof)
                    index_u = global_index + local_index-1 # location of u(0) or u(L) dof(torsion)
                    index_v = global_index + local_index+1 # location of v(0) or v(L) dof (depending if the node is origin of this beam or not)  in the global dof vector 
                    index_w = global_index + local_index*2+2 # location of w(0) or w(L) dof 
                    index_w_prime =  global_index + local_index*2 +3# location of w'(0) or w'(L) dof 
                    index_l = global_index + local_index*2+6 # location of l(0) or l(L) dof 
                    index_l_prime =  global_index + local_index*2 +7# location of l'(0) or l'(L) dof
                    c_vec[0,index_v] = 1.0
                    c_vec[1,index_w] = 1.0
                    c_vec[2,index_l] = 1.0
                    c_vec[3,index_w_prime] = 1.0
                    c_vec[4,index_l_prime] = 1.0
                    c_vec[5,index_u] = 1.0
                    C = np.vstack([C,c_vec])

        print("C size",C.shape)
        print("C rank: ",np.linalg.matrix_rank(C))
        nc = C.shape[0]

        Se = np.hstack([np.vstack([S,C]),np.vstack([C.transpose(),np.zeros((nc,nc))])])
        print("Se size",Se.shape)
        print("Se rank: ",np.linalg.matrix_rank(Se))

        print("S size",S.shape)
        print("S rank: ",np.linalg.matrix_rank(S))

        self.Se_matrix = Se
        self.S_matrix = S
        self.C_matrix = C
        self.RHS = RHS
        return Se, RHS

    def solve_system(self): # when the matrices and RHS are constructed, steady solution can be obtained 
        RHS2 = np.hstack([self.RHS,np.zeros(self.Se_matrix.shape[0]-len(self.RHS))])
        dof = np.linalg.solve(self.Se_matrix,RHS2)
        self.dof = dof
        return dof

    

class Node:
    def __init__(self,coord):
        self.index = None
        self.coordinates = coord # x,y coordinates in numpy array
        self.status = None # type of joint defined by a string(or character)
        self.beams = [] # reference to all "Beam objects" that meet at this node 
        self.force = np.array([[0.0,0.0,0.0],[0.0,0.0,0.0]]) # force that act on this node (if any)

class Beam:
    def __init__(self,nodes):
        self.index = None # number of this beam 
        self.nodes = nodes # reference to a "Node objects" that are at the ends of this beam
        self.offset = nodes[0].coordinates # coordinates of the origin of this beam in global ref frame
        self.direction = (nodes[1].coordinates - nodes[0].coordinates) 
        self.length = np.sqrt(np.sum(self.direction**2)) # scalar length of the beam
        self.direction = self.direction / self.length # unit vector in the direction of the beam in global frame
        self.angles = None



#filename = "test.txt"
filename = r'C:\Users\Volodymyr\Documents\Master_courses\project_numerical_analysis\beam-num-analysis\src\frame\test_3d.txt'
x = Structure(filename)

x.assemble_matrices()

x.solve_system()
print(x.dof)
#x.get_line(0,x.beams[0],10)
x.plot()
x.plot_deformed()
#print(x.dof) # vector of u and v as in script page 14 

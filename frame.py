
import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(precision=2)
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
        self.A = 3.0
        mode = 0
        f = open(filename, "r") # read off nodes and beams from the file and store them in lists 

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
                coord = np.array([float(content[0]), float(content[1])]) # x,y coordinates of a node in a global coordinate system
                status = content[2] # get status of the node 
                node = Node(node_index,coord,status) # create node object
                if status == "FORCE": # if there is a force applied to the node, store it in the node object
                    node.force = np.array([float(content[3]), float(content[4]),float(content[5])])
                self.nodes.append(node) # store this node in a list with all other nodes in the structure 
                node_index +=1
                print(node.coordinates,node.status,node.force)

            if mode == 2:
                content = line.strip().split()
                nodes = (self.nodes[int(content[0])], self.nodes[int(content[1])]) # 2 nodes that define a beam
                beam = Beam(beam_index,nodes) # create beam object
                nodes[0].beams.append(beam) # add a reference to this beam for every node that is a part of this beam 
                nodes[1].beams.append(beam)
                print(beam.index)
                self.beams.append(beam) # store this beam in a list with all other beams in the structure 
                beam_index +=1
        f.close()
    def plot(self):
        # draws a frame with all nodes, beams and force 
        plt.axes()
        for node in self.nodes :
            p = node.coordinates
            index = node.index
            circle = plt.Circle(p, radius=0.05, fc='r')
            text = plt.text(p[0]-0.02, p[1]-0.02, str(index),fontsize = 10,color = 'w')
            plt.gca().add_patch(circle)

        for beam in self.beams :
            p1 = beam.nodes[0].coordinates # [x,y] of first node 
            p2 = beam.nodes[1].coordinates # [x,y] of second node
            v = beam.direction * 0.1 *beam.length # scaled beam direction
            index = beam.index # number of the beam
            line = plt.Line2D((p1[0],p2[0]), (p1[1],p2[1]), lw=1)
            arrow = plt.arrow(p1[0], p1[1], v[0], v[1],head_width=0.02, color = 'r')
            text = plt.text((p1[0]+p2[0])/2, (p1[1]+p2[1])/2,str(index),backgroundcolor = 'w',fontsize = 7,color = 'b')
            plt.gca().add_line(line)
            plt.gca().add_patch(arrow)
            
       
        plt.axis('scaled')
        plt.show()

    def is_origin(self,node,beam):
        if beam.nodes[0] == node:
            return 1
        elif beam.nodes[1] == node:
            return 2
        else :
            print("node doesnt belong to this beam")
            return 0
        
    def assemble_matrices(self):
        n_beams = len(self.beams)
        S = np.zeros((n_beams*6,n_beams*6)) # global stiffness matrix 
        S1_loc = np.array([[1.0,-1.0],[-1.0,1.0]])
        S2_loc = np.array([[ 12. ,  6. ,-12. ,  6.],
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
            global_index = beam.index*6
            h = beam.length
            h_array = (mode1 * h**-3) +(mode2 * h**-2) +(mode3 * h**-1)  
            S[global_index:global_index+2,global_index:global_index+2] = self.E * self.A * S1_loc/h
            S[global_index+2:global_index+6,global_index+2:global_index+6] = self.E * self.I * S2_loc * h_array

        #loop through nodes and get constraints and forces
        RHS = np.zeros(n_beams*6)
        C = np.zeros((n_beams*6,0))
        
        for node in self.nodes:
            if node.status == "FREE" or node.status == "FORCE":
                
                anchor_beam = node.beams[0] # beam to which all other beams in this node will be pairwise attached
                cos_phi1 = anchor_beam.direction[0]
                sin_phi1 = anchor_beam.direction[1]
                local_index_1 = self.is_origin(node,anchor_beam)
                global_index_1 = anchor_beam.index*6 # where dof of this vector are located in the global dof vector (every beam has 6 dof)
                index_1_v = global_index_1 + local_index_1-1 # location of v(0) or v(L) dof (depending if the node is origin of this beam or not)  in the global dof vector 
                index_1_w = global_index_1 + local_index_1*2 # location of w(0) or w(L) dof (depending if the node is origin of this beam or not)  in the global dof vector
                index_1_w_prime =  global_index_1 + local_index_1*2 +1# location of w'(0) or w'(L) dof ................ in the global dof vector
                
                if node.status == "FORCE":
                    fx = node.force[0]
                    fy = node.force[1]
                    M = node.force[2]
                    RHS[index_1_v] += fx*cos_phi1 + fy*sin_phi1# I don't know if this is correct !!!!!! 
                    RHS[index_1_w] += -fx*sin_phi1 + fy*cos_phi1# I don't know if this is correct !!!!!! 
                    RHS[index_1_w_prime] += M# reserved for the moment
                
                for beam in node.beams[1:]:
                    
                    cos_phi2 = beam.direction[0]
                    sin_phi2 = beam.direction[1]
                    local_index_2 = self.is_origin(node,beam)
                    global_index_2 = beam.index*6 # where dof of this vector are located in the global dof vector (every beam has 6 dof)
                    index_2_v = global_index_2 + local_index_2-1 # location of v(0) or v(L) dof (depending if the node is origin of this beam or not)  in the global dof vector 
                    index_2_w = global_index_2 + local_index_2*2 # location of w(0) or w(L) dof (depending if the node is origin of this beam or not)  in the global dof vector
                    index_2_w_prime =  global_index_2 + local_index_2*2 +1# location of w'(0) or w'(L) dof ................ in the global dof vector
                    
                    constrain_vectror = np.zeros((6*n_beams,3)) # initialize constraint matrix (3 equations for every pair of beams: x,y compliance + stiff angle )
                    # first equation (x constraint)
                    constrain_vectror[index_1_v,0] += cos_phi1
                    constrain_vectror[index_1_w,0] += -sin_phi1
                    constrain_vectror[index_2_v,0] += -cos_phi2
                    constrain_vectror[index_2_w,0] += sin_phi2
                    # second equation(y constraint)
                    constrain_vectror[index_1_v,1] += sin_phi1
                    constrain_vectror[index_1_w,1] += cos_phi1
                    constrain_vectror[index_2_v,1] += -sin_phi2
                    constrain_vectror[index_2_w,1] += -cos_phi2
                    # third equation (stiff ange condition)
                    constrain_vectror[index_1_w_prime,2] += 1
                    constrain_vectror[index_2_w_prime,2] += -1
                    #print("FREE node ",constrain_vectror)
                    C = np.hstack([C,constrain_vectror])


            elif node.status == "FIXED":
                for beam in node.beams:
                    local_index = self.is_origin(node,beam)
                    global_index = beam.index*6 # where dof of this vector are located in the global dof vector (every beam has 6 dof)
                    index_v = global_index + local_index-1 # location of v(0) or v(L) dof in the global dof vector 
                    index_w = global_index + local_index*2 # location of w(0) or w(L) dof in the global dof vector
                    index_w_prime =  global_index + local_index*2 +1# location of w'(0) or w'(L) dof in the global dof vector
                    constrain_vectror = np.zeros((6*n_beams,3)) # initialize constraint matrix (3 equations for every beam: v, w fixed + w' fixed )
                    constrain_vectror[index_v,0] += 1
                    constrain_vectror[index_w,1] += 1
                    constrain_vectror[index_w_prime,2] += 1
                    C = np.hstack([C,constrain_vectror])

        Se = np.vstack([np.hstack([S,C]),np.hstack([C.T,np.zeros((C.shape[1],C.shape[1]))])])
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
    def __init__(self,index,coord,status):
        self.index = index
        self.coordinates = coord # x,y coordinates in numpy array
        self.status = status # type of joint defined by a string(or character)
        self.beams = [] # reference to all "Beam objects" that meet at this node 
        self.force = np.array([0.0,0.0,0.0]) # force that act on this node (if any)

class Beam:
    def __init__(self,index,nodes):
        self.index = index # number of this beam 
        self.nodes = nodes # reference to a "Node objects" that are at the ends of this beam
        self.offset = nodes[0].coordinates # coordinates of the origin of this beam in global ref frame 
        self.direction = (nodes[1].coordinates - nodes[0].coordinates) 
        self.length = np.sqrt(np.sum(self.direction**2)) # scalar length of the beam
        self.direction = self.direction / self.length # unit vector in the direction of the beam in global frame 




filename = "frame1.txt"

x = Structure(filename)

x.assemble_matrices()

x.solve_system()

print(x.dof) # vector of u and v as in script page 14 

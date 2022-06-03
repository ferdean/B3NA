
import numpy as np
import matplotlib.pyplot as plt

class Structure:
    def __init__(self,filename):
        # copy all information about structure from a text file to an object 
        self.nodes = []
        self.beams = []
        self.C_matric = None # constraints 
        self.S_matrix = None 
        self.M_matrix = None # to be done 
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
                    node.force = np.array([float(content[3]), float(content[4])])
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
    def assemble_matrices(self):
        #loop through nodes and get constraints and forces
        n_beams = len(self.beams)
        C = np.zeros((n_beams,0))
        for node in self.nodes:
            print(C)
        return 0

class Node:
    def __init__(self,index,coord,status):
        self.index = index
        self.coordinates = coord # x,y coordinates in numpy array
        self.status = status # type of joint defined by a string(or character)
        self.beams = [] # reference to all "Beam objects" that meet at this node 
        self.force = np.array([0.0,0.0]) # force that act on this node (if any)

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
#x.plot()
x.assemble_matrices()

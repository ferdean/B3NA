import numpy as np

class Structure:
    def __init__(self,filename):
        # copy all information about structure from a text file to an object 
        self.nodes = []
        self.beams = []
        mode = 0
        f = open(filename, "r") # read off nodes and beams from the file and store them in lists 

        for line in f:
            
            if line[0] == "#":
                continue

            if line.strip() == "NODES:":
                mode = 1
                continue

            if line.strip() == "BEAMS:":
                mode = 2
                continue

            if mode == 1:
                content = line.strip().split()
                coord = np.array([float(content[0]), float(content[1])]) # x,y coordinates of a node in a global coordinate system
                status = content[2] # get status of the node 
                node = Node(coord,status) # create node object
                if status == "FORCE": # if there is a force applied to the node, store it in the node object
                    node.forces = np.array([float(content[3]), float(content[4])])
                self.nodes.append(node) # store this node in a list with all other nodes in the structure 
                print(node.coordinates,node.status,node.forces)

            if mode == 2:
                content = line.strip().split()
                nodes = (int(content[0]), int(content[1])) # 2 nodes that define a beam
                beam = Beam(nodes) # create beam object
                self.nodes[nodes[0]].beams.append(beam) # add a reference to this beam for every node that is a part of this beam 
                self.nodes[nodes[1]].beams.append(beam)

                self.beams.append(beam) # store this beam in a list with all other beams in the structure 

        f.close()




class Node:
    def __init__(self,coord,status):
        self.coordinates = coord
        self.status = status
        self.beams = []
        self.forces = np.array([0.0,0.0])

class Beam:
    def __init__(self,nodes):
        self.nodes = nodes
        self.offset = nodes[0]
        self.direction = (nodes[1] - nodes[0])
        self.direction = self.direction / np.sqrt(np.sum(self.direction**2))

filename = "frame1.txt"

x = Structure(filename)


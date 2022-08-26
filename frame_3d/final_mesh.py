# script to turn .stl files into the standard BENA .txt file

import numpy as np
from stl import mesh
import itertools

#
#
#
# sellect stl file 
#stl_data = mesh.Mesh.from_file('radar_v1.stl')
stl_data = mesh.Mesh.from_file('real_radar_v1.stl')

def close_to(p1,p2,tol):
    if np.sum((p1-p2)**2) < tol:
        return True
    else:
        return False

p = stl_data.points
v = stl_data.vectors
print(v[1])
print(p[1])
print("qdwdqwd")


point_list = v.reshape((v.shape[0]*v.shape[1],3))
#print(v.shape)
#print(v)
point_list = np.unique(point_list, axis=0)


print(point_list)
print(point_list.shape)
#print(p)

beam_list1 = []


for triangle in v:
    tol = 10**-10
    p1 = triangle[0]
    p2 = triangle[1]
    p3 = triangle[2]
    for p,i in zip(point_list,range(len(point_list))):
        if close_to(p,p1,tol):
            id1 = i+1
        if close_to(p,p2,tol):
            id2 = i+1
        if close_to(p,p3,tol):
            id3 = i+1
    beam_list1.append([id1,id2])
    beam_list1.append([id2,id3])
    beam_list1.append([id3,id1])

beam_list = []

for i in range(len(beam_list1)):
    lst = beam_list1[i].copy()
    lst.sort()
    beam_list.append(lst)

beam_list.sort()
beam_list = list(beam_list for beam_list,_ in itertools.groupby(beam_list))
print(len(beam_list))
# generate file     
#
#
#
##
# sellect output file
f = open("radar.txt", "w")
f.write("# Final structure\n")
f.write("NODES\n")
for i in range(len(point_list)):
    f.write(str(i+1) +" "+ np.array2string(point_list[i])[1:-1] + "\n")

f.write("EDGES\n")

for i in beam_list:
    f.write(str(i[0]) + " " + str(i[1]) + "\n")

f.write("TYPE\n")

for i in range(len(point_list)):
    f.write(str(i+1) + " " + "FREE" + "\n")

f.write("PARAMETERS\n")
f.write("E 12.0\n")
f.write("I 22.0\n")
f.write("G 11.0\n")
f.write("J 42.0\n")
f.write("A 15.0\n")
f.write("mu 12.0\n")

f.close()
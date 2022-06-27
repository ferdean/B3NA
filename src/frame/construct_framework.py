import numpy as np
import matplotlib.pyplot as plt
import sys 

options = sys.argv[1:]
print(
    'Options:\n'
    ' \'whole\': Round clicked points to nearest integers.\n'
    ' \'0-index\': Uses 0-indexing for the points'
)
indexing = 0 if '0-index' in options else 1

################ Get points ################
fig = plt.figure()
ax = fig.add_subplot(111)
ax.title.set_text('Add nodes')
ax.plot(0,0)
ax.set(xlim=(0,10),ylim=(0,10))
ax.grid()
coords = []

def onclick(event):
    ix, iy = event.xdata, event.ydata

    if ix is None or iy is None:
        return
    if 'whole' in options:
        ix = round(ix); iy = round(iy)
    if (ix,iy) in coords:
        return
    coords.append((ix, iy))

    ax.plot(ix,iy,marker="o",color="black")
    ax.set(xlim=(0,10),ylim=(0,10))
    fig.canvas.draw()

cid = fig.canvas.mpl_connect('button_press_event', onclick)

plt.show()
if len(coords) == 0:
    exit()
############################################

################ Get edges #################
def find_closest(coords, point):
    min_length = None
    index = None
    for i,coord in enumerate(coords):
        length = np.linalg.norm(np.array(coord)-np.array(point))
        if min_length is None or length < min_length:
            min_length = length
            index = i
    return index

fig = plt.figure()
ax = fig.add_subplot(111)
ax.title.set_text('Draw edges')
x = [x[0] for x in coords]
y = [y[1] for y in coords]
ax.plot(x,y,"o",color="black")
ax.set(xlim=(0,10),ylim=(0,10))
for i,(x,y) in enumerate(coords):
    ax.annotate(i+indexing, (x+0.1,y+0.1))

edges = []
line_coords = []
line_indices = []
node_indices = []

def onclick(event):
    ix, iy = event.xdata, event.ydata
    try:
        node_index = find_closest(coords, (ix,iy))
        node_indices.append(node_index + indexing)
    except TypeError:
        return
    line_coords.append(coords[node_index])
    if len(line_coords) == 2:
        if sorted([x for x in node_indices]) not in [sorted(x) for x in edges]:
            ax.plot([line_coords[0][0],line_coords[1][0]],[line_coords[0][1],line_coords[1][1]],color="black")
            ax.set(xlim=(0,10),ylim=(0,10))
            edges.append((node_indices[0],node_indices[1]))
        
        line_coords.pop()
        line_coords.pop()
        node_indices.pop()
        node_indices.pop()

    fig.canvas.draw()

cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show(block=False)
############################################

############# Determine types ##############
print('Determine type of each node:\n 1 -> FIXED\n 2 -> FREE\n 3 -> FORCE\n 4 -> MOVABLE')
allowed_vals = ['1','2','3','4']
typedict = {'1':'FIXED','2':'FREE','3':'FORCE','4':'MOVABLE'}
typesstr = []
types = []
for i in range(indexing,len(coords)+indexing):
    nr = ''
    while nr not in allowed_vals:
        nr = input(f'Node {i}: ')
    if nr == '3':
        Fx = input('F_x = ')
        Fy = input('F_y = ')
        M = input('M = ')
        typesstr.append(f'{i} FORCE [{Fx} {Fy}] [{M}]\n')
    elif nr == '4':
        x = input('x-val = ')
        y = input('y-val = ')
        typesstr.append(f'{i} MOVABLE [{x} {y}]\n')
    else:
        typesstr.append(f'{i} {typedict[nr]}\n')
    types.append(typedict[str(nr)])
############################################


########### Determine parameters ###########
parameters = ''
print('Determine parameter functions for E,I,A and mu')
pars = ['E','I','A','mu']
for par in pars:
    fnc = input(f'{par}(x) = ')
    parameters += f'{par} {fnc}\n'
parameters = parameters[0:-1]
############################################

############### Save file ##################
save_file = input('Save to file in ../frame (.ne will be added): ')

with open(f'{save_file}.ne', 'w') as f:
    f.write('NODES\n')
    for coord in coords:
        f.write(f'{coord[0]} {coord[1]}\n')
    f.write('EDGES\n')
    for edge in edges:
        f.write(f'{edge[0]} {edge[1]}\n')
    f.write('TYPE\n')
    for type_ in typesstr:
        f.write(type_)
    f.write('PARAMETERS\n')
    f.write(parameters)
############################################

############ Plot with types ###############
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set(xlim=(0,10),ylim=(0,10))
ax.grid()

colors = {'FIXED':'black','FREE':'blue','FORCE':'red','MOVABLE':'green'}
fig.text(0.2,0.9,'FIXED',ha='center',va='bottom',size='large',color=colors['FIXED'])
fig.text(0.4,0.9,'FREE',ha='center',va='bottom',size='large',color=colors['FREE'])
fig.text(0.6,0.9,'FORCE',ha='center',va='bottom',size='large',color=colors['FORCE'])
fig.text(0.8,0.9,'MOVABLE',ha='center',va='bottom',size='large',color=colors['MOVABLE'])
for edge in edges:
    x1,y1 = coords[edge[0]-indexing]
    x2,y2 = coords[edge[1]-indexing]

    ax.plot([x1,x2],[y1,y2],color='black')
for (i,coord,type) in zip(range(len(coords)),coords,types):
    ax.plot(coord[0],coord[1],color=colors[type],marker='o')
    ax.annotate(i+indexing, (coord[0]+0.1,coord[1]+0.1),color='black')
plt.savefig(f'../frame/{save_file}.png')
############################################
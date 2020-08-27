import numpy as np
import matplotlib.pyplot as plt
import PyQt5
from scipy.interpolate import UnivariateSpline, CubicSpline
from scipy.integrate import quad
# from pp_utils import *

'''
Assumed inputs:
    1. IMS track description
    2. Current vehicle position in Frenet space (s,l)
    3. Array of obstacles -> just circles for now, defined as ((x,y), radius) tuples

This code does the following:
    1. Build Frenet-space lattice & generate splines
    2. Assign cost to each edge, using a curvature-based cost function
    3. Given current position and a fixed planning horizon, extract local segment of graph for search
    4. Collision checks for obstacles
        a.  Node deletion for larger obstacles (also deletes associated edges to/from node)
        b.  Edge deletion for smaller (static) obstacles
    5. Shortest path search based on edge costs
    6. C2 Trajectory generation

TO-DO:
    1. Implement a graph data structure to hold properties, coordinates, and coefficients at each node
    2. Fix initial headings on graph generation. All should not be the same at the same station
    3. Separate graph gen (offline) code from online portion
    4. Optimize collision checks; separate into path-like and static obstacles
    5. Consider different type of graph search (DFS, A*)
'''

# Load Track Data (meters)
track_data = np.genfromtxt("../data/ims_track_midline.csv", delimiter=",", names=True)

# Fit spline
s_map = np.linspace(0, 1.0, num=track_data.shape[0])
x_spline = UnivariateSpline(s_map, track_data["x"])
dxds     = x_spline.derivative()
y_spline = UnivariateSpline(s_map, track_data["y"])
dyds     = y_spline.derivative()

# Plot Track
x_tmp = np.array([x_i for i, x_i in enumerate(track_data["x"]) if i % 10 == 0])
y_tmp = np.array([y_i for i, y_i in enumerate(track_data["y"]) if i % 10 == 0])
s_tmp = np.linspace(0, 1.0, num=track_data.shape[0]*5)
plt.title('IMS Track (meters)')
plt.plot(x_tmp, y_tmp, "o", x_spline(s_map), y_spline(s_map), "-")
plt.axes().set_aspect('equal')
plt.show()

# Estimate length of track midline
# Should be close to 2.5 miles * 1609.34 = 4023.35 meters
n_interps = 50000
t_lg = np.linspace(0.0, 1.0, n_interps)
spl_coords = np.zeros((n_interps, 2))
spl_coords[:,0] = x_spline(t_lg)
spl_coords[:,1] = y_spline(t_lg)
tracklength = np.sum(np.sqrt(np.sum(np.power(np.diff(spl_coords, axis=0), 2), axis=1)))
print('Track Length: ', tracklength)

# Re-fit spline using estimated track length. Allows frenet station to correspond
# to actual distance along the raceline/midline 
s_map = np.linspace(0, tracklength, num=track_data.shape[0])
x_spline = UnivariateSpline(s_map, track_data["x"])
dxds     = x_spline.derivative()
y_spline = UnivariateSpline(s_map, track_data["y"])
dyds     = y_spline.derivative()

# generate frenet space
def p(s, l, x_ref, y_ref):
    dx = dxds(s)
    dy = dyds(s)
    dx = dx / (dx**2 + dy**2)**(0.5)
    dy = dy / (dx**2 + dy**2)**(0.5)
    n_x = -dy 
    n_y =  dx 
    return np.vstack((n_x*l + x_ref(s), n_y*l + y_ref(s)))

def t(s, l, x_ref, y_ref):
    dx = dxds(s)
    dy = dyds(s)
    dx = dx / (dx**2 + dy**2)**(0.5)
    dy = dy / (dx**2 + dy**2)**(0.5)
    return np.vstack((dx*l + x_ref(s), dy*l + y_ref(s)))

def heading(s, dxds, dyds):
    dx = dxds(s)
    dy = dyds(s)
    return np.arctan2(dy, dx)

##### Spline helpers

def spline(u, a0, a1, a2, a3):
    return a3*u**3 + a2*u**2 + a1*u + a0

def calc_curvature(x_coeff, y_coeff, n_interps=50):
    # Assumes s normalized (0,1)
    t = np.linspace(0.0, 1.0, n_interps)
    x_d = x_coeff[1] + 2 * x_coeff[2] * t + 3 * x_coeff[3] * np.power(t, 2)
    y_d = y_coeff[1] + 2 * y_coeff[2] * t + 3 * y_coeff[3] * np.power(t, 2)
    x_dd = 2 * x_coeff[2] + 6 * x_coeff[3] * t
    y_dd = 2 * y_coeff[2] + 6 * y_coeff[3] * t
    kappa = (x_d * y_dd - y_d * x_dd) / np.power(np.power(x_d, 2) + np.power(y_d, 2), 1.5)
    return kappa

def calc_spline_length(x_coeff, y_coeff, n_interps=50):
    t = np.linspace(0.0, 1.0, n_interps)
    x = x_coeff[0] + x_coeff[1] * t + x_coeff[2] * np.power(t,2) + x_coeff[3] * np.power(t,3)
    y = y_coeff[0] + y_coeff[1] * t + y_coeff[2] * np.power(t,2) + y_coeff[3] * np.power(t,3)
    spl_coords = np.hstack((x.reshape(n_interps,1), y.reshape(n_interps,1)))
    spl_len = np.sum(np.sqrt(np.sum(np.power(np.diff(spl_coords, axis=0), 2), axis=1)))
    return spl_len

def edge_cost_curv(layers, coeffs, lidx_end):

    spl_len = calc_spline_length(coeffs[0], coeffs[1])
    k = calc_curvature(coeffs[0], coeffs[1])
    k_avg = np.mean(np.abs(k))
    k_delta = np.max(k) - np.min(k)
    lat_offset = abs(lidx_end)
    # Using TUM weights for now
    w_len = 10
    w_kavg = 750
    w_kdelta = 1500
    w_offset = 5000
    cost = spl_len * (w_len + w_kavg*k_avg**2 + w_kdelta*k_delta**2 + w_offset*lat_offset)
    return cost


#### Build grid space
s_map = np.array([d for d in range(0, int(tracklength), 20)])
l_map = np.array([d for d in range(-10, 11)])
displacements = l_map

layers = []
for i in range(len(s_map)):
    station = s_map[i]
    layers.append(dict())

    for l in displacements:
        
        if -track_data["width_left"][i] < l < track_data["width_right"][i]:
            xy = p(station,l, x_spline, y_spline)
            layers[i][l] = (heading(station, dxds, dyds), xy[0][0], xy[1][0])  

#### connect layers
paths = []
for i in range(len(s_map) - 1):
    s_start = s_map[i]
    s_end   = s_map[i+1]

    paths.append(dict())

    for j, l_start in enumerate(layers[i]):

        paths[i][l_start] = dict()

        for l_end in layers[i+1]:

            heading_start, x_start, y_start = layers[i]  [l_start]
            heading_end,   x_end,   y_end   = layers[i+1][l_end]

            x_a0 = x_start
            x_a1 = np.cos(heading_start) # zeroth order length approximation
            x_a2 = 3*x_end - np.cos(heading_end) - 2*x_a1 - 3*x_a0
            x_a3 = x_end - x_a2 - x_a1 - x_a0

            y_a0 = y_start
            y_a1 = np.sin(heading_start) # zeroth order length approximation
            y_a2 = 3*y_end - np.sin(heading_end) - 2*y_a1 - 3*y_a0
            y_a3 = y_end - y_a2 - y_a1 - y_a0

            coeffs = (
                (x_a0, x_a1, x_a2, x_a3), 
                (y_a0, y_a1, y_a2, y_a3)
                )  

            ## Cost info
            cost = edge_cost_curv(layers, coeffs, l_end)
            
            paths[i][l_start][l_end] = (
                (x_a0, x_a1, x_a2, x_a3), 
                (y_a0, y_a1, y_a2, y_a3),
                cost
                )  

            

##### Local Graph Extraction
def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx

def nearest_node_idx(s_map, l_map, s, l):
    s_idx = find_nearest(s_map, s)
    s_idx = s_idx + 1 if (s_map[s_idx] - s) < 0 else s_idx
    l_idx = l_map[find_nearest(l_map, l)]
    return s_idx, l_idx

def extract_local_graph(layers, paths, plan_horizon, cur_s, cur_l):
    ''' 
    Given full graph (layers + paths) and current location in frenets space,
    extract the local graph, returning local layers and paths
    '''
    start_idxs = nearest_node_idx(s_map, l_map, cur_s, cur_l)
    end_idxs = nearest_node_idx(s_map, l_map, cur_s + plan_horizon, 0)
    # Extract graph between start and end idxs
    local_layers = layers[start_idxs[0]:end_idxs[0]+1]
    local_paths = paths[start_idxs[0]:end_idxs[0]]

    # Replace first layer with only the closest (start) node
    local_layers[0] = {start_idxs[1]: local_layers[0][start_idxs[1]]}

    # Remove all paths originating from old nodes in the first layer
    local_paths[0] = {start_idxs[1]: local_paths[0][start_idxs[1]]}
    return local_layers, local_paths


##### Obstacle Avoidance

def collision_check_circle(pt, obs):
    '''
    Checks for collisions by determining whether point pt is inside circle obs
    pt is a 2-tuple of form (x,y)
    Obs consists of a 2-tuple ((x,y),r) to define the circle
    Returns False for no collision; True indicates collision
    '''
    d = ((pt[0] - obs[0][0])**2 + (pt[1] - obs[0][1])**2)**0.5
    if d <= obs[1]:
        return True
    else:
        return False

def remove_node(layers, paths, s_idx, l_idx):
    # Delete node from layers
    if l_idx in layers[s_idx]:
        del layers[s_idx][l_idx]
    
    # Delete outgoing paths
    if l_idx in paths[s_idx]:
        del paths[s_idx][l_idx]
    
    # Delete incoming paths
    for snode, path in paths[s_idx - 1].items():
        if l_idx in path:
            del path[l_idx]
            
def remove_edge(paths, sidx_start, lidx_start, lidx_end):
    # Find edge in paths and delete
    if lidx_start in paths[sidx_start] and lidx_end in paths[sidx_start][lidx_start]:
        del paths[sidx_start][lidx_start][lidx_end]
        
        
def collision_checker(layers, paths, obs_list):
    delete_nodes = []
    delete_edges = []
    u = np.linspace(0,1,20)  # Last term determines # of sample points for each spline
    for s,layer in enumerate(layers):
        for node in layer:
            x,y = layer[node][1], layer[node][2]
            for obs in obs_list:
                if collision_check_circle((x,y), obs):  # True indicates collision
                    delete_nodes.append((s, node))
                    
    for s,layer in enumerate(paths):
        for startnode in layer:
            for dn in layer[startnode]:
                coeff = layer[startnode][dn]
                x_spline = spline(u, *coeff[0])
                y_spline = spline(u, *coeff[1])
                for i in range(len(x_spline)):
                    pt = (x_spline[i], y_spline[i])
                    for obs in obs_list:
                        if collision_check_circle(pt, obs):
                            delete_edges.append((s,startnode,dn))
    # Delete everything in delete lists
    for i in delete_nodes:
        remove_node(layers, paths, *i)
        
    for i in delete_edges:
        remove_edge(paths, *i)   

    # check for orphaned nodes with no parents or no children
    for s,layer in enumerate(paths):
        for node in layer:
            # Check if node has no children
            if layer[node] == {}:
                delete_nodes.append((s,node))
            else:
                # Check if node has parents
                if s != 0:
                    no_parents = True
                    for parentnode in paths[s-1]:
                        if node in paths[s-1][parentnode]:
                            no_parents = False
                            break
                    if no_parents:
                        delete_nodes.append((s,node))

    # Delete everything in delete lists
    for i in delete_nodes:
        remove_node(layers, paths, *i)
        
    for i in delete_edges:
        remove_edge(paths, *i)   

def plot_graph(ax, layers, paths):
    u = np.linspace(0,1)
    for layer in layers:
        for node in layer:
            x,y = layer[node][1], layer[node][2]
            ax.scatter(x,y, color='b')
    for layer in paths:
        for startnode in layer:
            for dn in layer[startnode]:
                coeff = layer[startnode][dn]
                x_spline = spline(u, *coeff[0])
                y_spline = spline(u, *coeff[1])
                ax.plot(x_spline, y_spline, color='skyblue')
    ax.set_aspect('equal')

def plot_obstacles(ax, obs_list):
    for obs in obs_list:
        obs_shape = plt.Circle(obs[0], obs[1], color = 'black', zorder=10)
        ax.add_artist(obs_shape)


##### Path Search

def find_min_cost_path(local_layers, local_paths):
    '''
    DP approach working backwards from end of local graph to find min cost path
    '''
    n_stations = len(local_layers)
    stations = range(n_stations - 1, -1, -1)
    segments = range(n_stations - 2, -1, -1)

    costs = {}  # mimics layers structure
    for s in stations:
        costs[s] = {}
    for l_start in local_layers[n_stations-1]:   # init last layer costs to use lateral offset
        costs[n_stations-1][l_start] = (abs(l_start), l_start)

    for s in segments:
        for l_start in local_paths[s]:
            # if s == n_stations - 1:   # initialize last layer costs
            #     costs[s][l_start] = 0  
            # else:
                
            for l_end in local_paths[s][l_start]:  # child nodes

                # Get edge cost
                ec = local_paths[s][l_start][l_end][2] 
                
                # Update parent node cost
                if l_start in costs[s]:
                    # update cost
                    if (ec + costs[s+1][l_end][0]) < costs[s][l_start][0]:
                        costs[s][l_start] = (ec + costs[s+1][l_end][0], l_end)
                else:
                    costs[s][l_start] = (ec + costs[s+1][l_end][0], l_end)  # initialize

    min_cost_path = []  # (station(s), node(l))
    # for s in stations:
    for s in range(n_stations):
        if s == 0:
            cur_node = list(costs[s].keys())[0]
            child_node = list(costs[s].values())[0][1]
        else:
            child_node = costs[s][cur_node][1]

        min_cost_path.append((s, cur_node))
        cur_node = child_node
    
    return min_cost_path


def plot_mcp(ax, layers, paths, mcp):
    u = np.linspace(0,1)
    # Plot nodes
    for i,p in enumerate(mcp):
        data = layers[p[0]][p[1]]
        x,y = data[1], data[2]
        ax.scatter(x, y, color = 'orange', zorder = 15)
    
    # Plot edges
    for i in range(len(mcp)-1):
        p = mcp[i+1]
        s = mcp[i][0]
        l = mcp[i][1]
        coeffs = paths[s][l][p[1]]
        x_spline = spline(u, *coeffs[0])
        y_spline = spline(u, *coeffs[1])
        ax.plot(x_spline, y_spline, color='purple')
        
##### Extract curvature-continuous path
def mcp_to_coords(layers, mcp):
    '''
    Takes in mcp and returns np arrays x and y for cartesian coordinates + start/end headings
    '''
    x = np.zeros((len(mcp), 1))
    y = np.zeros((len(mcp), 1))
    for sn in mcp:  # format (s, node)
        s = sn[0]
        n = sn[1]
        x[s] = layers[s][n][1]
        y[s] = layers[s][n][2]

    # Boundary conditions - headings
    bch_start = layers[0][mcp[0][1]][0]
    bch_end = layers[mcp[-1][0]][mcp[-1][1]][0]
    bch = (bch_start, bch_end)

    return x, y, bch


def gen_c2_spline(x, y, bc_headings, slen_start, slen_end):
    '''
    Generates a C2 continuous spline using scipy CubicSpline lib
    x: np.array of x-coordinate points
    y: np.array of y-coordinate points
    '''
    # define mu, a virtual path variable of length 1 for each spline segment
    assert(len(x) == len(y))
    mu = np.arange(0,len(x), 1.0)

    # build splines
    cs_x = CubicSpline(mu, x, 
                   bc_type=((1, slen_start * np.cos(bc_headings[0])), 
                            (1, slen_end * np.cos(bc_headings[1]))))
    cs_y = CubicSpline(mu, y, 
                   bc_type=((1, slen_start * np.sin(bc_headings[0])), 
                            (1, slen_end * np.sin(bc_headings[1]))))
    return cs_x, cs_y


def calc_c2_traj(x, y, bc_headings, eps = 0.005):
    '''
    Iteratively compute spline coefficients until spline length of first and last segment converges
    '''
    # Start with euclidean dist as slen approx for first and last segments
    slen_start = np.sqrt((x[1] - x[0])**2 + (y[1] - y[0])**2)
    slen_end = np.sqrt((x[-1] - x[-2])**2 + (y[-1] - y[-2])**2)
#     print(slen_start, slen_end)

    while True:
        cx, cy = gen_c2_spline(x, y, bc_headings, slen_start, slen_end)
        coeffs_x_start = np.flip(cx.c[:,0])
        coeffs_y_start = np.flip(cy.c[:,0])
        coeffs_x_end = np.flip(cx.c[:,-1])
        coeffs_y_end = np.flip(cy.c[:,-1])

        slen_start_new = calc_spline_length(coeffs_x_start, coeffs_y_start)
        slen_end_new = calc_spline_length(coeffs_x_end, coeffs_y_end)
#         print(slen_start_new, slen_end_new)

        if abs(slen_start_new - slen_start) < eps and abs(slen_end_new - slen_end) < eps:
            break
        else:
            slen_start = slen_start_new
            slen_end = slen_end_new
    return cx, cy


def plot_trajectory(ax, cx, cy, n_segments, n_interps = 50):
    t = np.linspace(0.0, n_segments, n_segments*n_interps)
    ax.plot(cx(t), cy(t), color='olive')


###################### Test Case ##################
# Current location of vehicle
cur_s = 1010.85
cur_l = 3.1
horizon = 200  # given along s midline (meters)
# Create some obstacles
# Small static - circle
obstacle1 = ((617.8, -493.1), 1.5)   # ((x,y), rad) -> circle
# path-like - series of circles
obstacle2 = ((647, -476), 4.2)
obstacle3 = ((673.0, -436.0), 5.0)
obs_list = [obstacle1, obstacle2, obstacle3]

# Extract local graph & plot
local_layers, local_paths = extract_local_graph(layers, paths, horizon, cur_s, cur_l)

### Plot raw extraction graph
fig, ax = plt.subplots(figsize=(10,6))
plot_graph(ax, local_layers, local_paths)
plt.scatter(*p(cur_s, cur_l, x_spline, y_spline), color='green')  # print current location
plt.show()

# Obstacle checking
collision_checker(local_layers, local_paths, obs_list)
fig, ax = plt.subplots(figsize=(10,6))
plot_graph(ax, local_layers, local_paths)
plt.scatter(*p(cur_s, cur_l, x_spline, y_spline), color='green')  # print current location
plot_obstacles(ax, obs_list)

# Min cost path
mcp = find_min_cost_path(local_layers, local_paths)
print(mcp)
plot_mcp(ax, local_layers, local_paths, mcp)
# plt.show()

# C2 trajectory generation
x, y, bch = mcp_to_coords(local_layers, mcp)
traj = calc_c2_traj(x, y, bch)
print(x)
plot_trajectory(ax, traj[0], traj[1], len(mcp)-1)
plt.show()





import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.interpolate import UnivariateSpline, CubicSpline
from utils import *

"""
Script to Generate Track Lattice and Cost Map

Inputs
    - trackfile: [x ,y, width_left, width_right] csv file describing track midline and boundary widths
    - parameters: See descriptions below in comments; set within this script

Outputs
    - nodes.json: output file in ./data with x,y,heading info for each node in track lattice
    - edges.json: output file in ./data with C1 spline coefficients and curvature-based cost value for each edge

Info

"""

##### PARAMETERS #####
trackfile             = "./data/ims_track_midline.csv"  # Data in meters. Defined as midline and normal distance to track boundaries
parsed_trackfile_out  = "./data/ims_track_parsed.csv"   # Location to write parsed track data
vehicle_width         = 1.6                             # Meters
lateral_safety_buffer = 1.0                             # Min distance from vehicle outer edge to track bound for node creation (meters)

long_sep_straight     = 25                              # Longitudinal separation between nodes on straights (m)
long_sep_curve        = 10                              # Longitudinal separation between nodes on curves (m)
lat_sep               = 1.0                             # Lateral separation between nodes (m)
curvature_thresh      = 0.001                           # If segment's max curvature > threshold, segment is considered a curve, else a straight

cost_weights = {                                        # Weights for curvature cost function that sets costs on each segment length
    'w_len': 10,
    'w_kavg': 750,
    'w_kdelta': 1500,
    'w_offset': 5000
}

##### PARSE TRACK GEOMETRY #####

# Load track data and close the loop w/ start point
track_data = np.genfromtxt(trackfile, delimiter=",", names=True)
track_data = np.append(track_data, track_data[0])
# TODO if/else append or replace last point depending on distance from orginal last to start

# Fit midline spline using provided discretization
s_map        = np.linspace(0, 1.0, num=track_data.shape[0])
x_spline_mid = CubicSpline(s_map, track_data["x"])
y_spline_mid = CubicSpline(s_map, track_data["y"])

# Estimate length of track midline
# Should be close to 2.5 miles * 1609.34 = 4023.35 meters
n_interps    = 50000
t_fine       = np.linspace(0.0, 1.0, n_interps)
spl_coords   = np.zeros((n_interps, 2))
spl_coords[:,0] = x_spline_mid(t_fine)
spl_coords[:,1] = y_spline_mid(t_fine)
tracklength     = np.sum(np.sqrt(np.sum(np.power(np.diff(spl_coords, axis=0), 2), axis=1)))
print('Track Length: ', tracklength)

# Re-fit spline using estimated track length. Allows frenet station to correspond
# to actual distance along the raceline/midline 
s_map        = np.linspace(0, tracklength, num=track_data.shape[0])
x_spline_mid = CubicSpline(s_map, track_data["x"])
y_spline_mid = CubicSpline(s_map, track_data["y"])

# Remap to finer meter-level interpolation
s_fine = s_fine = np.arange(0, tracklength, 1)
x_mid  = x_spline_mid(s_fine)
y_mid  = y_spline_mid(s_fine)

# Get normal and tangent unit vectors wrt refline at every station
mid_norms    = get_unit_norms(s_fine, x_spline_mid, y_spline_mid)
mid_tangents = get_unit_tangents(s_fine, x_spline_mid, y_spline_mid)

# Headings and curvature
mid_headings = get_headings(s_fine, x_spline_mid, y_spline_mid)
mid_curvature = get_curvature(s_fine, x_spline_mid, y_spline_mid)

# Interpolate track widths
left_widths  = np.interp(s_fine, s_map, track_data["width_left"])
right_widths = np.interp(s_fine, s_map, track_data["width_right"])

# Get track boundary points (normal to midline)
left_boundary_xy  = frenet2cart(s_fine, -left_widths, x_spline_mid, y_spline_mid)
right_boundary_xy = frenet2cart(s_fine, right_widths, x_spline_mid, y_spline_mid)

# Get headings at track boundary by fitting boundary points to spline
left_bound_spline_x  = CubicSpline(s_fine, left_boundary_xy[:,0])
left_bound_spline_y  = CubicSpline(s_fine, left_boundary_xy[:,1])
right_bound_spline_x = CubicSpline(s_fine, right_boundary_xy[:,0])
right_bound_spline_y = CubicSpline(s_fine, right_boundary_xy[:,1])
left_bound_headings  = get_headings(s_fine, left_bound_spline_x, left_bound_spline_y)
right_bound_headings = get_headings(s_fine, right_bound_spline_x, right_bound_spline_y)

# Write parsed track data out to file
parsed_track = np.column_stack((s_fine, x_mid, y_mid, mid_headings, mid_curvature,
                                mid_norms, mid_tangents, 
                                left_widths, right_widths,
                                left_boundary_xy, left_bound_headings,
                                right_boundary_xy, right_bound_headings))
col_names = "s, x_mid, y_mid, psi_mid, kappa_mid x_n, y_n, x_t, y_t, \
             w_left, w_right, x_left, y_left, psi_left, x_right, y_right, psi_right"
np.savetxt(parsed_trackfile_out, parsed_track, delimiter=",", header=col_names)


##### GENERATE TRACK LATTICE #####

# Generate station intervals for nodes based on curvature parameters
stations = [0]                                       # start at s = 0
while s_fine[-1] - stations[-1] > long_sep_curve:
    s_cur = stations[-1]
    s_next_straight = s_cur + long_sep_straight
    s_next_curve = s_cur + long_sep_curve
    
    curved = False
    for s in range(s_cur, s_next_straight):
        if mid_curvature[s] >= curvature_thresh:
            curved = True
    if curved:
        stations.append(s_next_curve)
    else:
        stations.append(s_next_straight)
    
if tracklength - stations[-1] < long_sep_curve/2:
    stations.pop()
stations_cl = stations + [stations[0]]

# Generate nodes dict and with x,y,psi values for each node
# Example format: nodes[s][l][x/y/psi]: nodes['240']['-3.5']['y']
lat_disp = np.arange(-10,11,lat_sep)                 # max possible range of lateral offsets (meters)
nodes = {}
for s in stations:
    nodes[str(s)] = {}

    for l in lat_disp:
        l_bound = -left_widths[s] + (vehicle_width/2 + lateral_safety_buffer)
        r_bound = right_widths[s] - (vehicle_width/2 + lateral_safety_buffer)
        if l_bound < l < r_bound:
            
            xy = frenet2cart(s, l, x_spline_mid, y_spline_mid)
            psi = np.interp(l,                      # Linearly interpolate between boundary and midline headings
                            [-left_widths[s], 0, right_widths[s]], 
                            [left_bound_headings[s], mid_headings[s], right_bound_headings[s]]
                           )
            
            nodes[str(s)][str(l)] = {
                'x':     xy[0][0],
                'y':     xy[0][1],
                'psi':   psi
            }

# Save nodes dict to file
with open('./data/nodes.json', 'w') as fp:
    json.dump(nodes, fp,  indent=4)

# Generate edges dict with C1 spline coefficients and costs for each pair of adjacent stations
edges = {}
for i in range(len(stations_cl) - 1):
    s_start = stations_cl[i]
    s_end = stations_cl[i+1]
    edges[str(s_start)] = {}
    
    for l_start, v_start in nodes[str(s_start)].items():
        edges[str(s_start)][str(l_start)] = {}
        x_start   = v_start['x']
        y_start   = v_start['y']
        psi_start = v_start['psi']
        
        for l_end, v_end in nodes[str(s_end)].items():
            x_end   = v_end['x']
            y_end   = v_end['y']
            psi_end = v_end['psi']
            
            x_a0 = x_start
            x_a1 = np.cos(psi_start)
            x_a2 = 3*x_end - np.cos(psi_end) - 2*x_a1 - 3*x_a0
            x_a3 = x_end - x_a2 - x_a1 - x_a0

            y_a0 = y_start
            y_a1 = np.sin(psi_start) 
            y_a2 = 3*y_end - np.sin(psi_end) - 2*y_a1 - 3*y_a0
            y_a3 = y_end - y_a2 - y_a1 - y_a0
            
            coeffs = ((x_a0, x_a1, x_a2, x_a3),
                      (y_a0, y_a1, y_a2, y_a3))
            
            ## Cost info
            cost = edge_cost_curv(coeffs, float(l_end), **cost_weights)
            
            edges[str(s_start)][str(l_start)][str(l_end)] = {
                'x_coef': coeffs[0],
                'y_coef': coeffs[1],
                'cost': cost
            }

# Save edges dict to file
with open('./data/edges.json', 'w') as fp:
    json.dump(edges, fp,  indent=4)


#### Plot Track ####
plot_trackfile(track_data, x_mid, y_mid, left_boundary_xy, right_boundary_xy)
plot_lattice_nodes(x_mid, y_mid, left_boundary_xy, right_boundary_xy, nodes)
plot_lattice_paths_partial(x_mid, y_mid, left_boundary_xy, right_boundary_xy, 
                           nodes, edges, stations, 100, 220)
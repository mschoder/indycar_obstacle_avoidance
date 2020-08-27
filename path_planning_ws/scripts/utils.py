import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline, CubicSpline



def frenet2cart(s, l, x_ref, y_ref):
    '''
    Inputs:
        - s: station (float, meters)
        - l: lateral offset (float, meters)
        - x_ref, y_ref: x and y splines defining the frenet reference path (spline object)
    Outputs:
        - np array of (x,y) tuples in cartesian coordinates
    '''
    norms = get_unit_norms(s, x_ref, y_ref)
    n_x = norms[:,0]
    n_y = norms[:,1]
    return np.column_stack((n_x*l + x_ref(s), n_y*l + y_ref(s)))


def get_unit_norms(s, x_ref, y_ref):
    dx = x_ref.derivative()(s)
    dy = y_ref.derivative()(s)
    dx = dx / (dx**2 + dy**2)**(0.5)
    dy = dy / (dx**2 + dy**2)**(0.5)
    n_x = dy 
    n_y = -dx
    return np.column_stack((n_x, n_y))

def get_unit_tangents(s, x_ref, y_ref):
    dx = x_ref.derivative()(s)
    dy = y_ref.derivative()(s)
    dx = dx / (dx**2 + dy**2)**(0.5)
    dy = dy / (dx**2 + dy**2)**(0.5)
    return np.column_stack((dx, dy))

def get_headings(s, x_ref, y_ref):
    dx = x_ref.derivative()(s)
    dy = y_ref.derivative()(s)
    return np.arctan2(dy, dx)

def get_curvature(s, x_ref, y_ref):
    dx  = x_ref.derivative(1)(s)
    dy  = y_ref.derivative(1)(s)
    d2x = x_ref.derivative(2)(s)
    d2y = y_ref.derivative(2)(s)
    kappa = (dx * d2y - dy * d2x) / np.power(np.power(dx, 2) + np.power(dy, 2), 1.5)
    return kappa


#### Cost function ####
def curvature_from_coeffs(x_coeff, y_coeff, n_interps=50):
    # Assumes s normalized (0,1)
    t = np.linspace(0.0, 1.0, n_interps)
    x_d = x_coeff[1] + 2 * x_coeff[2] * t + 3 * x_coeff[3] * np.power(t, 2)
    y_d = y_coeff[1] + 2 * y_coeff[2] * t + 3 * y_coeff[3] * np.power(t, 2)
    x_dd = 2 * x_coeff[2] + 6 * x_coeff[3] * t
    y_dd = 2 * y_coeff[2] + 6 * y_coeff[3] * t
    kappa = (x_d * y_dd - y_d * x_dd) / np.power(np.power(x_d, 2) + np.power(y_d, 2), 1.5)
    return kappa

def spline_length_from_coeffs(x_coeff, y_coeff, n_interps=50):
    t = np.linspace(0.0, 1.0, n_interps)
    x = x_coeff[0] + x_coeff[1] * t + x_coeff[2] * np.power(t,2) + x_coeff[3] * np.power(t,3)
    y = y_coeff[0] + y_coeff[1] * t + y_coeff[2] * np.power(t,2) + y_coeff[3] * np.power(t,3)
    spl_coords = np.hstack((x.reshape(n_interps,1), y.reshape(n_interps,1)))
    spl_len = np.sum(np.sqrt(np.sum(np.power(np.diff(spl_coords, axis=0), 2), axis=1)))
    return spl_len

def edge_cost_curv(coeffs, lidx_end, w_len, w_kavg, w_kdelta, w_offset):
    spl_len = spline_length_from_coeffs(coeffs[0], coeffs[1])
    k = curvature_from_coeffs(coeffs[0], coeffs[1])
    k_avg = np.mean(np.abs(k))
    k_delta = np.max(k) - np.min(k)
    lat_offset = abs(lidx_end)
    cost = spl_len * (w_len + w_kavg*k_avg**2 + w_kdelta*k_delta**2 + w_offset*lat_offset)
    return cost


def spline_eval(u, a0, a1, a2, a3):
    return a3*u**3 + a2*u**2 + a1*u + a0

##### PLOT UTILS #####

def plot_trackfile(track_data, x_mid, y_mid, left_boundary_xy, right_boundary_xy):
    x_tmp = np.array([x_i for i, x_i in enumerate(track_data["x"]) if i % 10 == 0])
    y_tmp = np.array([y_i for i, y_i in enumerate(track_data["y"]) if i % 10 == 0])
    plt.title('IMS Track (meters)')
    plt.plot(x_tmp, y_tmp, "o", x_mid, y_mid, "-")
    plt.plot(left_boundary_xy[:,0], left_boundary_xy[:,1], color="blue")
    plt.plot(right_boundary_xy[:,0], right_boundary_xy[:,1], color="purple")
    plt.axes().set_aspect('equal')
    plt.show()


def plot_lattice_nodes(x_mid, y_mid, left_boundary_xy, right_boundary_xy, nodes):
    plt.figure(figsize=(10,16))
    plt.title('IMS Track Nodes')
    plt.plot(x_mid, y_mid, "-", color='orange')
    plt.plot(left_boundary_xy[:,0], left_boundary_xy[:,1], color="black")
    plt.plot(right_boundary_xy[:,0], right_boundary_xy[:,1], color="black")
    node_x = [nodes[s][l]['x'] for s in nodes for l in nodes[s]]
    node_y = [nodes[s][l]['y'] for s in nodes for l in nodes[s]]
    plt.plot(node_x, node_y, "o")
    plt.axes().set_aspect('equal')
    plt.show()

def plot_lattice_paths_partial(x_mid, y_mid, left_boundary_xy, right_boundary_xy, 
                               nodes, edges, stations, s_start, s_end):
    s_part = [s for s in stations if s_start <= s <= s_end]

    plt.figure(figsize=(10,16))
    plt.title('IMS Track Nodes')
    plt.plot(x_mid, y_mid, "-", color='orange')
    plt.plot(left_boundary_xy[:,0], left_boundary_xy[:,1], color="black")
    plt.plot(right_boundary_xy[:,0], right_boundary_xy[:,1], color="black")

    node_x = [nodes[s][l]['x'] for s in nodes for l in nodes[s] if int(s) in s_part]
    node_y = [nodes[s][l]['y'] for s in nodes for l in nodes[s] if int(s) in s_part]
    plt.plot(node_x, node_y, "o", color="darkblue")

    u = np.linspace(0,1)
    coefs_x = [edges[s][ls][le]['x_coef'] for s in edges for ls in edges[s] for le in edges[s][ls] if int(s) in s_part[:-1]]
    coefs_y = [edges[s][ls][le]['y_coef'] for s in edges for ls in edges[s] for le in edges[s][ls] if int(s) in s_part[:-1]]
    edges_x = [spline_eval(u, *cx) for cx in coefs_x]
    edges_y = [spline_eval(u, *cy) for cy in coefs_y]
    for i in range(len(edges_x)):
        plt.plot(edges_x[i], edges_y[i], color="skyblue")
    plt.axes().set_aspect('equal')
    plt.show()

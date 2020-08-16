#include <vector>
#include <math.h>
#include <array>
#include <Trackgraph.h>
#include <Obstacle.h>
#include <matplotlibcpp.h>

#include "interpolation.h"
#include "ap.h"
#include "ap.cpp"
#include "alglibinternal.cpp"
#include "alglibmisc.cpp"
#include "integration.cpp"
#include "interpolation.cpp"
#include "linalg.cpp"
#include "optimization.cpp"
#include "solvers.cpp"
#include "specialfunctions.cpp"


using namespace std;
namespace plt = matplotlibcpp;

namespace pp_utils {

vector<double> splineEval(const vector<double> &u, const array<double, 4> &coefs) {
    vector<double> eval;
    for (auto &el : u) {
        eval.push_back(coefs[3]*pow(el,3) + coefs[2]*pow(el,2) + coefs[1]*el + coefs[0]);
    }
    return eval;
}

vector<double> LinearSpacedArray(double a, double b, std::size_t N) {
    double h = (b - a) / static_cast<double>(N-1);
    vector<double> xs(N);
    vector<double>::iterator x;
    double val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) {
        *x = val;
    }
    return xs;
}

pair<vector<double>, vector<double>> circle_points(StaticObstacle &obs) {
    pair<vector<double>, vector<double>> result;
    for (double angle=0; angle<=2*M_PI; angle+=0.1) {
        result.first.push_back(obs.x + obs.rad*cos(angle));
        result.second.push_back(obs.y + obs.rad*sin(angle));
    }
    return result;
}

bool collision_check_circle(const Node &node, const StaticObstacle &obs) {
    // Return true if collision; false otherwise
    return obs.rad >= sqrt(pow(node.x - obs.x, 2) + pow(node.y - obs.y, 2)) ? true : false;
}

bool collision_check_circle(const double x, const double y, const StaticObstacle &obs) {
    return obs.rad >= sqrt(pow(x - obs.x, 2) + pow(y - obs.y, 2)) ? true : false;
}

void collision_checker(Trackgraph &graph, Obstacles &obstacles) {
    // TODO: implement dynamic obstacles checking
    // Nodes
    for (const auto &s : graph.nodes) {
        for (const auto &[offset,node] : s.second) {
            for (const auto &obs : obstacles.sv) {
                if (collision_check_circle(node, obs)) {
                    graph.removeNode(s.first, offset);
                }
            }
        }
    }
    // Edges
    vector<double> u = LinearSpacedArray(0, 1, 50); // for spline eval
    bool break_two = false;
    vector<int> s_vals;
    vector<double> ls_vals, le_vals;
    s_vals = graph.getEdgeStations();
    for (auto &s : s_vals) {
        ls_vals = graph.getEdgeOffsets(s);
        for (auto &ls : ls_vals) {
            le_vals = graph.getEdgeDests(s, ls);
            for (auto &le : le_vals) {
                auto &edge = graph.getEdge(s, ls, le);
                vector<double> xspl = splineEval(u, edge.x_coef);
                vector<double> yspl = splineEval(u, edge.y_coef);
                // check each spline pt against each obs pt
                for (int i = 0; i != xspl.size(); ++i) {
                    for (auto &obs : obstacles.sv) {
                        if (collision_check_circle(xspl[i], yspl[i], obs)) {
                            graph.removeEdge(s, ls, le);
                            break_two = true;
                            break;
                        }
                    }
                    if (break_two) {break_two = false; break;}
                }
            }
        }
    }
    // Check for nodes with no children
    vector<double> lsp_vals; // parent start offsets at s-1
    bool removed = true;
    while (removed){
        removed = false;
        s_vals = graph.getEdgeStations();
        for (auto &s : s_vals) {
            ls_vals = graph.getEdgeOffsets(s);
            for (auto &ls : ls_vals) {
                if (graph.getEdgeDests(s, ls).empty()) { // if no destination edges for a given parent
                    graph.removeNode(s, ls);
                    removed = true;
                // And check for orphan nodes
                } else if (s != 0) { // First node will have no parents - OK
                    bool orphan = true;
                    lsp_vals = graph.getEdgeOffsets(s-1);
                    for (auto &lsp : lsp_vals) {
                        if (graph.edges.at(s-1).at(lsp).find(ls) != graph.edges.at(s-1).at(lsp).end()) {
                            orphan = false;
                            break;
                        }
                    }
                    if (orphan) {
                        graph.removeNode(s, ls);
                        removed = true;
                    }
                }
            }
        }
    }
}

void plotTrack(Trackgraph &graph) {
    // Plot the global_graph
    vector<double> x_nodes, y_nodes;
    for (auto &s : graph.getNodeStations()) {
        for (auto &el : graph.getNodeOffsets(s)) {
            auto node = graph.getNode(s, el);
            x_nodes.push_back(node.x);
            y_nodes.push_back(node.y); 
        }
    }
    plt::figure_size(780, 1200);
    plt::plot(x_nodes, y_nodes, "o");
    plt::show();
}

void plotLocalGraph(Trackgraph &graph, Obstacles &obstacles, 
        vector<pair<int,double>> &mcp, vector<double> &pltx, vector<double> &plty) {
    // pltx, plty are overloads for trajectory
    // Nodes
    plt::figure_size(780, 1200);
    vector<double> x_nodes, y_nodes;
    for (auto &s : graph.getNodeStations()) {
        for (auto &el : graph.getNodeOffsets(s)) {
            auto node = graph.getNode(s, el);
            x_nodes.push_back(node.x);
            y_nodes.push_back(node.y); 
        }
    }
    // Edges
    vector<double> u = LinearSpacedArray(0, 1, 50);
    for (const auto &s : graph.getEdgeStations()) {
        for (const auto &sel : graph.getNodeOffsets(s)) {
            for (const auto &eel : graph.getEdgeDests(s, sel)) {
                const auto &edge = graph.getEdge(s, sel, eel);
                vector<double> x_paths = splineEval(u, edge.x_coef);
                vector<double> y_paths = splineEval(u, edge.y_coef);
                plt::plot(x_paths, y_paths, "skyblue");
            }
        }
    }
    // Obstacles
    for (auto &el : obstacles.sv) {
        pair<vector<double>, vector<double>> obs_pts = circle_points(el);
        plt::plot(obs_pts.first, obs_pts.second, "k-");
    }
    plt::plot(x_nodes, y_nodes, "bo");
    // Min Cost Path
    pair<vector<double>, vector<double>> mcp_xy = mcp_to_xy(graph, mcp);
    plt::plot(mcp_xy.first, mcp_xy.second, "ro");
    // Trajectory
    plt::plot(pltx, plty);
    plt::show();
}

void plotLocalGraph(Trackgraph &graph, Obstacles &obstacles) {
    // Nodes
    plt::figure_size(780, 1200);
    vector<double> x_nodes, y_nodes;
    for (auto &s : graph.getNodeStations()) {
        for (auto &el : graph.getNodeOffsets(s)) {
            auto node = graph.getNode(s, el);
            x_nodes.push_back(node.x);
            y_nodes.push_back(node.y); 
        }
    }
    // Edges
    vector<double> u = LinearSpacedArray(0, 1, 50);
    for (const auto &s : graph.getEdgeStations()) {
        for (const auto &sel : graph.getNodeOffsets(s)) {
            for (const auto &eel : graph.getEdgeDests(s, sel)) {
                const auto &edge = graph.getEdge(s, sel, eel);
                vector<double> x_paths = splineEval(u, edge.x_coef);
                vector<double> y_paths = splineEval(u, edge.y_coef);
                plt::plot(x_paths, y_paths, "skyblue");
            }
        }
    }
    // Obstacles
    for (auto &el : obstacles.sv) {
        pair<vector<double>, vector<double>> obs_pts = circle_points(el);
        plt::plot(obs_pts.first, obs_pts.second, "k-");
    }
    plt::plot(x_nodes, y_nodes, "bo");
    plt::show();
}

// Trajectory generation -- TODO - format output as non-alglib, depending on msg needed
pair<alglib::spline1dinterpolant, alglib::spline1dinterpolant> trajectory_gen(
    Trackgraph &graph, vector<pair<int, double>> &mcp) {
    
    pair<vector<double>, vector<double>> mcp_xy = mcp_to_xy(graph, mcp);
    int n_interps = 50;
    int n_points = mcp_xy.first.size();
    vector<double> svec = LinearSpacedArray(0, 1, n_points);
    int sfirst = mcp.front().first;
    int slast = mcp.back().first;
    double lfirst = mcp.front().second;
    double llast = mcp.back().second;
    double heading_start = graph.nodes.at(sfirst).at(lfirst).psi;
    double heading_last = graph.nodes.at(slast).at(llast).psi;
    double hix = cos(heading_start);
    double hex = cos(heading_last);
    double hiy = sin(heading_start);
    double hey = sin(heading_last);

    alglib::real_1d_array s, x, y;
    s.setlength(n_points);
    x.setlength(n_points);
    y.setlength(n_points);
    for (int i = 0; i <= n_points-1; ++i) {
        s(i) = svec.at(i);
        x(i) = mcp_xy.first.at(i);
        y(i) = mcp_xy.second.at(i);
    }

    // initialize slen_start and slen_end
    double x0, x1, xlast, x2last, y0, y1, ylast, y2last;
    x0 = mcp_xy.first.at(0);
    x1 = mcp_xy.first.at(1);
    xlast = mcp_xy.first.at(n_points-1);
    x2last = mcp_xy.first.at(n_points-2);
    y0 = mcp_xy.second.at(0);
    y1 = mcp_xy.second.at(1);
    ylast = mcp_xy.second.at(n_points-1);
    y2last = mcp_xy.second.at(n_points-2);
    double slen_start = sqrt(pow((x1 - x0), 2) + pow((y1 - y0), 2));
    double slen_end = sqrt(pow((xlast - x2last), 2) + pow((ylast - y2last), 2));

    // Fit splines
    alglib::spline1dinterpolant xspline, yspline;
    alglib::ae_int_t bct = 1; // first derivative (heading)
    alglib::spline1dbuildcubic(s, x, n_points, bct, slen_start * hix, 
                                               bct, slen_end * hex, xspline);
    alglib::spline1dbuildcubic(s, y, n_points, bct, slen_start * hiy, 
                                               bct, slen_end * hey, yspline);
    pair<alglib::spline1dinterpolant, alglib::spline1dinterpolant> result;
    result = {xspline, yspline};
    return result;
}







} // namespace pp_utils
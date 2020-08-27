

#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <utility>
#include <nlohmann/json.hpp>
#include <matplotlibcpp.h>

#include "interpolation.h"
#include <ap.h>
#include "ap.cpp"
#include "alglibinternal.cpp"
#include "alglibmisc.cpp"
#include "integration.cpp"
#include "interpolation.cpp"
#include "linalg.cpp"
#include "optimization.cpp"
#include "solvers.cpp"
#include "specialfunctions.cpp"

#include <Trackgraph.h>
#include <Obstacle.h>

using namespace std;
using json = nlohmann::json;
namespace plt = matplotlibcpp;


Trackgraph parseGlobalGraph(string nodeFile, string edgeFile) {
    // Load json graph from file
    std::ifstream ifs_nodes(nodeFile);
    std::ifstream ifs_edges(edgeFile);
    json nodes_json = json::parse(ifs_nodes);
    json edges_json = json::parse(ifs_edges);
    
    // Get vector list of stations - nodes
    std::vector<int> global_stations;
    for (json::iterator it = nodes_json.begin(); it != nodes_json.end(); ++it) {
        global_stations.push_back(stoi(it.key()));
    }
    sort(global_stations.begin(), global_stations.end());

    // Get vector list of stations - edges
    std::vector<int> global_stations_e;
    for (json::iterator it = edges_json.begin(); it != edges_json.end(); ++it) {
        global_stations_e.push_back(stoi(it.key()));
    }
    sort(global_stations_e.begin(), global_stations_e.end());

    // Build graph
    Trackgraph graph;
    for (auto &s : global_stations) {
        nlohmann::json layer_start = edges_json.at(to_string(s));
        for (auto &[ls_str, val] : layer_start.items()) {
            nlohmann::json pxy = nodes_json.at(to_string(s)).at(ls_str);
            double ls_flt = stod(ls_str);
            Node n(pxy.at("psi"), pxy.at("x"), pxy.at("y"));
            graph.setNode(s, ls_flt, n);
            nlohmann::json layer_end = layer_start.at(ls_str);
            for (auto &[le_str, eg] : layer_end.items()) {
                auto le_flt = stod(le_str);
                Edge edge(eg.at("x_coef"), eg.at("y_coef"), eg.at("cost"));
                graph.setEdge(s, ls_flt, le_flt, edge);
            }
        }
    }
    return graph;
}

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

vector<pair<int, double>> min_cost_path_search(Trackgraph &graph) {
    // pair<vector<double>, vector<double>> mcp;
    map<int, map<double, pair<double, double>>> costs; // mimics edge structure; pair(cpst, parent offset)
    for (auto &s : graph.nodes){
        costs[s.first];
    }
    // Define iterators for graph
    auto rit_c = graph.edges.rbegin(); // child node
    auto rit_p = rit_c; ++rit_p; // parent node

    // Initialize last layer costs to use lateral offset
    for (auto &ls : rit_c->second) {
        costs[rit_c->first][ls.first] = pair<double, double> {abs(ls.first), ls.first};
    }

    for (; rit_p != graph.edges.rend(); ++rit_p, ++rit_c) {
        int s_child = rit_c->first;
        int s_parent = rit_p->first;
        for (auto &edge : rit_p->second) { // parent node
            double l_start = edge.first;
            for (auto &dest : edge.second) { // child node
                double l_end = dest.first;
                double ec = dest.second.cost; // edge from parent to child node

                // Update parent node cost
                // auto &s_cost = costs.at(rit->first);
                if (costs.at(s_parent).find(l_start) != costs.at(s_parent).end()) { // if parent node in costs
                    // Update parent cost if (child cost + ec) is lower
                    if (ec + costs.at(s_child).at(l_end).first < costs.at(s_parent).at(l_start).first) { 
                        costs.at(s_parent).at(l_start) = 
                            pair<double, double> {(ec + costs.at(s_child).at(l_end).first), l_end};
                    }
                } else { // initialize cost at parent node if it doesn't exist
                    costs.at(s_parent)[l_start] = 
                        pair<double, double> {(ec + costs.at(s_child).at(l_end).first), l_end};
                }
            }
        }
    }
    vector<pair<int, double>> mcp;  // (station(s), offset(l))
    double cur_offset, child_offset;
    for (auto &[s, node] : costs) {
        if (s == 0) {
            cur_offset = costs.at(s).begin()->first;
            child_offset = costs.at(s).begin()->second.second;
        } else {
            child_offset = costs.at(s).at(cur_offset).second;
        }
    mcp.push_back(pair<int, double> {s, cur_offset});
    cur_offset = child_offset;
    }
    return mcp;
}

pair<vector<double>, vector<double>> mcp_to_xy(Trackgraph &graph, vector<pair<int, double>> &mcp) {
    vector<double> mcp_x, mcp_y;
        for (auto &el : mcp) {
            double x = graph.nodes.at(el.first).at(el.second).x;
            double y = graph.nodes.at(el.first).at(el.second).y;
            mcp_x.push_back(x);
            mcp_y.push_back(y);
        }
    return pair<vector<double>, vector<double>> {mcp_x, mcp_y};
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


int main() {

    // Track lattice input files    
    string nodeFile = "./data/nodes.json";
    string edgeFile = "./data/edges.json";

    // Test data
    double s_pos = 3700.432;
    double l_pos = 1.79;
    double planning_horizon = 400.0;
    Obstacles obstacles;
    obstacles.sv.push_back(StaticObstacle{3.75, 1.46, 1.8});
    obstacles.sv.push_back(StaticObstacle{-0.37, 115.6, 2.0});
    obstacles.sv.push_back(StaticObstacle{-6.18, 147.2, 2.4});
    obstacles.sv.push_back(StaticObstacle{2.98, -49.9, 3.1});

    // Parse global graph json data
    Trackgraph global_graph = parseGlobalGraph(nodeFile, edgeFile);

    auto start = std::chrono::high_resolution_clock::now();

    cout << global_graph.getNode(1105, 1.0).x << endl;
    cout << global_graph.getEdge(1105, 2.0, 3.0).cost << endl;

    // Copy local segment of graph to work off of
    Trackgraph local_graph = global_graph.extractLocalGraph(s_pos, l_pos, planning_horizon);

    // collision checking
    collision_checker(local_graph,   obstacles);
    // plotLocalGraph(local_graph, obstacles);

    // min cost graph search
    vector<pair<int, double>> mcp = min_cost_path_search(local_graph);

    // c2 trajectory gen
    pair<vector<double>, vector<double>> mcp_xy = mcp_to_xy(local_graph, mcp);
    int n_interps = 50;
    int n_points = mcp_xy.first.size();
    vector<double> svec = LinearSpacedArray(0, 1, n_points);
    int sfirst = mcp.front().first;
    int slast = mcp.back().first;
    double lfirst = mcp.front().second;
    double llast = mcp.back().second;
    double heading_start = local_graph.nodes.at(sfirst).at(lfirst).psi;
    double heading_last = local_graph.nodes.at(slast).at(llast).psi;
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


    // Get coefficients
    alglib::real_2d_array tbl;
    alglib::ae_int_t nk;
    alglib::spline1dunpack(xspline, nk, tbl);
    vector<vector<double>> coeffs(tbl.rows(), vector<double>(tbl.cols()));
    for (int row = 0; row < tbl.rows(); ++row) {
        for (int col = 0; col < tbl.cols(); ++col) {
            coeffs[row][col] = (*tbl[row, col]);
        }
    }

    double c0 = coeffs[0][2];
    double c1 = coeffs[0][3];
    double c2 = coeffs[0][3];
    double c3 = coeffs[0][4];
    vector<double> pltx, plty;
    for (auto &s : svec) {
        double valx = alglib::spline1dcalc(xspline, s);
        double valy = alglib::spline1dcalc(yspline, s);
        pltx.push_back(valx);
        plty.push_back(valy);
    }


    // Check test
    double t = 0.25;
    double v;
    v = spline1dcalc(xspline, t);
    cout << "v: " << v << endl;
    cout << *tbl[0,0] << endl;
    cout << *tbl[0,1] << endl;
    cout << *tbl[0,2] << endl;
    cout << *tbl[0,3] << endl;
    cout << *tbl[0,4] << endl;
    cout << *tbl[0,5] << endl;

    cout << "________________________" << endl;
    cout << coeffs[3][3] << endl;
    

    /////////////////////////////////////////////////////////////

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << duration.count() << " microsec" << endl;

    // plotTrack(global_graph);
    plotLocalGraph(local_graph, obstacles, mcp, pltx, plty);


    return 0;
};





 
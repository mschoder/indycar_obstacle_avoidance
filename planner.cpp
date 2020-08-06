

#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <nlohmann/json.hpp>
#include <matplotlibcpp.h>
#include <gnuplot_i.hpp>
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
    vector<double> u = LinearSpacedArray(0, 1, 100); // for spline eval
    for (const auto &s : graph.edges) {
        for (const auto &node : s.second) {
            for (const auto &dn : node.second) {
                vector<double> xspl = splineEval(u, dn.second.x_coef);
                vector<double> yspl = splineEval(u, dn.second.y_coef);
                // check each spline pt against each obs pt
                for (int i = 0; i != xspl.size(); ++i) {
                    for (auto &obs : obstacles.sv) {
                        if (collision_check_circle(xspl[i], yspl[i], obs)) {
                            graph.removeEdge(s.first, node.first, dn.first);
                        }
                    }
                }
            }
        }
    }
    // Check for nodes with no children
    for (const auto &s : graph.edges) {
        for (const auto &edge : s.second) {
            if (edge.second.empty()) {
                graph.removeNode(s.first, edge.first);
            } else if (s.first != 0) { // First node will have no parents - OK
                bool orphan = true;
                for (const auto &parent : graph.edges.at(s.first - 1)) {
                    if (parent.second.find(edge.first) != parent.second.end()) { // if node is a dest of parent
                        orphan = false;
                        break;
                    }
                if (orphan) {graph.removeNode(s.first, edge.first);}
                }
            }
        }
    }
    
    // Check for nodes with no children
}

int main() {

    // Track lattice input files    
    string nodeFile = "./data/nodes.json";
    string edgeFile = "./data/edges.json";

    // Test data
    double s_pos = 3793.432;
    double l_pos = -4.935;
    double planning_horizon = 340.0;
    Obstacles obstacles;
    obstacles.sv.push_back(StaticObstacle{1.0, 1.0, 1.2});
    obstacles.sv.push_back(StaticObstacle{-4.48, 147.2, 2.1});
    obstacles.sv.push_back(StaticObstacle{2.98, -49.9, 1.6});

    // Parse global graph json data
    Trackgraph global_graph = parseGlobalGraph(nodeFile, edgeFile);

    auto start = std::chrono::high_resolution_clock::now();

    cout << global_graph.getNode(1105, 1.0).x << endl;
    cout << global_graph.getEdge(1105, 2.0, 3.0).cost << endl;

    // Copy local segment of graph to work off of
    Trackgraph local_graph = global_graph.extractLocalGraph(s_pos, l_pos, planning_horizon);

    // collision checking
    collision_checker(local_graph, obstacles);

    // min cost graph search

    // c2 trajectory gen

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << duration.count() << " microsec" << endl;

    // plotTrack(global_graph);
    plotLocalGraph(local_graph, obstacles);


    return 0;
};





 
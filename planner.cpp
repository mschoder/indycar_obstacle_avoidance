// Extract local graph from json track lattice

#include <iostream>
#include <chrono>
#include <fstream>
#include <math.h>
#include <vector>
#include <string>
#include <iterator>
#include <nlohmann/json.hpp>
#include <matplotlibcpp.h>

using namespace std;
using json = nlohmann::json;
namespace plt = matplotlibcpp;

class trackGraph {
    public:
        struct Node {
            double psi;
            double x;
            double y;
        };
        struct Edge {
            std::array<double, 4> x_coef;
            std::array<double, 4> y_coef;
            double cost;
        };

        std::map<int, std::map<double, Node>> nodes;
        std::map<int, std::map<double, std::map<double, Edge>>> edges;
    public:
        void setNode(int s, double l, Node node) {
            nodes[s][l].x = node.x;
            nodes[s][l].y = node.y;
            nodes[s][l].psi = node.psi;
        }
        void setEdge(int s, double l_start, double l_end, Edge edge) {
            edges[s][l_start][l_end].x_coef = edge.x_coef;
            edges[s][l_start][l_end].y_coef = edge.y_coef;
            edges[s][l_start][l_end].cost = edge.cost;
        }
        Node getNode(int s, double l) {
            if (nodes.find(s) != nodes.end() & nodes[s].find(l) != nodes[s].end()){
                return nodes[s][l];
            }
        }
        Edge getEdge(int s, double l_start, double l_end) {
            if (edges.find(s) != edges.end() & edges[s].find(l_start) != edges[s].end()
                    & edges[s][l_start].find(l_end) != edges[s][l_start].end()) {
                return edges[s][l_start][l_end];
            }
        }
        std::vector<int> getNodeStations() {
            std::vector<int> stations;
            for (auto &el : nodes) {
                stations.push_back(el.first);
            }
            sort(stations.begin(), stations.end());
            return stations;
        }
        std::vector<int> getEdgeStations() {
            std::vector<int> stations;
            for (auto &el : edges) {
                stations.push_back(el.first);
            }
            sort(stations.begin(), stations.end());
            return stations;
        }
        std::vector<double> getNodeOffsets(int s) {
            if (nodes.find(s) == nodes.end()) {
                throw "station s does not exist in trackGraph";
            }
            std::vector<double> offsets;
            for (auto &el : nodes[s]) {
                offsets.push_back(el.first);
            }
            sort(offsets.begin(), offsets.end());
            return offsets;
        }
        int FindNearestStation(double s_pos) {
            std::vector<int> vec = getNodeStations();
            if (s_pos > vec.back()) {
                return FindNearestStation(s_pos - vec.back());
            }
            auto const ub = std::lower_bound(vec.begin(), vec.end(), s_pos);
            auto lb = ub;
            lb--;
            if(s_pos - *lb <= *ub - s_pos) {return *lb;} else {return *ub;}
        }
        double FindNearestOffset(int s, double value) {
            std::vector<double> vec = getNodeOffsets(s);
            if (value > vec.back()) {
                return vec.back();
            } else if (value < vec.front()) {
                return vec.front();
            }
            auto const ub = std::lower_bound(vec.begin(), vec.end(), value);
            auto lb = ub;
            lb--;
            if(value - *lb <= *ub - value) {return *lb;} else {return *ub;}
        }
        trackGraph extractLocalGraph(double s_pos, double l_pos, double planning_horizon) {
            // Find nearest station & offset
            int s_start = this->FindNearestStation(s_pos); 
            int s_end   = this->FindNearestStation(s_pos + planning_horizon);
            double l_start = this->FindNearestOffset(s_start, l_pos);


            trackGraph local_graph;
            for (auto &s : this->nodes) {
                if (s.first == s_start) { // only copy single node in first layer
                    local_graph.nodes[s_start][l_start] = this->nodes[s_start][l_start];
                    local_graph.edges[s_start][l_start] = this->edges[s_start][l_start];
                } else if (s.first > s_start & s.first <= s_end) {
                    local_graph.nodes[s.first] = this->nodes[s.first];
                    local_graph.edges[s.first] = this->edges[s.first];
                }
            }
            return local_graph;
        } 


        void removeNode(const int s, const double l) {
            if (nodes.find(s) != nodes.end() & nodes[s].find(l) != nodes[s].end()) {
                nodes[s].erase(l); // erase node
            }
            if (edges.find(s) != edges.end() & edges[s].find(l) != edges[s].end()) {
                edges[s].erase(l); // erase outgoing edges
            }
            auto prev = edges.find(s);
            prev--;
            std::cout << "prev: " << prev->first << endl;
            if (prev != edges.end()) {
                for (auto &sel : prev->second) {
                    auto &pnode = sel.first;
                    for (auto &lel : edges[prev->first][pnode]) {
                        if (lel.first == l) {
                            edges[prev->first][pnode].erase(l);  // erase incoming edges
                        }
                    }
                }
            }
        }

        void removeEdge(const int s, const double l_start, const double l_end) {
            if (edges.find(s) != edges.end() & edges[s].find(l_start) != edges[s].end()
                    & edges[s][l_start].find(l_end) != edges[s][l_start].end()) {
                edges[s][l_start].erase(l_end); // erase outgoing edges
            }
        }
};

// template <typename T, typename TT>
// int FindNearestStation(const std::vector<T>& vec, TT value) {
//         if (value > vec.back()) {
//             return FindNearestStation(vec, value - vec.back());
//         }
//         auto const ub = std::lower_bound(vec.begin(), vec.end(), value);
//         auto lb = ub;
//         lb--;
//         if(value - *lb <= *ub - value) {return *lb;} else {return *ub;}
//     }

// template <typename T, typename TT>
// int FindNearestOffset(const std::vector<T>& vec, TT value) {
//         if (value > vec.back()) {
//             return vec.back();
//         } else if (value < vec.front()) {
//             return vec.front();
//         }
//         auto const ub = std::lower_bound(vec.begin(), vec.end(), value);
//         auto lb = ub;
//         lb--;
//         if(value - *lb <= *ub - value) {return *lb;} else {return *ub;}
//     }

trackGraph parseGlobalGraph(string nodeFile, string edgeFile) {
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

    cout << "last node s: " << global_stations.back() << endl;
    cout << "last edge s: " << global_stations_e.back() << endl;

    // Build graph
    trackGraph graph;
    for (auto &s : global_stations) {
        nlohmann::json layer_start = edges_json.at(to_string(s));
        for (auto &[ls_str, val] : layer_start.items()) {
            nlohmann::json pxy = nodes_json.at(to_string(s)).at(ls_str);
            double ls_flt = stof(ls_str);
            trackGraph::Node n;
            n.x = pxy.at("x");
            n.y = pxy.at("y");
            n.psi = pxy.at("psi");
            graph.setNode(s, ls_flt, n);
            nlohmann::json layer_end = layer_start.at(ls_str);
            for (auto &[le_str, egvals] : layer_end.items()) {
                auto le_flt = stof(le_str);
                trackGraph::Edge edge;
                edge.x_coef = egvals.at("x_coef");
                edge.y_coef = egvals.at("y_coef");
                edge.cost = egvals.at("cost");
                graph.setEdge(s, ls_flt, le_flt, edge);
            }
        }
    }
    return graph;
}


// trackGraph extractLocalGraph(double s_pos, double l_pos, double planning_horizon, trackGraph &global_graph) {

//     // Find nearest station to curent pos
//     std::vector<int> global_stations;
//     for (auto &el : global_graph.nodes) {
//         global_stations.push_back(el.first);
//     }
//     sort(global_stations.begin(), global_stations.end());
    
//     // Find nearest station in graph to current pos
//     int s_start = FindNearestStation(global_stations, s_pos);
//     int s_end = FindNearestStation(global_stations, s_pos + planning_horizon);

//     // Find nearest offset to current pos
//     std::vector<double> layer_offsets;
//     for (auto &el : global_graph.nodes[s_start]) {
//         layer_offsets.push_back(el.first);
//     }
//     sort(layer_offsets.begin(), layer_offsets.end());
//     double l_local = FindNearestOffset(layer_offsets, l_pos);

//     cout << "Cur Pos (s,l): " << s_pos << " , " << l_pos << endl;
//     cout << "Start Station: " << s_start << endl;
//     cout << "Start Offset: " << l_local << endl;
//     cout << "Start S + PH: " << s_start + planning_horizon << endl;
//     cout << "End Station: " << s_end << endl;

//     trackGraph local_graph;
//     for (auto &s : global_stations) {
//         if (s == s_start) { // only copy single node in first layer
//             local_graph.nodes[s_start][l_local] = global_graph.nodes[s_start][l_local];
//             local_graph.edges[s_start][l_local] = global_graph.edges[s_start][l_local];
//         } else if (s > s_start & s <= s_end) {
//             local_graph.nodes[s] = global_graph.nodes[s];
//             local_graph.edges[s] = global_graph.edges[s];
//         }
//     }
//     return local_graph;
// }

void plotTrack(trackGraph &global_graph) {
    std::vector<double> x_nodes, y_nodes;
    for (auto &s : global_graph.getNodeStations()) {
        for (auto &el : global_graph.getNodeOffsets(s)) {
            auto node = global_graph.getNode(s, el);
            x_nodes.push_back(node.x);
            y_nodes.push_back(node.y); 
        }
    }
    plt::figure_size(780, 1200);
    plt::plot(x_nodes, y_nodes, "o");
    plt::show();
}

std::vector<double> splineEval(std::vector<double> &u, std::array<double, 4> &coefs) {
    std::vector<double> eval;
    for (auto &el : u) {
        eval.push_back(coefs[3]*pow(el,3) + coefs[2]*pow(el,2) + coefs[1]*el + coefs[0]);
    }
    return eval;
}

std::vector<double> LinearSpacedArray(double a, double b, std::size_t N) {
    double h = (b - a) / static_cast<double>(N-1);
    std::vector<double> xs(N);
    std::vector<double>::iterator x;
    double val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) {
        *x = val;
    }
    return xs;
}

// void plotLocalGraph(trackGraph &local_graph) {
//     plt::figure_size(780, 1200);
//     std::vector<double> x_nodes, y_nodes;
//     for (auto &s : local_graph.nodes) {
//         for (auto &el : s.second) {
//             auto &node = el.second;
//             x_nodes.push_back(node.x);
//             y_nodes.push_back(node.y); 
//         }
//     }
//     std::vector<double> u = LinearSpacedArray(0, 1, 50);
//     for (auto &s : local_graph.edges) {
//         for (auto &sel : s.second) {
//             for (auto &eel : sel.second) {
//                 auto &edge = eel.second;
//                 std::vector<double> x_paths = splineEval(u, edge.x_coef);
//                 std::vector<double> y_paths = splineEval(u, edge.y_coef);
//                 plt::plot(x_paths, y_paths, "skyblue");
//             }
//         }
//     }
//     plt::plot(x_nodes, y_nodes, "o");
//     plt::show();
// }



int main() 
{
    // Tmp test values
    double s_pos = 1093.432;
    double l_pos = -4.935;
    double planning_horizon = 340.0;
    string nodeFile = "./data/nodes.json";
    string edgeFile = "./data/edges.json";

    trackGraph global_graph = parseGlobalGraph(nodeFile, edgeFile);

    auto start = std::chrono::high_resolution_clock::now();

    // trackGraph local_graph = extractLocalGraph(s_pos, l_pos, planning_horizon, global_graph);
    trackGraph local_graph = global_graph.extractLocalGraph(s_pos, l_pos, planning_horizon);

    // obstacle collision checking

    // divide graph into primitives

    // sp search

    // c2 trajectory


    // local_graph.removeNode(1185, 1);
    // cout << local_graph.nodes[1185][1].psi << endl;
    // cout << local_graph.edges[1185][1][2].cost << endl; 
    // cout << local_graph.edges[1185][-3][1].x_coef[2] << endl;

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << duration.count() << " microsec" << endl;

    plotTrack(global_graph);
    // plotLocalGraph(local_graph);
    
}


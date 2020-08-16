
#ifndef PARSE_GRAPH_H
#define PARSE_GRAPH_H


#include <vector>
#include <fstream>
#include <nlohmann/json.hpp>
#include <Trackgraph.h>

using namespace std;
using json = nlohmann::json;

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

#endif
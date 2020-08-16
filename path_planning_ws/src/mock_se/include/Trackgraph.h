/* Trackgraph class definition header
 * Structure holds graph for frenet lattice discretization of racetrack
 * 
 * mschoder@mit.edu | 7/28/2020
 * 
 */

#ifndef TRACK_GRAPH_H
#define TRACK_GRAPH_H

#include <iostream>
#include <vector>
#include <map>
#include <algorithm> 
#include <math.h>

using namespace std;

struct Node {
    double psi;
    double x;
    double y;
    Node() = default;
    Node(double _psi, double _x, double _y): 
        psi(_psi), x(_x), y(_y) {} 
};

struct Edge {
    array<double, 4> x_coef;
    array<double, 4> y_coef;
    double cost;
    Edge() = default;
    Edge(array<double,4> _xc, array<double,4> _yc, double _cost): 
        x_coef(_xc), y_coef(_yc), cost(_cost) {}
};

class Trackgraph {

    public:
        map<int, map<double, Node>> nodes;
        map<int, map<double, map<double, Edge>>> edges;
        const double TRACKLENGTH = 4023.35;

        void setNode(const int& s, const double& l, Node& node);
        void setEdge(int &s, double &l_start, double &l_end, Edge &edge);

        Node& getNode(int s, double l);
        Edge& getEdge(int s, double l_start, double l_end);

        vector<int> getNodeStations();
        vector<int> getEdgeStations();
        vector<double> getNodeOffsets(int s);
        vector<double> getEdgeOffsets(int s);
        vector<double> getEdgeDests(int s, double l_start);

        int FindNearestStation(double s_pos);
        double FindNearestOffset(int s, double value);

        void removeEdge(const int s, const double l_start, const double l_end);
        void removeNode(const int s, const double l);

        pair<int,double> xy2frenet(double x, double y);

        Trackgraph extractLocalGraph(double s_pos, double l_pos, double planning_horizon);

        vector<pair<int, double>> min_cost_path_search();
};

pair<vector<double>, vector<double>> mcp_to_xy(
    Trackgraph &graph, vector<pair<int, double>> &mcp);

#endif
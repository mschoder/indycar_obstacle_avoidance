/* Trackgraph class definition header
 * Structure holds graph for frenet lattice discretization of racetrack
 * 
 * mschoder | 7/28/2020
 * 
 * TODO:
 * header guards
 * 
 */

#include <iostream>
#include <vector>
#include <map>
#include <algorithm> 

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

    public:
        const double TRACKLENGTH = 4023.35;


        void setNode(const int& s, const double& l, Node& node) {
            nodes[s][l] = Node {node.psi, node.x, node.y};
        }

        void setEdge(int& s, double& l_start, double& l_end, Edge& edge) {
            edges[s][l_start][l_end] = Edge {edge.x_coef, edge.y_coef, edge.cost};
        }

        Node& getNode(int s, double l) {
            if (nodes.find(s) != nodes.end() && nodes[s].find(l) != nodes[s].end()) {
                return nodes[s][l];
            } else throw std::invalid_argument( "No node at given indices" );
        }

        Edge& getEdge(int s, double l_start, double l_end) {
            if (edges.find(s) != edges.end() && edges[s].find(l_start) != edges[s].end()
                    && edges[s][l_start].find(l_end) != edges[s][l_start].end()) {
                return edges[s][l_start][l_end];
            } else throw std::invalid_argument( "No edge at given indices" );
        }
        
        vector<int> getNodeStations() {
            vector<int> stations;
            for (auto &el : nodes) {
                stations.push_back(el.first);
            }
            sort(stations.begin(), stations.end());
            return stations;
        }

        vector<int> getEdgeStations() {
            vector<int> stations;
            for (auto &el : edges) {
                stations.push_back(el.first);
            }
            sort(stations.begin(), stations.end());
            return stations;
        }

        vector<double> getNodeOffsets(int s) {
            if (nodes.find(s) == nodes.end()) {
                throw std::invalid_argument("station s does not exist in Trackgraph");
            }
            std::vector<double> offsets;
            for (auto &el : nodes[s]) {
                offsets.push_back(el.first);
            }
            sort(offsets.begin(), offsets.end());
            return offsets;
        }

        vector<double> getEdgeOffsets(int s) {
            if (edges.find(s) == edges.end()) {
                throw std::invalid_argument("station s does not exist in Trackgraph");
            }
            std::vector<double> offsets;
            for (auto &el : edges[s]) {
                offsets.push_back(el.first);
            }
            sort(offsets.begin(), offsets.end());
            return offsets;
        }

        vector<double> getEdgeDests(int s, double l_start) {
            if (edges.find(s) != edges.end() && edges[s].find(l_start) != edges[s].end()) {
                std::vector<double> dests;
                for (auto &el : edges[s][l_start]) {
                    dests.push_back(el.first);
                }
                return dests;
            } else throw std::invalid_argument( "No edge at given indices" );
        }

        int FindNearestStation(double s_pos) {
            std::vector<int> vec = getNodeStations();
            if (s_pos > vec.back()) {
                return FindNearestStation(s_pos - TRACKLENGTH);
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

        void removeEdge(const int s, const double l_start, const double l_end) {
            if (edges.find(s) != edges.end() && edges[s].find(l_start) != edges[s].end()
                    && edges[s][l_start].find(l_end) != edges[s][l_start].end()) {
                edges[s][l_start].erase(l_end); // erase outgoing edge
            }
        }

        void removeNode(const int s, const double l) {
            if (nodes.find(s) != nodes.end() && nodes[s].find(l) != nodes[s].end()) {
                nodes[s].erase(l); // erase node
            }
            if (edges.find(s) != edges.end() && edges[s].find(l) != edges[s].end()) {
                edges[s].erase(l); // erase all outgoing edges
            }
            auto prev = edges.find(s);
            prev--;
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

        Trackgraph extractLocalGraph(double s_pos, double l_pos, double planning_horizon) {
            // Find nearest station & offset
            int s_start = this->FindNearestStation(s_pos); 
            int s_end   = this->FindNearestStation(s_pos + planning_horizon);
            double l_start = this->FindNearestOffset(s_start, l_pos);

            // Get iterators
            auto ss_it = nodes.find(s_start);
            auto se_it = nodes.find(s_end);
            auto ss2_it = ss_it; ++ss2_it;  // one past the start node
            auto se2_it = se_it; ++se2_it;  // one past end node
            auto sfirst_it = nodes.begin(); // first node in track
            auto send_it = nodes.end();     // one past last node in track

            // Copy relevant nodes
            Trackgraph local_graph;
            int s_local = 0;

            // only copy single node in first layer
            local_graph.nodes[s_local][l_start] = this->nodes[s_start][l_start];
            local_graph.edges[s_local][l_start] = this->edges[s_start][l_start];
            ++s_local;

            if (s_start < s_end) {  // local segment doesn't cross start point
                for (auto s = ss2_it; s != se2_it; ++s) {
                    local_graph.nodes[s_local] = this->nodes[s->first];
                    local_graph.edges[s_local] = this->edges[s->first];
                    ++s_local;
                }
            } else {                // local segment does cross start point
                for (auto s = ss2_it; s != send_it; ++s) { // current loc to track start point
                    local_graph.nodes[s_local] = this->nodes[s->first];
                    local_graph.edges[s_local] = this->edges[s->first];
                    ++s_local;
                }
                for (auto s = sfirst_it; s != se2_it; ++s) { // track start pt to horizon end
                    local_graph.nodes[s_local] = this->nodes[s->first];
                    local_graph.edges[s_local] = this->edges[s->first];
                    ++s_local;
                }
            }
            return local_graph;
        }

};


/* trackGraph class definition header
 * Structure holds graph for frenet lattice discretization of racetrack
 * 
 * mschoder | 7/28/2020
 */

#include <iostream>
#include <vector>
#include <map>
#include <bits/stdc++.h> 

using namespace std;

struct Node {
    double psi;
    double x;
    double y;
};

struct Edge {
    array<double, 4> x_coef;
    array<double, 4> y_coef;
    double cost;
};

class trackGraph {

    private:
        map<int, map<double, Node>> nodes;
        map<int, map<double, map<double, Edge>>> edges;

    public:
        void setNode(const int& s, const double& l, Node& node) {
            nodes[s][l].x   = node.x;
            nodes[s][l].y   = node.y;
            nodes[s][l].psi = node.psi;
        }

        void setEdge(int& s, double& l_start, double& l_end, Edge& edge) {
            edges[s][l_start][l_end].x_coef = edge.x_coef;
            edges[s][l_start][l_end].y_coef = edge.y_coef;
            edges[s][l_start][l_end].cost   = edge.cost;
        }

        Node getNode(int s, double l) {
            if (nodes.find(s) != nodes.end() & nodes[s].find(l) != nodes[s].end()) {
                return nodes[s][l];
            } else throw std::invalid_argument( "No node at given indices" );
        }

        Edge getEdge(int s, double l_start, double l_end) {
            if (edges.find(s) != edges.end() & edges[s].find(l_start) != edges[s].end()
                    & edges[s][l_start].find(l_end) != edges[s][l_start].end()) {
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
                throw std::invalid_argument("station s does not exist in trackGraph");
            }
            std::vector<double> offsets;
            for (auto &el : nodes[s]) {
                offsets.push_back(el.first);
            }
            sort(offsets.begin(), offsets.end());
            return offsets;
        }

        vector<double> getEdgeDests(int s, double l_start) {
            if (edges.find(s) != edges.end() & edges[s].find(l_start) != edges[s].end()) {
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

        void removeEdge(const int s, const double l_start, const double l_end) {
            if (edges.find(s) != edges.end() & edges[s].find(l_start) != edges[s].end()
                    & edges[s][l_start].find(l_end) != edges[s][l_start].end()) {
                edges[s][l_start].erase(l_end); // erase outgoing edges
            }
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

        trackGraph extractLocalGraph(double s_pos, double l_pos, double planning_horizon) {
            // Find nearest station & offset
            int s_start = this->FindNearestStation(s_pos); 
            int s_end   = this->FindNearestStation(s_pos + planning_horizon);
            double l_start = this->FindNearestOffset(s_start, l_pos);

            // Copy relevant nodes
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

};


/* Trackgraph class implementation
*  mschoder@mit.edu | Modified 8/15/2020
*/

#include <Trackgraph.h>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm> 

using namespace std;

void Trackgraph::setNode(const int& s, const double& l, Node& node) {
    nodes[s][l] = Node {node.psi, node.x, node.y};
}

void Trackgraph::setEdge(int &s, double &l_start, double &l_end, Edge &edge) {
    edges[s][l_start][l_end] = Edge {edge.x_coef, edge.y_coef, edge.cost};
}

Node& Trackgraph::getNode(int s, double l) {
    if (nodes.find(s) != nodes.end() && nodes[s].find(l) != nodes[s].end()) {
        return nodes[s][l];
    } else throw std::invalid_argument( "No node at given indices" );
}

Edge& Trackgraph::getEdge(int s, double l_start, double l_end) {
    if (edges.find(s) != edges.end() && edges[s].find(l_start) != edges[s].end()
         && edges[s][l_start].find(l_end) != edges[s][l_start].end()) {
        return edges[s][l_start][l_end];
    } else throw std::invalid_argument( "No edge at given indices" );
}

vector<int> Trackgraph::getNodeStations() {
    vector<int> stations;
    for (auto &el : nodes) {
        stations.push_back(el.first);
    }
    sort(stations.begin(), stations.end());
    return stations;
}

vector<int> Trackgraph::getEdgeStations() {
    vector<int> stations;
    for (auto &el : edges) {
        stations.push_back(el.first);
    }
    sort(stations.begin(), stations.end());
    return stations;
}

vector<double> Trackgraph::getNodeOffsets(int s) {
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

vector<double> Trackgraph::getEdgeOffsets(int s) {
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

vector<double> Trackgraph::getEdgeDests(int s, double l_start) {
    if (edges.find(s) != edges.end() && edges[s].find(l_start) != edges[s].end()) {
        std::vector<double> dests;
        for (auto &el : edges[s][l_start]) {
            dests.push_back(el.first);
        }
        return dests;
    } else throw std::invalid_argument( "No edge at given indices" );
}

int Trackgraph::FindNearestStation(double s_pos) {
    std::vector<int> vec = getNodeStations();
    if (s_pos > vec.back()) {
        return FindNearestStation(s_pos - TRACKLENGTH);
    }
    auto const ub = std::lower_bound(vec.begin(), vec.end(), s_pos);
    auto lb = ub;
    lb--;
    if(s_pos - *lb <= *ub - s_pos) {return *lb;} else {return *ub;}
}

double Trackgraph::FindNearestOffset(int s, double value) {
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

void Trackgraph::removeEdge(const int s, const double l_start, const double l_end) {
    if (edges.find(s) != edges.end() && edges[s].find(l_start) != edges[s].end()
            && edges[s][l_start].find(l_end) != edges[s][l_start].end()) {
    edges[s][l_start].erase(l_end); // erase outgoing edge
    }
}

void Trackgraph::removeNode(const int s, const double l) {
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

pair<int,double> Trackgraph::xy2frenet(double x, double y) {
    // TODO FIX: deal with crossing over start line
    pair<int,double> result;
    std::vector<int> stations = getNodeStations();
    double nearest_d = 10000.;
    double next_nearest_d = 10000.;
    int nearest_s, next_nearest_s;
    pair<double,double> nearest, next_nearest;
    for (auto &s : stations) {
        auto &node = getNode(s,0);
        double xt = node.x;
        double yt = node.y;
        double dist = sqrt(pow((x-xt),2) + pow((y-yt),2));
        if (dist < nearest_d) {
            next_nearest_d = nearest_d;
            next_nearest = {nearest.first, nearest.second};
            next_nearest_s = nearest_s;
            nearest_d = dist;
            nearest = {xt, yt};
            nearest_s = s;
        }
        else if (dist < next_nearest_d) {
            next_nearest_d = dist;
            next_nearest = {xt, yt};
            next_nearest_s = s;
        }
    }
    // order s1 and s2
    if (nearest_s > next_nearest_s) {
        swap(nearest_s, next_nearest_s);
        swap(nearest_d, next_nearest_d);
        swap(nearest, next_nearest);
    }
    // interpolate between s and s_prev            
    double fract = nearest_d / (nearest_d + next_nearest_d);
    double s_est = nearest_s + (next_nearest_s - nearest_s)*fract;
    double x_est = nearest.first + (next_nearest.first - nearest.first)*fract;
    double y_est = nearest.second + (next_nearest.second - nearest.second)*fract;
    double l_est = sqrt(pow((x_est - x),2) + pow((y_est - y), 2));
    // TODO FIX: if l is to the left, it's negative

    result = {s_est, l_est};
    // cout << "x vs x_est: " << x << " " << x_est << " y vs y_est: " << y 
    //         << " " << y_est << endl;
    return result;
}


Trackgraph Trackgraph::extractLocalGraph(double s_pos, double l_pos, double planning_horizon) {
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
    } else { // local segment does cross start point
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

vector<pair<int, double>> Trackgraph::min_cost_path_search() {
    map<int, map<double, pair<double, double>>> costs; // mimics edge structure; pair(cpst, parent offset)
    for (auto &s : nodes){
        costs[s.first];
    }
    // Define iterators for graph
    auto rit_c = edges.rbegin(); // child node
    auto rit_p = rit_c; ++rit_p; // parent node

    // Initialize last layer costs to use lateral offset
    for (auto &ls : rit_c->second) {
        costs[rit_c->first][ls.first] = pair<double, double> {abs(ls.first), ls.first};
    }

    for (; rit_p != edges.rend(); ++rit_p, ++rit_c) {
        int s_child = rit_c->first;
        int s_parent = rit_p->first;
        for (auto &edge : rit_p->second) { // parent node
            double l_start = edge.first;
            for (auto &dest : edge.second) { // child node
                double l_end = dest.first;
                double ec = dest.second.cost; // edge from parent to child node

                // Update parent node cost
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
/* Defines classes for static (stationary) and moving obtacles
*  for use with the path planner
*
*  mschoder | 8/6/20
*/

#include <iostream>
#include <vector>
#include <algorithm> 
#include <utility>
#include <math.h>

using namespace std;


struct StaticObstacle {
    double x;      // center coordintates (m)
    double y;
    double rad;    // radius (m)
};

struct DynamicObstacle {
    double x; 
    double y;
    double x_dot;
    double y_dot;
    double rad;
};

struct Obstacles {
    vector<StaticObstacle> sv;
    vector<DynamicObstacle> dv;
};

pair<vector<double>, vector<double>> circle_points(StaticObstacle &obs) {
    pair<vector<double>, vector<double>> result;
    for (double angle=0; angle<=2*M_PI; angle+=0.1) {
        result.first.push_back(obs.x + obs.rad*cos(angle));
        result.second.push_back(obs.y + obs.rad*sin(angle));
    }
    return result;
}

bool collision_check_circle(const Node node, const StaticObstacle &obs) {
    // Return true if collision; false otherwise
    return obs.rad >= sqrt(pow(node.x - obs.x, 2) + pow(node.y - obs.y, 2)) ? true : false;
}

bool collision_check_circle(const double x, const double y, const StaticObstacle &obs) {
    return obs.rad >= sqrt(pow(x - obs.x, 2) + pow(y - obs.y, 2)) ? true : false;
}
    

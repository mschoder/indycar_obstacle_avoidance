/* Defines classes for static (stationary) and moving obtacles
*  for use with the path planner
*
*  mschoder | 8/6/20
*/

#ifndef OBSTACLE_H
#define OBSTACLE_H

#include <iostream>
#include <vector>
#include <algorithm> 
#include <utility>
#include <math.h>
#include <Trackgraph.h>

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
    

#endif
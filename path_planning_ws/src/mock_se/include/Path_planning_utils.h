
#ifndef PATH_PLANNING_UTILS
#define PATH_PLANNING_UTILS

#include <vector>
#include <array>
#include <Trackgraph.h>
#include <Obstacle.h>
#include <matplotlibcpp.h>
#include "interpolation.h"

using namespace std;

namespace pp_utils {

vector<double> splineEval(const vector<double> &u, const array<double, 4> &coefs);
vector<double> LinearSpacedArray(double a, double b, std::size_t N);
void collision_checker(Trackgraph &graph, Obstacles &obstacles);
pair<vector<double>, vector<double>> circle_points(StaticObstacle &obs);
bool collision_check_circle(const Node &node, const StaticObstacle &obs);
bool collision_check_circle(const double x, const double y, const StaticObstacle &obs);


void plotTrack(Trackgraph &graph);
void plotLocalGraph(Trackgraph &graph, Obstacles &obstacles);
void plotLocalGraph(Trackgraph &graph, Obstacles &obstacles, 
        vector<pair<int,double>> &mcp, vector<double> &pltx, vector<double> &plty);

// Trajectory generation
pair<alglib::spline1dinterpolant, alglib::spline1dinterpolant> trajectory_gen(
    Trackgraph &graph, vector<pair<int, double>> &mcp);

} // namespace pp_utils

#endif

#ifndef CUBIC_SPLINE_H
#define CUBIC_SPLINE_H

#include<vector>

using namespace std;
using vec = vector<double>;

struct SplineSet{
    double a;
    double b;
    double c;
    double d;
    double x;
};

vector<SplineSet> cubic_spline(vec &x, vec &y);

#endif
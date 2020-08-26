/* 
cubic spline generation; curvature-continuous
https://en.wikipedia.org/w/index.php?title=Spline_%28mathematics%29&oldid=288288033#Algorithm_for_computing_natural_cubic_splines

spline form: y_j(x_j) = a_j + b_j(x - x_j) + c_j(x - x_j)^2 + d_j(x - x_j)^3
 ** j splines, where j = n - 1 and n = number of knot points
 ** C2, second derivatives of endpoints are set to zero curvature
*/

#include<vector>
#include<algorithm>
#include<cmath>
#include<cubic_spline.h>

using namespace std;
using vec = vector<double>;

vector<SplineSet> cubic_spline(vec &x, vec &y)
{
    int n = x.size()-1;
    vec a;
    a.insert(a.begin(), y.begin(), y.end());
    vec b(n);
    vec d(n);
    vec h;

    for(int i = 0; i < n; ++i)
        h.push_back(x[i+1]-x[i]);

    vec alpha;
    alpha.push_back(0);
    for(int i = 1; i < n; ++i)
        alpha.push_back( 3*(a[i+1]-a[i])/h[i] - 3*(a[i]-a[i-1])/h[i-1]  );

    vec c(n+1);
    vec l(n+1);
    vec mu(n+1);
    vec z(n+1);
    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;

    for(int i = 1; i < n; ++i)
    {
        l[i] = 2 *(x[i+1]-x[i-1])-h[i-1]*mu[i-1];
        mu[i] = h[i]/l[i];
        z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
    }

    l[n] = 1;
    z[n] = 0;
    c[n] = 0;

    for(int j = n-1; j >= 0; --j)
    {
        c[j] = z [j] - mu[j] * c[j+1];
        b[j] = (a[j+1]-a[j])/h[j]-h[j]*(c[j+1]+2*c[j])/3;
        d[j] = (c[j+1]-c[j])/3/h[j];
    }

    vector<SplineSet> output_set(n);
    for(int i = 0; i < n; ++i)
    {
        output_set[i].a = a[i];
        output_set[i].b = b[i];
        output_set[i].c = c[i];
        output_set[i].d = d[i];
        output_set[i].x = x[i];
    }
    return output_set;
}
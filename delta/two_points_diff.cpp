// Practical course on computing derivatives in code: https://dl.acm.org/doi/10.1145/3305366.3328073

#include <iostream>
#include <cmath>
using namespace std;
// f(x+h)= f(x) + hf'(x) + O(h)
// O(h) ~ (h/2*f''(x))
double f(double x)
{
    return 1.0 / x;
}

void compute(double x, double &dx, double &dxx)
{
    dx = -1.0 / (x * x);
    dxx = 2.0 / (x * x * x);
}

int main()
{
    // f(x) = 1/x
    // exact
    double x, dx, dxx;
    x = 2;
    compute(x, dx, dxx);
    // differences to approximate
    double h = 0.1;
    double dx_ = (f(x + h) - f(x)) / h;
    printf("%f %f\n", fabs(dx_ - dx), fabs(dxx * 0.5 * h));

    return 0;
}
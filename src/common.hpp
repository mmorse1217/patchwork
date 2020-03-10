#ifndef __COMMON_HPP__
#define __COMMON_HPP__
// Junk that doesn't have a home yet.
typedef pair<double, double> Interval;
typedef pair<Interval, Interval> Rectangle;
typedef tuple<Interval, Interval, Interval> Cube;
typedef pair<int, int> Index;
namespace Common {
    enum BasisType {
        BEZIER = 0,
        CHEBYSHEV =1,
        BSPLINE =2
    };

    enum Direction{
        U = 0,
        V = 1
    };
    enum Periodicity{
        PERIODIC_XY = 0
    };
};
#endif

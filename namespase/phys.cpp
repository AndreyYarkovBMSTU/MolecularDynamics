#include "phys.h"

double phys::getClausiusMossotti(double eps_p, double eps_s)
{
    return (eps_p - eps_s) / (eps_p + 2 * eps_s);
}

// Length
double phys::m()
{
    return 1.0;
}
double phys::dm()
{
    return 1e-1;
}
double phys::sm()
{
    return 1e-2;
}
double phys::mm()
{
    return 1e-3;
}
double phys::um()
{
    return 1e-6;
}
double phys::nm()
{
    return 1e-9;
}
double phys::pm()
{
    return 1e-12;
}
double phys::fm()
{
    return 1e-15;
}
double phys::am()
{
    return 1e-18;
}

// Mass
double phys::kg()
{
    return 1.0;
}
double phys::g()
{
    return 1e-2;
}
double phys::mg()
{
    return 1e-3;
}
double phys::ug()
{
    return 1e-6;
}
double phys::ng()
{
    return 1e-9;
}
double phys::pg()
{
    return 1e-12;
}
double phys::fg()
{
    return 1e-15;
}
double phys::ag()
{
    return 1e-18;
}

// Time
double phys::s()
{
    return 1.0;
}
double phys::min()
{
    return 60.0;
}
double phys::h()
{
    return min() * 60.0;
}
double phys::ms()
{
    return 1e-3;
}
double phys::us()
{
    return 1e-6;
}
double phys::ns()
{
    return 1e-9;
}
double phys::ps()
{
    return 1e-12;
}
double phys::fs()
{
    return 1e-15;
}
double phys::as()
{
    return 1e-18;
}

// Temperature
double phys::K()
{
    return 1.0;
}

// Force
double phys::N()
{
    return kg() * m() / pow(s(), 2);
}

// Energy, work, heat
double phys::J()
{
    return N() * m();
}

// Pressure
double phys::Pa()
{
    return N() / pow(m(), 2);
}

// Gauss distribution
double phys::getGaussian(double _t, double _tm, double _sigma)
{
    return (exp(-pow(_t - _tm, 2) / (2 * pow(_sigma, 2))) / (_sigma * pow(2 * M_PI, 0.5)));
}

#ifndef PHYS_H
#define PHYS_H

#include <math.h>

#pragma once

namespace phys
{
    double getClausiusMossotti(double eps_p, double eps_s)
    {
        return (eps_p - eps_s) / (eps_p + 2 * eps_s);
    }

    // Length
    double m()
    {
        return 1.0;
    }
    double dm()
    {
        return 1e-1;
    }
    double sm()
    {
        return 1e-2;
    }
    double mm()
    {
        return 1e-3;
    }
    double um()
    {
        return 1e-6;
    }
    double nm()
    {
        return 1e-9;
    }
    double pm()
    {
        return 1e-12;
    }
    double fm()
    {
        return 1e-15;
    }
    double am()
    {
        return 1e-18;
    }

    // Mass
    double kg()
    {
        return 1.0;
    }
    double g()
    {
        return 1e-2;
    }
    double mg()
    {
        return 1e-3;
    }
    double ug()
    {
        return 1e-6;
    }
    double ng()
    {
        return 1e-9;
    }
    double pg()
    {
        return 1e-12;
    }
    double fg()
    {
        return 1e-15;
    }
    double ag()
    {
        return 1e-18;
    }

    // Time
    double s()
    {
        return 1.0;
    }
    double min()
    {
        return 60.0;
    }
    double h()
    {
        return min() * 60.0;
    }
    double ms()
    {
        return 1e-3;
    }
    double us()
    {
        return 1e-6;
    }
    double ns()
    {
        return 1e-9;
    }
    double ps()
    {
        return 1e-12;
    }
    double fs()
    {
        return 1e-15;
    }
    double as()
    {
        return 1e-18;
    }

    // Temperature
    double K()
    {
        return 1.0;
    }

    // Force
    double N()
    {
        return kg() * m() / pow(s(), 2);
    }

    // Energy, work, heat
    double J()
    {
        return N() * m();
    }

    // Pressure
    double Pa()
    {
        return N() / pow(m(), 2);
    }

    // Gauss distribution
    double getGaussian(double _t, double _tm, double _sigma)
    {
        return (exp(-pow(_t - _tm, 2) / (2 * pow(_sigma, 2))) / (_sigma * pow(2 * M_PI, 0.5)));
    }
}

#endif // PHYS_H

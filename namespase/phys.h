#ifndef PHYS_H
#define PHYS_H

#include <math.h>

#pragma once

namespace phys
{
    double getClausiusMossotti(double eps_p, double eps_s);

    // Length
    double m();
    double dm();
    double sm();
    double mm();
    double um();
    double nm();
    double pm();
    double fm();
    double am();

    // Mass
    double kg();
    double g();
    double mg();
    double ug();
    double ng();
    double pg();
    double fg();
    double ag();

    // Time
    double s();
    double min();
    double h();
    double ms();
    double us();
    double ns();
    double ps();
    double fs();
    double as();

    // Temperature
    double K();

    // Force
    double N();

    // Energy, work, heat
    double J();

    // Pressure
    double Pa();

    // Gauss distribution
    double getGaussian(double _t, double _tm, double _sigma);
    double getDistribFunc(double _t, double _tm, double _sigma);
}

#endif // PHYS_H

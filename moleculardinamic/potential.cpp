#include "potential.h"

double IPL::getPotential(Vector _r_1, Vector _r_2)
{
    return pow(2.0 / ((_r_1 - _r_2).norm() * scale), n);
}

Vector IPL::getGradPotential(Vector _r_1, Vector _r_2)
{
    return - n * getPotential(_r_1, _r_2) * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 2) * (_r_1 - _r_2) * scale;
}

double LJ::getPotential(Vector _r_1, Vector _r_2)
{
        return pow(2.0 / ((_r_1 - _r_2).norm() * scale), 12) - pow(2.0 / ((_r_1 - _r_2).norm() * scale), 6);
}


Vector LJ::getGradPotential(Vector _r_1, Vector _r_2)
{
        return - 12 * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 12) * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 2) * (_r_1 - _r_2) * scale + 6 * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 6) * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 2) * (_r_1 - _r_2) * scale;
}

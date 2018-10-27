#include "potential.h"

double IPL::getPotential(Vector _r_1, Vector _r_2)
{
    return pow(2.0 / ((_r_1 - _r_2).norm() * scale), n);
}

Vector IPL::getGradPotential(Vector _r_1, Vector _r_2)
{
    return - n * getPotential(_r_1, _r_2) * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 2) * (_r_1 - _r_2) * scale;
}


double IPL::getPotential_dipole(Vector _r_1, Vector _r_2, Vector _dipolemoment_1, Vector _dipolemoment_2)
{
    return pow(2.0 / ((_r_1 - _r_2).norm() * scale), n);
}

Vector IPL::getGradPotential_dipole(Vector _r_1, Vector _r_2, Vector _dipolemoment_1, Vector _dipolemoment_2)
{
    return - n * getPotential(_r_1, _r_2) * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 2) * (_r_1 - _r_2) * scale;
}

double LJ::getPotential(Vector _r_1, Vector _r_2)
{
        return 4 * (pow(2.0 / ((_r_1 - _r_2).norm() * scale), 12)
                    - pow(2.0 / ((_r_1 - _r_2).norm() * scale), 6));
}

Vector LJ::getGradPotential(Vector _r_1, Vector _r_2)
{
        return 4 * (- 12 * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 12) * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 2) * (_r_1 - _r_2) * scale);
                    //+ 6 * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 6) * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 2) * (_r_1 - _r_2) * scale);
}

double LJ::getPotential_dipole(Vector _r_1, Vector _r_2, Vector _dipolemoment_1, Vector _dipolemoment_2)
{
        return 4 * (pow(2.0 / ((_r_1 - _r_2).norm() * scale), 12)
                    - pow(2.0 / ((_r_1 - _r_2).norm() * scale), 6));
}


Vector LJ::getGradPotential_dipole(Vector _r_1, Vector _r_2, Vector _dipolemoment_1, Vector _dipolemoment_2)
{
        return 4 * (- 12 * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 12) * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 2) * (_r_1 - _r_2) * scale
                    + 6 * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 6) * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 2) * (_r_1 - _r_2) * scale);
}

double StockMayer::getPotential(Vector _r_1, Vector _r_2)
{
        return 4 * (pow(2.0 / ((_r_1 - _r_2).norm() * scale), 12)
                    - pow(2.0 / ((_r_1 - _r_2).norm() * scale), 6));
}


Vector StockMayer::getGradPotential(Vector _r_1, Vector _r_2)
{
        return 4 * (- 12 * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 12) * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 2) * (_r_1 - _r_2) * scale
                    + 6 * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 6) * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 2) * (_r_1 - _r_2) * scale);
}

double StockMayer::getPotential_dipole(Vector _r_1, Vector _r_2, Vector _dipolemoment_1, Vector _dipolemoment_2)
{
        return 4 * (pow(2.0 / ((_r_1 - _r_2).norm() * scale), 12)
                    - pow(2.0 / ((_r_1 - _r_2).norm() * scale), 6))
                - _dipolemoment_1.dot(_dipolemoment_2) / pow((_r_1 - _r_2).norm() * scale, 3);
}


Vector StockMayer::getGradPotential_dipole(Vector _r_1, Vector _r_2, Vector _dipolemoment_1, Vector _dipolemoment_2)
{
        return 4 * (- 12 * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 12) * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 2) * (_r_1 - _r_2) * scale
                    + 6 * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 6) * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 2) * (_r_1 - _r_2) * scale)
                + 3 * (_dipolemoment_1.dot(_dipolemoment_2) / pow((_r_1 - _r_2).norm() * scale, 3)) * (_r_1 - _r_2) * scale;
}

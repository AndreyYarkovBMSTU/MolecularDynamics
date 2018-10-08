#include "thermostat.h"

Vector Langevin::getForce()
{
    return Vector(phys::getGaussian((-300 + std::rand() % 600) / 100, 0.0, 1.0), phys::getGaussian((-300 + std::rand() % 600) / 100, 0.0, 1.0), 0.0);
//        return Vector((-2000 + std::rand() % 4000), (-2000 + std::rand() % 4000), 0.0) * 1e3 * prop->k * prop->temperature / prop->radius - Vector((frictionkoef *  _velocity(0)), (frictionkoef * _velocity(1)), 0.0);
}


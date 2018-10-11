#include "thermostat.h"

Vector Langevin::getForce()
{
    std::mt19937 generator(rd());
    std::normal_distribution<> distribution(0, 1);

    return Vector(distribution(generator), distribution(generator), 0.0) * prop->koef_Brown;

//        return Vector((-2000 + std::rand() % 4000), (-2000 + std::rand() % 4000), 0.0) * 1e3 * prop->k * prop->temperature / system->particles[0]->radius - Vector((frictionkoef *  _velocity(0)), (frictionkoef * _velocity(1)), 0.0);
}


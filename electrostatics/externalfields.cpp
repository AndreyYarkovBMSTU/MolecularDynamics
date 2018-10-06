#include "externalfields.h"

Vector DirectedField::getElectricField(int _numAngle = 0)
{
    return electricfield;
}

Vector RotatingField::getElectricField(int _numAngle)
{
    return Vector(electricfield.norm() * cos(_numAngle * omega * prop->timestep), electricfield.norm() * sin(_numAngle *omega * prop->timestep), 0.0);
}

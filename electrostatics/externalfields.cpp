#include "externalfields.h"

Vector DirectedField::getElectricField(int _numAngle = 0, double _timestep = 0.0)
{
    return electricfield;
}

Vector RotatingField::getElectricField(int _numAngle, double _timestep)
{
    return Vector(electricfield.norm() * cos(_numAngle * omega * _timestep), electricfield.norm() * sin(_numAngle * omega * _timestep), 0.0);
}

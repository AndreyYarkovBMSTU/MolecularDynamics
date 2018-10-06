#ifndef INTERACTION_H
#define INTERACTION_H

#include "header.h"
#include "particle/particlesystem.h"
#include "moleculardinamic/potential.h"
#include "methods.h"

class Interaction
{
public:
    ParticleSystem* system;
    Methods* method;
    Properties* prop;

    Interaction(ParticleSystem* _system, Methods* _method, Properties* _prop) :
        system(_system), method(_method), prop(_prop)
    {

    }

    Vector getElectricForce(int _nParticle);

    void recordPotentials(int _numPoints, std::string _energytype);

    void recordEnergy(int _nParticle, std::string _energytype);

    double getAverageEnergy(int _nParticle, std::string _energytype);
};

#endif // INTERACTION_H

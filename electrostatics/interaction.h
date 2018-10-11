#ifndef INTERACTION_H
#define INTERACTION_H

#include "header.h"
#include "particle/particlesystem.h"
#include "moleculardinamic/potential.h"
#include "methods.h"

class Interaction
{
private:
    int n;
    double average_energy;
    double dr;
    double U;
    double U_;
    double phi;
    double Ex_;
    double Ey_;
    Vector r_;
    Vector force;
    Vector electricfield_;
    std::string path;
    std::string path0;
    std::ofstream out;
    std::ofstream out0;
public:
    ParticleSystem* system;
    Methods* method;
    Properties* prop;

    Interaction(ParticleSystem* _system, Methods* _method, Properties* _prop) :
        system(_system), method(_method), prop(_prop)
    {
        dr = 1e-1;
    }

    double getGradEnergy(int _nParticle, int _nCoord);

    Vector getElectricForce(int _nParticle);

    void recordPotentials(int _numPoints, std::string _energytype);

    void recordEnergy(int _nParticle, std::string _energytype);

    double getAverageEnergy(int _nParticle, std::string _energytype);
};

#endif // INTERACTION_H

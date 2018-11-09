#ifndef INTERACTION_H
#define INTERACTION_H

#include "header.h"
#include "particle/particlesystem.h"
#include "moleculardinamic/potential.h"
#include "methods.h"

class Interaction
{
private:
    int numpoints;
    int numParticles;
    int m;
    int count;
    double average_energy;
    double dr;
    double U;
    double U_;
    double phi;
    double amplitude;
    double Ex_;
    double Ey_;
    Vector r_;
    Vector force;
    Vector electricfield_;
    Matrix Potentials;
    Matrix ElectricField;
    std::string path;
    std::string path0;
    std::ofstream out;
    std::ofstream out0;
    std::ifstream in;
    std::vector<Vector> points;
    std::vector<State*> states;
    std::vector<Particle*> cluster;
    ExternalFields* externalfield0;
    Environment* environment0;
    Methods* method0;
    ParticleSystem* system0;
public:
    ParticleSystem* system;
    Methods* method;
    Properties* prop;
    std::string accuracyEnergy;

    Interaction(ParticleSystem* _system, Methods* _method, Properties* _prop, std::string _accuracyEnergy) :
        system(_system), method(_method), prop(_prop), accuracyEnergy(_accuracyEnergy)
    {
        dr = 1e-1;
    }

    double getGradEnergy(int _nParticle, int _nCoord);
    double getGradAverageEnergy(int _nParticle, int _nCoord);

    Vector getElectricForce(int _nParticle);
    Vector getElectricForceDipole(int _nParticle);

    void recordPotentials(int _numPoints, std::string _energytype);
    void recordEnergy(int _nParticle, std::string _energytype);

    double getAverageEnergy(int _nParticle, std::string _energytype);

    void recordPlane(std::string _pathinParticles, std::string _pathinPoints, std::string _pathoutPot, std::string _pathoutEl);
};

#endif // INTERACTION_H

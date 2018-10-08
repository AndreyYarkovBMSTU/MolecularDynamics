#ifndef METHODS_H
#define METHODS_H

#include "header.h"
#include "particle/particlesystem.h"

class Methods
{
private:

public:
    ParticleSystem* system;
    Properties* prop;

    Methods(ParticleSystem* _system, Properties* _prop) :
        system(_system), prop(_prop)
    {

    }

    virtual void setDipoleMoment() = 0;
};

class NonIteractingDipoles : public Methods
{
private:
    Vector dipolemoment;
public:
    NonIteractingDipoles(ParticleSystem* _system, Properties* _prop) :
        Methods(_system, _prop)
    {

    }

    void setDipoleMoment();
};

class SelfConsistentDipoles : public Methods
{
private:
    int a, b, k;
    double kronec_ab;
    double kronec_ij;
    Vector n_ab;
    Matrix Block;
    Matrix Block_Matrix;
    vector<Matrix> B;
    Matrix E;
    Matrix Dipolemoment;
public:
    SelfConsistentDipoles(ParticleSystem* _system, Properties* _prop) :
        Methods(_system, _prop)
    {

    }

    Matrix getBlock(Particle* particle1, Particle* particle2);

    Matrix getInteractionsMatrix();

    void setDipoleMoment();
};

#endif // METHODS_H

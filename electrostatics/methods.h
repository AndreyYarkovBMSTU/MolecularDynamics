#ifndef METHODS_H
#define METHODS_H

#include "header.h"
#include "particle/particlesystem.h"

class Methods
{
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
public:
    NonIteractingDipoles(ParticleSystem* _system, Properties* _prop) :
        Methods(_system, _prop)
    {

    }

    void setDipoleMoment();
};

class SelfConsistentDipoles : public Methods
{
public:
    SelfConsistentDipoles(ParticleSystem* _system, Properties* _prop) :
        Methods(_system, _prop)
    {

    }

    Matrix getBlock(int _a, int _b, Vector _ra, Vector _rb);

    Matrix getInteractionsMatrix();

    void setDipoleMoment();
};

#endif // METHODS_H

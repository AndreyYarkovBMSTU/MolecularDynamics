#include "interaction.h"

double Interaction::getGradEnergy(int _nParticle, int _nCoord)
{
    system->particles[_nParticle]->state->r(_nCoord) = r_(_nCoord) + dr;
    U_ = system->getInteractionEnergy(_nParticle).energy + system->getSelfEnergy(_nParticle).energy;
    system->particles[_nParticle]->state->r(_nCoord) = r_(_nCoord) - dr;
    U_ = U_ - system->getInteractionEnergy(_nParticle).energy + system->getSelfEnergy(_nParticle).energy;
    system->particles[_nParticle]->state->r(_nCoord) = r_(_nCoord);

    return (U_ - U) / (2 * dr);
}

Vector Interaction::getElectricForce(int _nParticle)
{
    method->setDipoleMoment();
    r_ = system->particles[_nParticle]->getCoordinate();

    return Vector(- getGradEnergy(_nParticle, 0), - getGradEnergy(_nParticle, 1), 0.0);
}

Vector Interaction::getElectricForceDipole(int _nParticle)
{
    method->setDipoleMoment();
    force = Vector(0.0, 0.0, 0.0);
    for (int i = 0; i < system->numParticles; i++)
    {
        if (i != _nParticle)
        {
            force += Vector(0.0, 0.0, 0.0);
        }
    }

    return force;
}

void Interaction::recordPotentials(int _numPoints, std::string _energytype)
{
    r_ = system->particles[1]->getCoordinate();
    path0 = prop->path + "Self" + prop->filetype;
    out0.open(path0, std::ios::out | std::ios::binary);
    out0 << system->dipolemoment0.dot(system->dipolemoment0) / (4 * pow(2 * system->particles[0]->radius, 3));
    out0.close();
    for (int i = 0; i < _numPoints; i++)
    {
        method->setDipoleMoment();
        if (i < 10)
        {
            system->particles[1]->setCoordinate(Vector(pow(10, log10(2.01 + i * 0.01)), r_(1), r_(2)));
        }
        else
        {
            system->particles[1]->setCoordinate(Vector(pow(10, log10(2.1 + (i - 10) * 0.1)),  r_(1), r_(2))); // 2.5304
        }

        path = prop->path + _energytype + "Average" + std::to_string(i) + prop->filetype;
        out.open(path, std::ios::out | std::ios::binary);
        out << system->particles[1]->getCoordinate()(0) / (2 * system->particles[0]->radius) << " " << getAverageEnergy(0, _energytype);
        out.close();

        if (system->particles[1]->getCoordinate()(0) / (2 * system->particles[0]->radius) >= 10.0)
        {
            system->particles[1]->setCoordinate(Vector(1e10, r_(1), r_(2)));
            std::cout << _energytype << " " << getAverageEnergy(1, _energytype) << std::endl;
            system->particles[1]->setCoordinate(r_);
            break;
        }
        //system->particles[1]->state->r(0) += 10.0 *  2 * system->particles[0]->radius / _numPoints;
    }
}

void Interaction::recordEnergy(int _nParticle, std::string _energytype)
{
    electricfield_ = system->environment->externalfield->electricfield;

    if (system->environment->externalfield->omega != 0.0)
    {
        for (int i = 0; i <= prop->numAngles; i++)
        {
            phi = i * M_PI / 2 / prop->numAngles;
            Ex_ = electricfield_.norm() * cos(phi);
            Ey_ = electricfield_.norm() * sin(phi);
            system->environment->externalfield->electricfield = Vector(Ex_, Ey_, 0.0);
            method->setDipoleMoment();

            path = prop->path + _energytype + std::to_string(i) + prop->filetype;
            out.open(path, std::ios::out | std::ios::binary);
            if (_energytype == "induction")
            {
                out << phi << " " << system->getInductionEnergy(_nParticle).energy;
            }
            else if (_energytype == "self")
            {
                out << phi << " " << system->getSelfEnergy(_nParticle).energy;
            }
            else
            {
                out << phi << " " << 0.0;
            }
            out.close();
        }
    }
    system->environment->externalfield->electricfield = electricfield_;
}

double Interaction::getAverageEnergy(int _nParticle, std::string _energytype)
{
    electricfield_ = system->environment->externalfield->electricfield;
    if (system->environment->externalfield->omega != 0.0)
    {
        for (int i = 0; i < prop->numAngles; i++)
        {
            phi = i * 2 * M_PI / prop->numAngles;
            Ex_ = electricfield_.norm() * cos(phi);
            Ey_ = electricfield_.norm() * sin(phi);
            system->environment->externalfield->electricfield = Vector(Ex_, Ey_, 0.0);
            method->setDipoleMoment();
            if (_energytype == "induction")
            {
                average_energy = system->getInductionEnergy(_nParticle).energy;
            }
            else if (_energytype == "interaction")
            {
                average_energy = system->getInteractionEnergy(_nParticle).energy;
            }
            else if (_energytype == "self")
            {
                average_energy = system->getSelfEnergy(_nParticle).energy;
            }
        }
        n = prop->numAngles;
    }

    if (_energytype == "ipl3")
    {
        average_energy = - system->particles[1]->dipolemoment.dot(system->particles[1]->dipolemoment) * pow((1.0 / (system->particles[1]->getCoordinate() - system->particles[0]->getCoordinate()).norm()), 3) / 4.0;
        n = 1;
    }
    system->environment->externalfield->electricfield = electricfield_;

    return average_energy / n;
}

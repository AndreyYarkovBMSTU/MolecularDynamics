#include "interaction.h"

Vector Interaction::getElectricForce(int _nParticle)
{
    method->setDipoleMoment();
    Vector r_ = system->particles[_nParticle]->state->r;
    double dr = prop->radius * 1e-1;
//! ERROR!!!!
    double U = system->getInteractionEnergy(_nParticle).energy + system->getSelfEnergy(_nParticle).energy;
    system->particles[_nParticle]->state->r(0) = r_(0) + dr;
    double Ux = system->getInteractionEnergy(_nParticle).energy + system->getSelfEnergy(_nParticle).energy;
    system->particles[_nParticle]->state->r(0) = r_(0) - dr;
    Ux = Ux - system->getInteractionEnergy(_nParticle).energy + system->getSelfEnergy(_nParticle).energy;
    system->particles[_nParticle]->state->r(0) = r_(0);
    system->particles[_nParticle]->state->r(1) = r_(1) + dr;
    double Uy = system->getInteractionEnergy(_nParticle).energy + system->getSelfEnergy(_nParticle).energy;
    system->particles[_nParticle]->state->r(1) = r_(1) - dr;
    Uy = Uy - system->getInteractionEnergy(_nParticle).energy + system->getSelfEnergy(_nParticle).energy;
    system->particles[_nParticle]->state->r(1) = r_(1);

//        system->particles[_nParticle]->state->r(2) = r_(2) + dr;
//        double Uz = system->getInteractionEnergy(_nParticle).energy + system->getSelfEnergy(_nParticle).energy;
//        system->particles[_nParticle]->state->r(2) = r_(2) - dr;
//        Uz = Uz - system->getInteractionEnergy(_nParticle).energy + system->getSelfEnergy(_nParticle).energy;
//        system->particles[_nParticle]->state->r(2) = r_(2);

    return Vector(- (Ux - U) / (2 * dr), - (Uy - U) / (2 * dr), 0.0);
}

void Interaction::recordPotentials(int _numPoints, std::string _energytype)
{
    Vector r_ = system->particles[1]->state->r;
    std::string path0 = prop->path + "Self" + prop->filetype;
    std::ofstream out0;
    out0.open(path0, std::ios::out | std::ios::binary);
    out0 << system->dipolemoment0.dot(system->dipolemoment0) / (4 * pow(2 * prop->radius, 3));
    out0.close();
    for (int i = 0; i < _numPoints; i++)
    {
        method->setDipoleMoment();
        if (i < 10)
        {
            system->particles[1]->state->r(0) = pow(10, log10(2.01 + i * 0.01));
        }
        else
        {
            system->particles[1]->state->r(0) = pow(10, log10(2.1 + (i - 10) * 0.1)); // 2.5304
        }

        std::string path = prop->path + _energytype + "Average" + std::to_string(i) + prop->filetype;
        std::ofstream out;
        out.open(path, std::ios::out | std::ios::binary);
        out << system->particles[1]->state->r(0) / (2 * prop->radius) << " " << getAverageEnergy(0, _energytype);
        out.close();

        if (system->particles[1]->state->r(0) / (2 * prop->radius) >= 10.0)
        {
            system->particles[1]->state->r(0) = 1e10;
            std::cout << _energytype << " " << getAverageEnergy(1, _energytype) << std::endl;
            system->particles[1]->state->r = r_;
            break;
        }
        //system->particles[1]->state->r(0) += 10.0 *  2 * prop->radius / _numPoints;
    }
}

void Interaction::recordEnergy(int _nParticle, std::string _energytype)
{
    Vector electricfield_ = system->environment->externalfield->electricfield;

    if (system->environment->externalfield->omega != 0.0)
    {
        for (int i = 0; i <= prop->numAngles; i++)
        {
            double phi = i * M_PI / 2 / prop->numAngles;
            double Ex_ = electricfield_.norm() * cos(phi);
            double Ey_ = electricfield_.norm() * sin(phi);
            system->environment->externalfield->electricfield = Vector(Ex_, Ey_, 0.0);
            method->setDipoleMoment();
            std::string path = prop->path + _energytype + std::to_string(i) + prop->filetype;
            std::ofstream out;
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
    int n;
    double average_energy;
    Vector electricfield_ = system->environment->externalfield->electricfield;
    if (system->environment->externalfield->omega != 0.0)
    {
        for (int i = 0; i < prop->numAngles; i++)
        {
            double phi = i * 2 * M_PI / prop->numAngles;
            double Ex_ = electricfield_.norm() * cos(phi);
            double Ey_ = electricfield_.norm() * sin(phi);
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
        average_energy = - system->particles[1]->dipolemoment.dot(system->particles[1]->dipolemoment) * pow((1.0 / (system->particles[1]->state->r - system->particles[0]->state->r).norm()), 3) / 4.0;
        n = 1;
    }
    system->environment->externalfield->electricfield = electricfield_;

    return average_energy / n;
}
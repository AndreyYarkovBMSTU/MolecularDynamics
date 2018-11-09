#include "interaction.h"

double Interaction::getGradEnergy(int _nParticle, int _nCoord)
{
    system->particles[_nParticle]->state->r(_nCoord) = r_(_nCoord) + dr;
//    U_ = system->getInteractionEnergy(_nParticle).energy + system->getSelfEnergy(_nParticle).energy;
    U_ = system->getSelfEnergy(_nParticle).energy;
    system->particles[_nParticle]->state->r(_nCoord) = r_(_nCoord) - dr;
//    U_ = U_ - (system->getInteractionEnergy(_nParticle).energy + system->getSelfEnergy(_nParticle).energy);
    U_ = U_ - system->getSelfEnergy(_nParticle).energy;
    system->particles[_nParticle]->state->r(_nCoord) = r_(_nCoord);

    return U_ / (2 * dr);
}

double Interaction::getGradAverageEnergy(int _nParticle, int _nCoord)
{
    system->particles[_nParticle]->state->r(_nCoord) = r_(_nCoord) + dr;
//    U_ = (getAverageEnergy(_nParticle, "interaction") + getAverageEnergy(_nParticle, "self")) * prop->numAngles;
    U_ = getAverageEnergy(_nParticle, "self") * prop->numAngles;
    system->particles[_nParticle]->state->r(_nCoord) = r_(_nCoord) - dr;
//    U_ = U_ - (getAverageEnergy(_nParticle, "interaction") + getAverageEnergy(_nParticle, "self")) * prop->numAngles;
    U_ = U_ - getAverageEnergy(_nParticle, "self") * prop->numAngles;
    system->particles[_nParticle]->state->r(_nCoord) = r_(_nCoord);

    return U_ / (2 * dr);
}

Vector Interaction::getElectricForce(int _nParticle)
{
    if (accuracyEnergy == "average")
    {
        r_ = system->particles[_nParticle]->getCoordinate();

//        return Vector(- getGradAverageEnergy(_nParticle, 0), - getGradAverageEnergy(_nParticle, 1), 0.0);
        return (Vector(- getGradAverageEnergy(_nParticle, 0), - getGradAverageEnergy(_nParticle, 1), 0.0) + system->getForce(_nParticle));
//        return system->getForce(_nParticle);
    }
    if (accuracyEnergy == "exact")
    {
        method->setDipoleMoment();
        r_ = system->particles[_nParticle]->getCoordinate();

//        return Vector(- getGradEnergy(_nParticle, 0), - getGradEnergy(_nParticle, 1), 0.0);
        return (Vector(- getGradEnergy(_nParticle, 0), - getGradEnergy(_nParticle, 1), 0.0) + system->getForce(_nParticle));
//        return system->getForce(_nParticle);
    }
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
    path0 = prop->path + "deltaenergy1/" + "Self" + prop->filetype;
    out0.open(path0, std::ios::out | std::ios::binary);
//    out0 << system->dipolemoment0.dot(system->dipolemoment0) / (4 * pow(2 * system->particles[0]->radius, 3));
    out0 << system->dipolemoment0_.dot(system->dipolemoment0_) / (4 * pow(2, 3));
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

        path = prop->path + "deltaenergy1/" + _energytype + "Average" + std::to_string(i) + prop->filetype;
        out.open(path, std::ios::out | std::ios::binary);
        out << system->particles[1]->getCoordinate()(0) / 2 << " " << getAverageEnergy(0, _energytype);
        out.close();

//        if (system->particles[1]->getCoordinate()(0) / (2 * system->particles[0]->radius) >= 10.0)
        if (system->particles[1]->getCoordinate()(0) / 2 >= 10.0)
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
            system->environment->externalfield->electricfield = Vector(electricfield_.norm() * cos(phi), electricfield_.norm() * sin(phi), 0.0);
            method->setDipoleMoment();

            path = prop->path + "energy/" + _energytype + std::to_string(i) + prop->filetype;
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

    average_energy = 0;
    if (system->environment->externalfield->omega != 0.0)
    {
        for (int i = 0; i < prop->numAngles; i++)
        {
            phi = i * M_PI / prop->numAngles;
//            amplitude = 1.0 - cos(2.0 * phi + M_PI) / 3.0;
            amplitude = 1.0;

            system->environment->externalfield->electricfield = Vector(amplitude * electricfield_.norm() * cos(phi), amplitude * electricfield_.norm() * sin(phi), 0.0);
            system->setProperties();

            method->setDipoleMoment();
            if (_energytype == "induction")
            {
                average_energy += system->getInductionEnergy(_nParticle).energy;
            }
            else if (_energytype == "interaction")
            {
                average_energy += system->getInteractionEnergy(_nParticle).energy;
            }
            else if (_energytype == "self")
            {
                average_energy += system->getSelfEnergy(_nParticle).energy;
            }
            else if (_energytype == "ipl3")
            {
                system->setProperties();
                average_energy += - system->dipolemoment0_.dot(system->dipolemoment0_) * pow((1.0 / (system->particles[1]->getCoordinate() - system->particles[0]->getCoordinate()).norm()), 3) / 4;
            }
        }
        numpoints = prop->numAngles;
    }
    system->environment->externalfield->electricfield = electricfield_;

    return average_energy / numpoints;
}

void Interaction::recordPlane(std::string _pathinParticles, std::string _pathinPoints, std::string _pathoutPot, std::string _pathoutEl)
{
    in.open(_pathinParticles, std::ios::in);
    count = 0;
    if (in.is_open())
    {
        char *str = new char [1024];
        while (!in.eof())
            {
                in.getline(str, 1024, '\n');
                count++;
            }
        delete[] str;
        in.close();

        std::cout << "Succes Particles" << std::endl;
     }
    else
    {
        std::cout << "Error Particles" << std::endl;
    }

    numParticles = count - 1;
    states.resize(numParticles);
    cluster.resize(numParticles);

    in.open(_pathinParticles, std::ios::in);
    if (in.is_open())
    {
        m = 3;
        double **Part;
        Part = new double*[numParticles];

        for (int i = 0; i < numParticles; i++)
        {
            Part[i] = new double[m];
            for (int j = 0; j < m; j++)
            {
                in >> Part[i][j];
            }
        }

        for (int i = 0; i < numParticles; i++)
        {
            states[i] = new State(i + 1, Vector(Part[i][0], Part[i][1], 0.0), Vector(0.0, 0.0, 0.0));
            cluster[i] = new Dipoloid(states[i], prop->particlematerial, prop->obj0);
            std::cout << "Particle " << i+1 << " is created" << std::endl;
        }

        externalfield0 = new DirectedField(prop, Vector(1.0, 0.0, 0.0));
        environment0 = new Environment(prop->solvent, externalfield0);
        system0 = new ParticleSystem(environment0, prop);
        for (int i = 0; i < numParticles; i++)
        {
            system0->setParticle(cluster[i]);
        }

       method0 = new SelfConsistentDipoles(system0, prop);
       method0->setDipoleMoment();

        for (int i = 0; i < numParticles; i++)
        {
            delete[] Part[i];
        }
        delete[] Part;

        in.close();
    }

    in.open(_pathinPoints, std::ios::in);
    count = 0;
    if (in.is_open())
    {
        char *str = new char [1024];
        while (!in.eof())
            {
                in.getline(str, 1024, '\n');
                count++;
                std::cout << "Point " << count << std::endl;
            }
        delete[] str;
        in.close();

        std::cout << "Succes Points" << std::endl;
     }
    else
    {
        std::cout << "Error Points" << std::endl;
    }
    numpoints = count - 1;;
    m = 3;
    points.resize(numpoints);

    in.open(_pathinPoints, std::ios::in);
    if (in.is_open())
    {
        double **Cent;
        Cent = new double*[numpoints];
        for (int i = 0; i < numpoints; i++)
        {
            Cent[i] = new double[m];
        }

        for (int i = 0; i < numpoints; i++)
        {
            for (int j = 0; j < m; j++)
            {
                in >> Cent[i][j];
            }
            points[i] = Vector(Cent[i][0], Cent[i][1], Cent[i][2]);
        }

        Potentials.resize(numpoints, 1);
        Potentials.setZero(numpoints, 1);
        ElectricField.resize(numpoints, 3);
        ElectricField.setZero(numpoints, 3);

        for (int i = 0; i < numpoints; i++)
        {
            for (int nParticle = 0; nParticle < numParticles; nParticle++)
            {
                Potentials(i, 0) += cluster[nParticle]->getPotential(points[i]);
                ElectricField(i, 0) += cluster[nParticle]->getElectricField(points[i])(0);
                ElectricField(i, 1) += cluster[nParticle]->getElectricField(points[i])(1);
                ElectricField(i, 2) += cluster[nParticle]->getElectricField(points[i])(2);
            }
            std::cout << "N = " << i << " Potential = " <<Potentials(i, 0) << std::endl;
        }

        out.open(_pathoutPot, std::ios::out | std::ios::binary);
        out << Potentials;
        out.close();

        out.open(_pathoutEl, std::ios::out | std::ios::binary);
        out << ElectricField;
        out.close();

        for (int i = 0; i < numpoints; i++)
        {
            delete[] Cent[i];
        }
        delete[] Cent;

        in.close();
    }
}

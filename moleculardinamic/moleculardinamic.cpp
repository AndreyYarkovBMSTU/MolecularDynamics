#include "moleculardinamic.h"

void MolecularDinamic::dump()
{
    for (auto &iParticle : system->particles)
    {
        std::cout << iParticle->state->number << "\t";
        std::cout << iParticle->name << std::endl;
    }
}

void MolecularDinamic::computeKoef()
{
    koef_demping = -18.0 / system->reinolds;
//    koef_demping = 0;
    koef_LenJon = prop->koef_LenJon / (2 * prop->k * prop->temperature);
//    koef_LenJon = 0;
    koef_randForce = prop->koef_Brown * pow(18 / system->reinolds, 0.5);
//    koef_randForce = pow(2 * system->friction * prop->k * prop->temperature, 0.5) / system->particles[0]->mass;
    koef_dipole = system->dipolemoment0.dot(system->dipolemoment0) / (2.0 * prop->k * prop->temperature * pow(2 * system->particles[0]->radius, 3));
//    koef_dipole = 0;
}

void MolecularDinamic::computer(int _nParticle, int _nFrame)
{
    computeKoef();

    r_.resize(system->numParticles);
    rp_.resize(system->numParticles);
    R_.resize(system->numParticles, 3);

    for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
    {
        rp_[nParticle] = system->particles[nParticle]->getCoordinate();
        r_[nParticle] = system->particles[nParticle]->getCoordinate();
    }

    rad = 0.0;
    srand(time(NULL));
    for (int frame = 0; frame < _nFrame; frame++)
    {
        if (accuracyEnergy != "average")
        {
            system->environment->externalfield->electricfield = system->environment->externalfield->getElectricField(frame);
        }

        for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
        {
            R_(nParticle, 0) = system->particles[nParticle]->getCoordinate()(0);
            R_(nParticle, 1) = system->particles[nParticle]->getCoordinate()(1);
            R_(nParticle, 2) = system->particles[nParticle]->getCoordinate()(2);
        }

        if (nameThermostat == "langevin")
        {
            setCoordinateLangevin();
            setVelocityLangevin();
        }
        else if (nameThermostat == "brownian")
        {
            setVelocityBrownian();
            setCoordinateBrownian();
        }
        else
        {
            std::cout << "NameThermostatError" << std::endl;
        }

        rad /= system->numParticles;
    }
    //        std::cout << system->particles[_nParticle]->state << std::endl;
    //std::cout << prop->mass << std::endl;
    //std::cout << thermostat->frictionkoef << std::endl;
    //std::cout << 4 * system->diffusion * prop->timestep / pow(2 * system->particles[0]->radius, 2) << std::endl;
    //std::cout << pow(rad / _nFrame, 2) << std::endl;
}

void MolecularDinamic::record(std::string _path)
{
    computeKoef();

    r_.resize(system->numParticles);            // Радиус-вектор частицы в настоящий момент времени
    rp_.resize(system->numParticles);           // Радиус-вектор частицы в предыдущий момент времени
    R_.resize(system->numParticles, 3);

    for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
    {
        rp_[nParticle] = system->particles[nParticle]->getCoordinate();
        r_[nParticle] = system->particles[nParticle]->getCoordinate();
    }

//    path = _path + std::to_string(0) + prop->filetype;
    path = _path + prop->filetype;
    out.open(path, std::ios::out | std::ios::binary);
    for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
    {
        out << system->particles[nParticle]->state;
    }

    rad = 0.0;
    srand(time(NULL));
    for (int i = 1; i < numFrames; i++)
    {
//        k = (0.016 / prop->timestep * 1e-4);            // проверка для записи определённых файлов
        k = 1;

        if (accuracyEnergy != "average")
        {
            system->environment->externalfield->electricfield = system->environment->externalfield->getElectricField(i, tau);
        }

//        path = _path + std::to_string(i) + prop->filetype;
//        out.open(path, std::ios::out | std::ios::binary);

        for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
        {
            R_(nParticle, 0) = system->particles[nParticle]->getCoordinate()(0);
            R_(nParticle, 1) = system->particles[nParticle]->getCoordinate()(1);
            R_(nParticle, 2) = system->particles[nParticle]->getCoordinate()(2);
        }

        if (nameThermostat == "langevin")
        {
            setCoordinateLangevin();
            setVelocityLangevin();
        }
        else if (nameThermostat == "brownian")
        {
            setVelocityBrownian();
            setCoordinateBrownian();
        }
        else
        {
            std::cout << "NameThermostatError" << std::endl;
        }

        rad /= system->numParticles;

//        for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
//        {
//            if (bool(i % k == 0))
//            {
//                out << system->particles[nParticle]->state;
//            }
//            else
//            {
//                path = _path + std::to_string(i) + prop->filetype;
//                c = path.c_str();
//                remove(c);
//            }
//        }
//        out.close();
        for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
        {
            out << system->particles[nParticle]->state;
        }

        if (i == numFrames / 4.0)
        {
            std::cout << "0.25" << std::endl;
        }
        if (i == numFrames / 2.0)
        {
            std::cout << "0.5" << std::endl;
        }
        if (i == 3.0 * numFrames / 4.0)
        {
            std::cout << "0.75" << std::endl;
        }
        std::cout << "N = " << i << "  v1 = " << system->particles[0]->getVelocity().norm() << "  v2 = " << system->particles[1]->getVelocity().norm() <<"  v3 = " << system->particles[2]->getVelocity().norm() << std::endl;
        if (system->particles[0]->getVelocity().norm() > 100 || system->particles[1]->getVelocity().norm() > 100 || system->particles[2]->getVelocity().norm() > 100)
        {
            break;
        }
    }
    out.close();
}

void MolecularDinamic::recordtest(std::string _path)
{
    computeKoef();

    r_.resize(system->numParticles);            // Радиус-вектор частицы в настоящий момент времени
    rp_.resize(system->numParticles);           // Радиус-вектор частицы в предыдущий момент времени
    phi0.resize(system->numParticles);

    R_.resize(system->numParticles, 3);

    for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
    {
        rp_[nParticle] = system->particles[nParticle]->getCoordinate();
        r_[nParticle] = system->particles[nParticle]->getCoordinate();
        phi0[nParticle] = atan(system->particles[nParticle]->getCoordinate()(1) / system->particles[nParticle]->getCoordinate()(0));
        if (system->particles[nParticle]->getCoordinate()(1) == 0 && system->particles[nParticle]->getCoordinate()(0) < 0)
        {
            phi0[nParticle] = M_PI;
        }
    }

    path = _path + prop->filetype;
    out.open(path, std::ios::out | std::ios::binary);
    for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
    {
        out << system->particles[nParticle]->state;
    }

    double omega = 0.0563218;
    srand(time(NULL));
    for (int i = 1; i < numFrames; i++)
    {
        k = 1;

        for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
        {
            system->particles[nParticle]->setCoordinate(Vector(system->particles[nParticle]->getCoordinate().norm() * cos(phi0[nParticle] + omega * tau * i), system->particles[nParticle]->getCoordinate().norm() * sin(phi0[nParticle] + omega * tau * i), 0.0));
            system->particles[nParticle]->setVelocity(system->particles[nParticle]->getCoordinate() - r_[nParticle]);
            r_[nParticle] = system->particles[nParticle]->getCoordinate();
        }

        for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
        {
            out << system->particles[nParticle]->state;
        }
    }
    out.close();
}

void MolecularDinamic::recordVMD(std::string _path)
{
    computeKoef();

    r_.resize(system->numParticles);            // Радиус-вектор частицы в настоящий момент времени
    rp_.resize(system->numParticles);           // Радиус-вектор частицы в предыдущий момент времени
    R_.resize(system->numParticles, 3);

    for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
    {
        rp_[nParticle] = system->particles[nParticle]->getCoordinate();
        r_[nParticle] = system->particles[nParticle]->getCoordinate();
    }

    path = _path + ".lammpstrj";
    out.open(path, std::ios::out | std::ios::binary);
    outputline = "ITEM: TIMESTEP\n" + std::to_string(0) + "\nITEM: NUMBER OF ATOMS\n" + std::to_string(system->numParticles) + "\nITEM: BOX BOUNDS pp pp\n-50 50\n-50 50\n0.0 0.0\nITEM: ATOMS id x y\n";
    out << outputline;
    for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
    {
        outputline = std::to_string(system->particles[nParticle]->state->number) + " " + std::to_string(system->particles[nParticle]->state->r(0)) + " " + std::to_string(system->particles[nParticle]->state->r(1)) + "\n";
        out << outputline;
    }

    rad = 0.0;
    srand(time(NULL));
    for (int i = 1; i < numFrames; i++)
    {
//        k = (0.016 / prop->timestep * 1e-4);            // проверка для записи определённых файлов
        k = 20;

        if (bool(i % k == 0))
        {
            outputline = "ITEM: TIMESTEP\n" + std::to_string(i) + "\nITEM: NUMBER OF ATOMS\n" + std::to_string(system->numParticles) + "\nITEM: BOX BOUNDS pp pp pp\n-50 50\n-50 50\n0.0 0.0\nITEM: ATOMS id x y\n";
            out << outputline;
        }

        if (accuracyEnergy != "average")
        {
            system->environment->externalfield->electricfield = system->environment->externalfield->getElectricField(i);
        }

        for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
        {
            R_(nParticle, 0) = system->particles[nParticle]->getCoordinate()(0);
            R_(nParticle, 1) = system->particles[nParticle]->getCoordinate()(1);
            R_(nParticle, 2) = system->particles[nParticle]->getCoordinate()(2);
        }

        if (nameThermostat == "langevin")
        {
            setCoordinateLangevin();
            setVelocityLangevin();
        }
        else if (nameThermostat == "brownian")
        {
            setVelocityBrownian();
            setCoordinateBrownian();
        }
        else
        {
            std::cout << "NameThermostatError" << std::endl;
        }

        rad /= system->numParticles;

        for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
        {
            if (bool(i % k == 0))
            {
                outputline = std::to_string(system->particles[nParticle]->state->number) + " " + std::to_string(system->particles[nParticle]->state->r(0)) + " " + std::to_string(system->particles[nParticle]->state->r(1)) + "\n";
                out << outputline;
            }
        }

        if (i == numFrames / 4.0)
        {
            std::cout << "0.25" << std::endl;
        }
        if (i == numFrames / 2.0)
        {
            std::cout << "0.5" << std::endl;
        }
        if (i == 3.0 * numFrames / 4.0)
        {
            std::cout << "0.75" << std::endl;
        }
        std::cout << "N = " << i << "  v1 = " << system->particles[0]->getVelocity().norm() << "  v2 = " << system->particles[1]->getVelocity().norm() <<"  v3 = " << system->particles[2]->getVelocity().norm() << std::endl;
        if (system->particles[0]->getVelocity().norm() > 100 || system->particles[1]->getVelocity().norm() > 100 || system->particles[2]->getVelocity().norm() > 100)
        {
            break;
        }
    }
    out.close();
}

void MolecularDinamic::recordmove(int _numfile, int _nFrame)
{
    computeKoef();

    numParticles = system->numParticles;
    system->numParticles = 1;

    path = prop->path + "move/" + std::to_string(_numfile) + prop->filetype;
    out.open(path, std::ios::out | std::ios::binary);
    outputline = std::to_string(0) + "  " + std::to_string(0) + "\n";

    r_.resize(system->numParticles);
    rp_.resize(system->numParticles);
    R_.resize(system->numParticles, 3);

    for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
    {
        rp_[nParticle] = system->particles[nParticle]->getCoordinate();
        r_[nParticle] = system->particles[nParticle]->getCoordinate();
    }

    rad = 0.0;
    srand(time(NULL));
    for (int frame = 0; frame < _nFrame; frame++)
    {
        for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
        {
            R_(nParticle, 0) = system->particles[nParticle]->getCoordinate()(0);
            R_(nParticle, 1) = system->particles[nParticle]->getCoordinate()(1);
            R_(nParticle, 2) = system->particles[nParticle]->getCoordinate()(2);
        }

        if (nameThermostat == "langevin")
        {
            setCoordinateLangevin();
            setVelocityLangevin();
        }
        else if (nameThermostat == "brownian")
        {
            setVelocityBrownian();
            setCoordinateBrownian();
        }
        else
        {
            std::cout << "NameThermostatError" << std::endl;
        }

        if (frame % 100 == 0)
        {
            std::cout << outputline << std::endl;
        }
    }
    out << outputline;
    out.close();

    system->numParticles = numParticles;
}

void MolecularDinamic::setVelocityLangevin()
{
    for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
    {
        f_LenJon = Vector(0.0, 0.0, 0.0);
        f_dipole = Vector(0.0, 0.0, 0.0);
        for (int j = 0; j < system->numParticles; j++)
        {
            if (j != nParticle)
            {
                f_LenJon += - potential->getGradPotential(Vector(R_(nParticle, 0), R_(nParticle, 1), R_(nParticle, 2)),
                                                          Vector(R_(j, 0), R_(j, 1), R_(j, 2)));
//                                            if ((system->particles[nParticle]->state->r - system->particles[j]->state->r).norm() > 2)
//                                            {
//                                                f_dipole += interaction->getElectricForce(nParticle);
//                                            }
            }
        }
        a_ = (koef_demping * system->particles[nParticle]->getVelocity() + koef_LenJon * f_LenJon + koef_randForce * thermostat->getForce());

        system->particles[nParticle]->setVelocity(numEq->getVelocity(system->particles[nParticle]->getVelocity(), a_, a));

        rp_[nParticle] = r_[nParticle];

        rad+=(system->particles[nParticle]->getCoordinate() - r_[nParticle]).dot((system->particles[nParticle]->getCoordinate() - r_[nParticle]));
    }
}

void MolecularDinamic::setCoordinateLangevin()
{
    for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
    {
        f_LenJon = Vector(0.0, 0.0, 0.0);
        f_dipole = Vector(0.0, 0.0, 0.0);
        for (int j = 0; j < system->numParticles; j++)
        {
            if (j != nParticle)
            {
                f_LenJon += - potential->getGradPotential(Vector(R_(nParticle, 0), R_(nParticle, 1), R_(nParticle, 2)),
                                                          Vector(R_(j, 0), R_(j, 1), R_(j, 2)));
            }
        }
        f_dipole = interaction->getElectricForce(nParticle);

        a = (koef_demping * system->particles[nParticle]->getVelocity() + koef_LenJon * f_LenJon + koef_randForce * thermostat->getForce() + koef_dipole * f_dipole);

        r_[nParticle] = system->particles[nParticle]->getCoordinate();

        system->particles[nParticle]->setCoordinate(numEq->getCoordinates(system->particles[nParticle]->getCoordinate(), rp_[nParticle], a));
    }
}

void MolecularDinamic::setVelocityBrownian()
{
    for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
    {
        f_LenJon = Vector(0.0, 0.0, 0.0);
        f_dipole = Vector(0.0, 0.0, 0.0);
        for (int j = 0; j < system->numParticles; j++)
        {
            if (j != nParticle)
            {
                f_LenJon += - potential->getGradPotential(Vector(R_(nParticle, 0), R_(nParticle, 1), R_(nParticle, 2)),
                                                          Vector(R_(j, 0), R_(j, 1), R_(j, 2)));
            }
        }

        f_dipole = interaction->getElectricForce(nParticle);

//        v = system->reinolds * (prop->koef_LenJon / (36 * prop->k * prop->temperature)) * f_LenJon + thermostat->getForce() * pow(system->reinolds / 18, 0.5) + koef_dipole * f_dipole * system->reinolds / 18;
        v = (koef_LenJon * f_LenJon + koef_randForce * thermostat->getForce() + koef_dipole * f_dipole) * (system->reinolds / 18);

        system->particles[nParticle]->setVelocity(v);

        rp_[nParticle] = r_[nParticle];
        r_[nParticle] = system->particles[nParticle]->getCoordinate();
    }
}

void MolecularDinamic::setCoordinateBrownian()
{
    for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
    {
        system->particles[nParticle]->setCoordinate(system->particles[nParticle]->getCoordinate() + system->particles[nParticle]->getVelocity() * t0);

        rad+=(system->particles[nParticle]->getCoordinate() - r_[nParticle]).dot((system->particles[nParticle]->getCoordinate() - r_[nParticle]));
    }
}

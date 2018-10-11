#include "moleculardinamic.h"

void MolecularDinamic::dump()
{
    for (auto &iParticle : system->particles)
    {
        std::cout << iParticle->state->number << "\t";
        std::cout << iParticle->name << std::endl;
    }
}

    void MolecularDinamic::computer(int _nParticle, int _nFrame)
    {
        r_.resize(system->numParticles);
        rp_.resize(system->numParticles);
        R_.resize(system->numParticles, 3);

        for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
        {
            rp_[nParticle] = system->particles[nParticle]->getCoordinate();
        }

        rad = 0.0;
        srand(time(NULL));
        for (int frame = 0; frame < _nFrame; frame++)
        {
            system->environment->externalfield->electricfield = system->environment->externalfield->getElectricField(frame);
            for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
            {
                R_(nParticle, 0) = system->particles[nParticle]->getCoordinate()(0);
                R_(nParticle, 1) = system->particles[nParticle]->getCoordinate()(1);
                R_(nParticle, 2) = system->particles[nParticle]->getCoordinate()(2);
            }

            f = Vector(0.0, 0.0, 0.0);
            // Расчёт координат частицы
            for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
            {
                f_LenJon = Vector(0.0, 0.0, 0.0);
                f_dipole = Vector(0.0, 0.0, 0.0);
                for (int j = 0; j < system->numParticles; j++)
                {
                    if (j != nParticle)
                    {
                        f_LenJon += - potential->getGradPotential(Vector(R_(nParticle, 0), R_(nParticle, 1), R_(nParticle, 2)),
                                                                  Vector(R_(j, 0), R_(j, 1), R_(j, 2))) * prop->koef_LenJon;
    //                        if ((system->particles[nParticle]->state->r - system->particles[j]->state->r).norm() > 2 * system->particles[0]->radius)
    //                        {
    //                            f_dipole += interaction->getElectricForce(nParticle);
    //                        }
                    }
                }
                f = ((-18.0 * system->particles[nParticle]->getVelocity() / system->reinolds) - (1.0 / (2 * prop->k * prop->temperature)) * f_LenJon + pow(system->reinolds, 0.5) * thermostat->getForce()) * system->particles[nParticle]->mass;
                //Vector f = f_therm + f_LenJon + f_dipole;
                r_[nParticle] = system->particles[nParticle]->getCoordinate();
                system->particles[nParticle]->setCoordinate(numEq->getCoordinates(system->particles[nParticle]->getCoordinate(), rp_[nParticle], f, system->particles[nParticle]->mass));
            }
            for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
            {
                R_(nParticle, 0) = system->particles[nParticle]->getCoordinate()(0);
                R_(nParticle, 1) = system->particles[nParticle]->getCoordinate()(1);
                R_(nParticle, 2) = system->particles[nParticle]->getCoordinate()(2);
            }

            // Расчёт скорости частицы
            for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
            {
                f_LenJon = Vector(0.0, 0.0, 0.0);
                f_dipole = Vector(0.0, 0.0, 0.0);
                for (int j = 0; j < system->numParticles; j++)
                {
                    if (j != nParticle)
                    {
                        f_LenJon += - potential->getGradPotential(Vector(R_(nParticle, 0), R_(nParticle, 1), R_(nParticle, 2)),
                                                                  Vector(R_(j, 0), R_(j, 1), R_(j, 2))) * prop->koef_LenJon;
                        //                        if ((system->particles[nParticle]->state->r - system->particles[j]->state->r).norm() > 2 * system->particles[0]->radius)
                        //                        {
                        //                            f_dipole_for_numeq += interaction->getElectricForce(nParticle);
                        //                        }
                    }
                }
                f_ = ((-18.0 * system->particles[nParticle]->getVelocity() / system->reinolds) - (1.0 / (2 * prop->k * prop->temperature)) * f_LenJon + pow(system->reinolds, 0.5) * thermostat->getForce()) * system->particles[0]->mass;
                system->particles[nParticle]->setVelocity(numEq->getVelocity(system->particles[nParticle]->getVelocity(), f_, f, system->particles[nParticle]->mass));
                rp_[nParticle] = r_[nParticle];

                rad+=(system->particles[nParticle]->getCoordinate() - r_[nParticle]).dot(system->particles[nParticle]->getCoordinate() - r_[nParticle]) * system->particles[0]->radius;
                //std::cout << f << std::endl;
                //std::cout << nParticle << ": " << f_electrictest << std::endl;
                //std::cout << nParticle << ": " << f_electric << std::endl;
            }
                rad /= system->numParticles;
        }
//        std::cout << system->particles[_nParticle]->state << std::endl;
        //std::cout << prop->mass << std::endl;
        //std::cout << thermostat->frictionkoef << std::endl;
        std::cout << sqrt(4 * system->diffusion * t0) << std::endl;
        std::cout << sqrt(rad / _nFrame) << std::endl;
    }

void MolecularDinamic::record()
{
    r_.resize(system->numParticles);            // Радиус-вектор частицы в настоящий момент времени
    rp_.resize(system->numParticles);           // Радиус-вектор частицы в предыдущий момент времени
    R_.resize(system->numParticles, 3);

    for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
    {
        rp_[nParticle] = system->particles[nParticle]->getCoordinate();
    }

    path = prop->path + std::to_string(0) + prop->filetype;
    out.open(path, std::ios::out | std::ios::binary);
    for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
    {
        out << system->particles[nParticle]->state;
    }
    out.close();

    srand(time(NULL));
    for (int i = 1; i < prop->numFrames; i++)
    {
//        k = (0.016 / prop->timestep * 1e-4);            // проверка для записи определённых файлов
//        k = 1;
        k = 160;
        system->environment->externalfield->electricfield = system->environment->externalfield->getElectricField(i);

        path = prop->path + std::to_string(i) + prop->filetype;
        out.open(path, std::ios::out | std::ios::binary);

        for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
        {
            R_(nParticle, 0) = system->particles[nParticle]->getCoordinate()(0);
            R_(nParticle, 1) = system->particles[nParticle]->getCoordinate()(1);
            R_(nParticle, 2) = system->particles[nParticle]->getCoordinate()(2);
        }

        f = Vector(0.0, 0.0, 0.0);
        // Расчёт координат частицы
        for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
        {
            f_LenJon = Vector(0.0, 0.0, 0.0);
            f_dipole = Vector(0.0, 0.0, 0.0);
            for (int j = 0; j < system->numParticles; j++)
            {
                if (j != nParticle)
                {
                    f_LenJon += - potential->getGradPotential(Vector(R_(nParticle, 0), R_(nParticle, 1), R_(nParticle, 2)),
                                                              Vector(R_(j, 0), R_(j, 1), R_(j, 2))) * prop->koef_LenJon;
//                        if ((system->particles[nParticle]->state->r - system->particles[j]->state->r).norm() > 2)
//                        {
//                            f_dipole += interaction->getElectricForce(nParticle);
//                        }
                }
            }
            f = ((-18.0 * system->particles[nParticle]->getVelocity() / system->reinolds) + (1.0 / (2 * prop->k * prop->temperature)) * f_LenJon + pow(system->reinolds, 0.5) * thermostat->getForce()) * system->particles[0]->mass;
            //Vector f = f_therm + f_LenJon + f_dipole;
            if (i == 160)
            {
                std::cout << thermostat->getForce() << std::endl;
            }
            r_[nParticle] = system->particles[nParticle]->getCoordinate();
            system->particles[nParticle]->setCoordinate(numEq->getCoordinates(system->particles[nParticle]->getCoordinate(), rp_[nParticle], f, system->particles[nParticle]->mass));
        }
        for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
        {
            R_(nParticle, 0) = system->particles[nParticle]->getCoordinate()(0);
            R_(nParticle, 1) = system->particles[nParticle]->getCoordinate()(1);
            R_(nParticle, 2) = system->particles[nParticle]->getCoordinate()(2);
        }

        // Расчёт скорости частицы
        for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
        {
            f_LenJon = Vector(0.0, 0.0, 0.0);
            f_dipole = Vector(0.0, 0.0, 0.0);
            for (int j = 0; j < system->numParticles; j++)
            {
                if (j != nParticle)
                {
                    f_LenJon += - potential->getGradPotential(Vector(R_(nParticle, 0), R_(nParticle, 1), R_(nParticle, 2)),
                                                              Vector(R_(j, 0), R_(j, 1), R_(j, 2))) * prop->koef_LenJon;
//                                            if ((system->particles[nParticle]->state->r - system->particles[j]->state->r).norm() > 2)
//                                            {
//                                                f_dipole += interaction->getElectricForce(nParticle);
//                                            }
                }
            }
            f_ = ((-18.0 * system->particles[nParticle]->getVelocity() / system->reinolds) - (1.0 / (2 * prop->k * prop->temperature)) * f_LenJon + pow(system->reinolds, 0.5) * thermostat->getForce()) * system->particles[0]->mass;
            system->particles[nParticle]->setVelocity(numEq->getVelocity(system->particles[nParticle]->getVelocity(), f_, f, system->particles[nParticle]->mass));
            rp_[nParticle] = r_[nParticle];

            if (bool(i % k == 0))
            {
                out << system->particles[nParticle]->state;
            }
            else
            {
                path = prop->path + std::to_string(i) + prop->filetype;
                c = path.c_str();
                remove(c);
            }
        }
            out.close();
    }
}

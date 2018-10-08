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
    //                Vector f_therm = thermostat->getForce(system->particles[nParticle]->getVelocity());
                f_therm = Vector(0.0, 0.0, 0.0);
                f_LenJon = Vector(0.0, 0.0, 0.0);
                f_dipole = Vector(0.0, 0.0, 0.0);
                for (int j = 0; j < system->numParticles; j++)
                {
                    if (j != nParticle)
                    {
                        f_LenJon += - potential->getGradPotential(Vector(R_(nParticle, 0), R_(nParticle, 1), R_(nParticle, 2)),
                                                                  Vector(R_(j, 0), R_(j, 1), R_(j, 2))) * prop->koef_LenJon;
    //                        if ((system->particles[nParticle]->state->r - system->particles[j]->state->r).norm() > 2 * prop->radius)
    //                        {
    //                            f_dipole += interaction->getElectricForce(nParticle);
    //                        }
                    }
                }
                f = - prop->prandtl * (system->particles[nParticle]->getVelocity() + (f_dipole / prop->prandtlEl) - pow(2.0, 0.5) * prop->knudsen * thermostat->getForce()) * prop->mass;
                //Vector f = f_therm + f_LenJon + f_dipole;
                r_[nParticle] = system->particles[nParticle]->getCoordinate();
                system->particles[nParticle]->setCoordinate(numEq->getCoordinates(system->particles[nParticle]->getCoordinate(), rp_[nParticle], f, prop->mass));

                R_(nParticle, 0) = system->particles[nParticle]->getCoordinate()(0);
                R_(nParticle, 1) = system->particles[nParticle]->getCoordinate()(1);
                R_(nParticle, 2) = system->particles[nParticle]->getCoordinate()(2);
            }

            // Расчёт скорости частицы
            for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
            {
    //                Vector f_therm_for_verle = thermostat->getForce(system->particles[nParticle]->state->v);
                f_LenJon = Vector(0.0, 0.0, 0.0);
                f_dipole = Vector(0.0, 0.0, 0.0);
                for (int j = 0; j < system->numParticles; j++)
                {
                    if (j != nParticle)
                    {
                        f_LenJon += - potential->getGradPotential(Vector(R_(nParticle, 0), R_(nParticle, 1), R_(nParticle, 2)),
                                                                  Vector(R_(j, 0), R_(j, 1), R_(j, 2))) * prop->koef_LenJon;
                        //                        if ((system->particles[nParticle]->state->r - system->particles[j]->state->r).norm() > 2 * prop->radius)
                        //                        {
                        //                            f_dipole_for_numeq += interaction->getElectricForce(nParticle);
                        //                        }
                    }
                }
                f_ = - prop->prandtl * (system->particles[nParticle]->getVelocity() + f_dipole / prop->prandtlEl - pow(2.0, 0.5) * prop->knudsen * thermostat->getForce()) * prop->mass;
                system->particles[nParticle]->setVelocity(numEq->getVelocity(system->particles[nParticle]->getVelocity(), f_, f, prop->mass));
                rp_[nParticle] = r_[nParticle];

                rad+=(system->particles[nParticle]->getCoordinate() - r_[nParticle]).dot(system->particles[nParticle]->getCoordinate() - r_[nParticle]);
                //std::cout << f << std::endl;
                //std::cout << nParticle << ": " << f_electrictest << std::endl;
                //std::cout << nParticle << ": " << f_electric << std::endl;
            }
                rad /= system->numParticles;
        }
//        std::cout << system->particles[_nParticle]->state << std::endl;
        //std::cout << prop->mass << std::endl;
        //std::cout << thermostat->frictionkoef << std::endl;
        //std::cout << sqrt(4 * thermostat->diffusionkoef * prop->timestep) << std::endl;
        std::cout << sqrt(rad / prop->numFrames) << std::endl;
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
//            k = (0.016 / prop->timestep * 1e-4);            // проверка для записи определённых файлов
        k = 1;
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
//                Vector f_therm = thermostat->getForce(system->particles[nParticle]->getVelocity());
            f_therm = Vector(0.0, 0.0, 0.0);
            f_LenJon = Vector(0.0, 0.0, 0.0);
            f_dipole = Vector(0.0, 0.0, 0.0);
            for (int j = 0; j < system->numParticles; j++)
            {
                if (j != nParticle)
                {
                    f_LenJon += - potential->getGradPotential(Vector(R_(nParticle, 0), R_(nParticle, 1), R_(nParticle, 2)),
                                                              Vector(R_(j, 0), R_(j, 1), R_(j, 2))) * prop->koef_LenJon;
//                        if ((system->particles[nParticle]->state->r - system->particles[j]->state->r).norm() > 2 * prop->radius)
//                        {
//                            f_dipole += interaction->getElectricForce(nParticle);
//                        }
                }
            }
            f = - prop->prandtl * (system->particles[nParticle]->getVelocity() + (f_dipole / prop->prandtlEl) - pow(2.0, 0.5) * prop->knudsen * thermostat->getForce()) * prop->mass;
            //Vector f = f_therm + f_LenJon + f_dipole;
            r_[nParticle] = system->particles[nParticle]->getCoordinate();
            system->particles[nParticle]->setCoordinate(numEq->getCoordinates(system->particles[nParticle]->getCoordinate(), rp_[nParticle], f, prop->mass));

            R_(nParticle, 0) = system->particles[nParticle]->getCoordinate()(0);
            R_(nParticle, 1) = system->particles[nParticle]->getCoordinate()(1);
            R_(nParticle, 2) = system->particles[nParticle]->getCoordinate()(2);
        }

        // Расчёт скорости частицы
        for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
        {
//                Vector f_therm_for_verle = thermostat->getForce(system->particles[nParticle]->state->v);
            f_LenJon = Vector(0.0, 0.0, 0.0);
            f_dipole = Vector(0.0, 0.0, 0.0);
            for (int j = 0; j < system->numParticles; j++)
            {
                if (j != nParticle)
                {
                    f_LenJon += - potential->getGradPotential(Vector(R_(nParticle, 0), R_(nParticle, 1), R_(nParticle, 2)),
                                                              Vector(R_(j, 0), R_(j, 1), R_(j, 2))) * prop->koef_LenJon;
                    //                        if ((system->particles[nParticle]->state->r - system->particles[j]->state->r).norm() > 2 * prop->radius)
                    //                        {
                    //                            f_dipole_for_numeq += interaction->getElectricForce(nParticle);
                    //                        }
                }
            }
            f_ = - prop->prandtl * (system->particles[nParticle]->getVelocity() + f_dipole / prop->prandtlEl - pow(2.0, 0.5) * prop->knudsen * thermostat->getForce()) * prop->mass;
            system->particles[nParticle]->setVelocity(numEq->getVelocity(system->particles[nParticle]->getVelocity(), f_, f, prop->mass));
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

#include "moleculardinamic.h"

void MolecularDinamic::dump()
{
    for (auto &iParticle : system->particles)
    {
        std::cout << iParticle->state->number << "\t";
        std::cout << iParticle->name << std::endl;
    }
}

//    void MolecularDinamic::computer(int _nParticle, int _nFrame)
//    {
//        Vector r_[system->numParticles];
//        Vector rp_[system->numParticles];

//        for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
//        {
//            rp_[nParticle] = system->particles[nParticle]->state->r;
//        }

//        double rad = 0.0;
//        srand(time(NULL));
//        for (int frame = 0; frame < _nFrame; frame++)
//        {
//            if (system->environment->externalfield->omega != 0.0)
//            {
//                system->environment->externalfield->electricfield = Vector(system->environment->externalfield->electricfield.norm() * cos(frame * system->environment->externalfield->omega * prop->timestep),
//                                                                           system->environment->externalfield->electricfield.norm() * sin(frame * system->environment->externalfield->omega * prop->timestep),
//                                                                           0.0);
//            }
//            for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
//            {
//                //Vector f_therm = Vector(0.0, 0.0, 0.0);
//                Vector f_LenJon = Vector(0.0, 0.0, 0.0);
//                Vector f_dipole = Vector(0.0, 0.0, 0.0);
////                for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
////                {
////                    if (nParticle != nParticle)
////                    {
////                        if ((system->particles[nParticle]->state->r - system->particles[nParticle]->state->r).norm() <= 2 * um())
////                        {
////                            f_solid += - potential->getGradPotential(system->particles[nParticle]->state->r, system->particles[nParticle]->state->r) * 1e-11;
////                        }
//////                        if ((system->particles[nParticle]->state->r - system->particles[i]->state->r).norm() > 2 * um())
//////                        {
//////                            f_electrictest += ipl3->getGradPotential(system->particles[nParticle]->state->r, system->particles[i]->state->r)* 1e-9;
//////                        }
////                    }
////                }
//                for (int i = 0; i < system->numParticles; i++)
//                {
//                    if (i != nParticle)
//                    {
//                        f_LenJon += - potential->getGradPotential(system->particles[nParticle]->state->r, system->particles[i]->state->r) * 1e-9;
////                        if ((system->particles[nParticle]->state->r - system->particles[i]->state->r).norm() > 2 * prop->radius)
////                        {
////                            f_dipole += interaction->getElectricForce(nParticle);
////                        }
//                    }
//                }
//                Vector f = f_therm + f_LenJon + f_dipole;
//                r_[nParticle] = system->particles[nParticle]->state->r;
//                system->particles[nParticle]->state->r = numEq->getCoordinates(system->particles[nParticle]->state->r, rp_[nParticle], f, prop->mass);

//                // Расчёт скорости частицы
//                Vector f_therm_for_verle = Vector(0.0, 0.0, 0.0);
//                Vector f_LenJon_for_verle = Vector(0.0, 0.0, 0.0);
//                Vector f_dipole_for_verle = Vector(0.0, 0.0, 0.0);
//                for (int i = 0; i < system->numParticles; i++)
//                {
//                    if (i != nParticle)
//                    {
//                        f_LenJon_for_verle += - potential->getGradPotential(system->particles[nParticle]->state->r, system->particles[i]->state->r) * 1e-9;
////                        if ((system->particles[nParticle]->state->r - system->particles[i]->state->r).norm() > 2 * prop->radius)
////                        {
////                            f_dipole_for_verle += interaction->getElectricForce(nParticle);
////                        }
//                    }
//                }
//                Vector f_for_verle = f_therm_for_verle + f_LenJon_for_verle + f_dipole_for_verle;
    //      ERROR!!!!!
/*
     system->particles[nParticle]->setVelocity(numEq->getVelocity(system->particles[nParticle]->state->v, f_for_verle, f, prop->mass));
*/
//                system->particles[nParticle]->state->v = numEq->getVelocity(system->particles[nParticle]->state->v, f_for_verle, f, prop->mass);
//                rp_[nParticle] = r_[nParticle];

//                rad+=(system->particles[nParticle]->state->r - r_[nParticle]).dot(system->particles[nParticle]->state->r - r_[nParticle]);
//                //std::cout << f << std::endl;
//                //std::cout << nParticle << ": " << f_electrictest << std::endl;
//                //std::cout << nParticle << ": " << f_electric << std::endl;
//            }
//                rad /= system->numParticles;
//        }
//        //std::cout << system->particles[_nParticle]->state << std::endl;
//        //std::cout << prop->mass << std::endl;
//        //std::cout << thermostat->frictionkoef << std::endl;
//        //std::cout << sqrt(4 * thermostat->diffusionkoef * prop->timestep) << std::endl;
//        std::cout << sqrt(rad / prop->numFrames) << std::endl;
//    }

void MolecularDinamic::record()
{
    Vector r_[system->numParticles];            // Радиус-вектор частицы в настоящий момент времени
    Vector rp_[system->numParticles];           // Радиус-вектор частицы в предыдущий момент времени

    for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
    {
        rp_[nParticle] = system->particles[nParticle]->state->r;
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

        Matrix rad_vectors(system->numParticles, 3);
        for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
        {
            rad_vectors(nParticle, 0) = system->particles[nParticle]->getCoordinate()(0);
            rad_vectors(nParticle, 1) = system->particles[nParticle]->getCoordinate()(1);
            rad_vectors(nParticle, 2) = system->particles[nParticle]->getCoordinate()(2);
        }

        f = Vector(0.0, 0.0, 0.0);
        // Расчёт координат частицы
        for (int nParticle = 0; nParticle < system->numParticles; nParticle++)
        {
//                Vector f_therm = thermostat->getForce(system->particles[nParticle]->state->v);
            f_therm = Vector(0.0, 0.0, 0.0);
            f_LenJon = Vector(0.0, 0.0, 0.0);
            f_dipole = Vector(0.0, 0.0, 0.0);
            for (int j = 0; j < system->numParticles; j++)
            {
                if (j != nParticle)
                {
                    //      ERROR!!!!!
                    f_LenJon += - potential->getGradPotential(Vector(rad_vectors(nParticle, 0), rad_vectors(nParticle, 1), rad_vectors(nParticle, 2)),
                                                              Vector(rad_vectors(j, 0), rad_vectors(j, 1), rad_vectors(j, 2))) * prop->koef_LenJon;
//                        if ((system->particles[nParticle]->state->r - system->particles[j]->state->r).norm() > 2 * prop->radius)
//                        {
//                            f_dipole += interaction->getElectricForce(nParticle);
//                        }
                }
            }
            f = - prop->prandtl * (system->particles[nParticle]->getVelocity() + (f_dipole / prop->prandtlEl) - pow(2.0, 0.5) * prop->knudsen * thermostat->getForce()) * prop->mass;
            //Vector f = f_therm + f_LenJon + f_dipole;
            r_[nParticle] = system->particles[nParticle]->getCoordinate();
            system->particles[nParticle]->state->r = numEq->getCoordinates(system->particles[nParticle]->getCoordinate(), rp_[nParticle], f, prop->mass);

            rad_vectors(nParticle, 0) = system->particles[nParticle]->getCoordinate()(0);
            rad_vectors(nParticle, 1) = system->particles[nParticle]->getCoordinate()(1);
            rad_vectors(nParticle, 2) = system->particles[nParticle]->getCoordinate()(2);
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
                    f_LenJon += - potential->getGradPotential(Vector(rad_vectors(nParticle, 0), rad_vectors(nParticle, 1), rad_vectors(nParticle, 2)),
                                                              Vector(rad_vectors(j, 0), rad_vectors(j, 1), rad_vectors(j, 2))) * prop->koef_LenJon;
                    //                        if ((system->particles[nParticle]->state->r - system->particles[j]->state->r).norm() > 2 * prop->radius)
                    //                        {
                    //                            f_dipole_for_numeq += interaction->getElectricForce(nParticle);
                    //                        }
                }
            }
            f_ = - prop->prandtl * (system->particles[nParticle]->getVelocity() + f_dipole / prop->prandtlEl - pow(2.0, 0.5) * prop->knudsen * thermostat->getForce()) * prop->mass;
            system->particles[nParticle]->state->v = numEq->getVelocity(system->particles[nParticle]->getVelocity(), f_, f, prop->mass);
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

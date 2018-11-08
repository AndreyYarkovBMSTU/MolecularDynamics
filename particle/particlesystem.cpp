#include "particlesystem.h"

void ParticleSystem::setParticle(Particle* _particle)
{
    particles.push_back(_particle);
    numParticles++;
}

void ParticleSystem::setProperties()
{
    v_thermal = pow(2 * prop->k * prop->temperature / particles[0]->mass, 0.5);
    reinolds = 2 * particles[0]->radius * v_thermal / environment->material->kinviscosity;
    ksi = phys::getClausiusMossotti(particles[0]->material->epsilon, environment->material->epsilon);
    friction = 6.0 * M_PI * particles[0]->radius * environment->material->viscosity;
    diffusion = prop->k * prop->temperature / friction;
    dipolemoment0 = prop->koef_dipole * environment->externalfield->electricfield * ksi * pow(2 * particles[0]->radius, 3) / 8.0;
    dipolemoment0_ = environment->externalfield->electricfield * ksi;
}

/*!
 * Получение частицы
 * \return Частицу
 */
Particle* ParticleSystem::getParticle(int _nParticle)
{
    return particles[_nParticle];
}

Vector ParticleSystem::getElectricField(int _nParticle)
{
    extra_electricfield = Vector(0.0, 0.0, 0.0);

    for (int b = 0; b < numParticles; b++)
    {
        extra_electricfield += particles[b]->getElectricField(particles[_nParticle]->getCoordinate());
    }

    return extra_electricfield;
}

Vector ParticleSystem::getForce(int _nParticle)
{
    f.resize(3);
    f[0] = 0.0;
    f[1] = 0.0;
    f[2] = 0.0;
    for (int nParticle = 0; nParticle < numParticles; nParticle++)
    {
        if (nParticle != _nParticle)
        {
            r_ = particles[nParticle]->getCoordinate() - particles[_nParticle]->getCoordinate();
            for (int i = 0; i < 3; i++)
            {
                for (int k = 0; k < 3; k++)
                {
                    kronec_ik = mathematics::getKronec(i, k);
                    for (int j = 0; j < 3; j++)
                    {
                        kronec_ij = mathematics::getKronec(i, j);
                        kronec_jk = mathematics::getKronec(j, k);

                        f[j] += 3 * particles[_nParticle]->dipolemoment(i) * (particles[_nParticle]->dipolemoment(k) * (r_(j) * kronec_ik + r_(k) * kronec_ij - (5 * r_(i) * r_(j) * r_(k) / pow(r_.norm(), 2)))
                                                                   + particles[_nParticle]->dipolemoment(j) * r_(i)) / pow(r_.norm(), 5);
                    }
                }
            }
        }
    }
    return Vector(f[0], f[1], f[2]);
}

Energy ParticleSystem::getSelfEnergy(int _nParticle)
{
    setProperties();

//    return (particles[_nParticle]->dipolemoment.dot(particles[_nParticle]->dipolemoment) - dipolemoment0.dot(dipolemoment0)) / (2 * pow(particles[0]->radius, 3));
    return ((particles[_nParticle]->dipolemoment.dot(particles[_nParticle]->dipolemoment)) - dipolemoment0_.dot(dipolemoment0_)) / 2.0;
}

Energy ParticleSystem::getInteractionEnergy(int _nParticle)
{
    setProperties();

    return - (particles[_nParticle]->dipolemoment).dot(getElectricField(_nParticle)) / 2.0;
}

Energy ParticleSystem::getInductionEnergy(int _nParticle)
{
    setProperties();

    return - (environment->externalfield->electricfield.dot(particles[_nParticle]->dipolemoment - dipolemoment0_)) / 2.0;
}
//Energy ParticleSystem::getInductionEnergy(int _nParticle)
//{
//    dipolemoment0 = environment->externalfield->electricfield * ksi * pow(2 * particles[0]->radius, 3) / 8.0;

//    return - (environment->externalfield->electricfield.dot(particles[_nParticle]->dipolemoment) - environment->externalfield->electricfield.dot(dipolemoment0)) / 2.0;
//}

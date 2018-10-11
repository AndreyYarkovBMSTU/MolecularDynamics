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
    dipolemoment0 = environment->externalfield->electricfield * ksi * pow(2 * particles[0]->radius, 3) / 8.0;
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

Energy ParticleSystem::getSelfEnergy(int _nParticle)
{
    return (particles[_nParticle]->dipolemoment.dot(particles[_nParticle]->dipolemoment) - dipolemoment0.dot(dipolemoment0)) / (2 * pow(particles[0]->radius, 3));
}

Energy ParticleSystem::getInteractionEnergy(int _nParticle)
{
    return - particles[_nParticle]->dipolemoment.dot(getElectricField(_nParticle)) / 2.0;
}

Energy ParticleSystem::getInductionEnergy(int _nParticle)
{
    dipolemoment0 = environment->externalfield->electricfield * ksi * pow(2 * particles[0]->radius, 3) / 8.0;

    return - (environment->externalfield->electricfield.dot(particles[_nParticle]->dipolemoment) - environment->externalfield->electricfield.dot(dipolemoment0)) / 2.0;
}

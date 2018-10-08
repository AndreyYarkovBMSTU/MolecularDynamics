#include "particlesystem.h"

void ParticleSystem::setParticle(Particle* _particle)
{
    particles.push_back(_particle);
    numParticles++;
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
    return (particles[_nParticle]->dipolemoment.dot(particles[_nParticle]->dipolemoment) - dipolemoment0.dot(dipolemoment0)) / (2 * pow(prop->radius, 3));
}

Energy ParticleSystem::getInteractionEnergy(int _nParticle)
{
    return - particles[_nParticle]->dipolemoment.dot(getElectricField(_nParticle)) / 2.0;
}

Energy ParticleSystem::getInductionEnergy(int _nParticle)
{
    dipolemoment0 = environment->externalfield->electricfield * prop->ksi * pow(2 * prop->radius, 3) / 8.0;

    return - (environment->externalfield->electricfield.dot(particles[_nParticle]->dipolemoment) - environment->externalfield->electricfield.dot(dipolemoment0)) / 2.0;
}

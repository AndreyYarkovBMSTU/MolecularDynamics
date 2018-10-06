#ifndef PARTICLESYSTEM_H
#define PARTICLESYSTEM_H

#include "header.h"
#include "particle.h"
#include "environment/environment.h"

struct Energy
{
    double energy;
    Energy(double _energy) :
        energy(_energy)
    {

    }
};

/*!
 * \brief Класс систем частиц

 * С помощью этого класса созданные частицы объединяются в систему.
 */
struct ParticleSystem
{
    int numParticles;                       ///< Число частиц
    Vector dipolemoment0;                    ///< Дипольный момент уединённой частицы
    std::vector<Particle*> particles;       ///< Контейнер частиц
    Environment* environment;               ///< Среда
    Properties* prop;
    /*!
     * Конструктор класса ParticleSystem
     */
    ParticleSystem(Environment* _environment, Properties* _prop) :
        environment(_environment), prop(_prop)
    {
        numParticles = 0;
        dipolemoment0 = environment->externalfield->electricfield * phys::getClausiusMossotti(prop->particlematerial->epsilon, environment->material->epsilon) * pow(2 * prop->radius, 3) / 8.0;
        prop->computePrandtlEl(dipolemoment0);
    }

    /*!
     * Помещение частицы в контейнер
     */
    void setParticle(Particle* _particle)
    {
        particles.push_back(_particle);
        numParticles++;
    }

    /*!
     * Получение частицы
     * \return Частицу
     */
    Particle* getParticle(int _nParticle)
    {
        return particles[_nParticle];
    }

    Vector getElectricField(int _nParticle)
    {
        Vector extra_electricfield = Vector(0.0, 0.0, 0.0);

        for (int b = 0; b < numParticles; b++)
        {
            extra_electricfield += particles[b]->getElectricField(particles[_nParticle]->state->r);
        }

        return extra_electricfield;
    }

    Energy getSelfEnergy(int _nParticle)
    {
        return (particles[_nParticle]->dipolemoment.dot(particles[_nParticle]->dipolemoment) - dipolemoment0.dot(dipolemoment0)) / (2 * pow(prop->radius, 3));
    }

    Energy getInteractionEnergy(int _nParticle)
    {
        return - particles[_nParticle]->dipolemoment.dot(getElectricField(_nParticle)) / 2.0;
    }

    Energy getInductionEnergy(int _nParticle)
    {
        dipolemoment0 = environment->externalfield->electricfield * phys::getClausiusMossotti(prop->particlematerial->epsilon, environment->material->epsilon) * pow(2 * prop->radius, 3) / 8.0;

        return - (environment->externalfield->electricfield.dot(particles[_nParticle]->dipolemoment) - environment->externalfield->electricfield.dot(dipolemoment0)) / 2.0;
    }
};


#endif // PARTICLESYSTEM_H

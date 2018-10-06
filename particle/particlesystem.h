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
    void setParticle(Particle* _particle);

    /*!
     * Получение частицы
     * \return Частицу
     */
    Particle* getParticle(int _nParticle);

    Vector getElectricField(int _nParticle);

    Energy getSelfEnergy(int _nParticle);

    Energy getInteractionEnergy(int _nParticle);

    Energy getInductionEnergy(int _nParticle);
};


#endif // PARTICLESYSTEM_H

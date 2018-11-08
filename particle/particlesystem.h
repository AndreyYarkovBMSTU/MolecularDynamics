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
private:
    std::vector<double> f;
    double kronec_ij;
    double kronec_ik;
    double kronec_jk;
    Vector r_;
public:
    int numParticles;                       ///< Число частиц
    double reinolds;
    double v_thermal;                       ///< Тепловая скорость
    double ksi;                              ///< Фактор Клаузиуса-Моссотти
    double friction;                         ///< Коэффициент внутреннего трения
    double diffusion;
    Vector dipolemoment0;                    ///< Дипольный момент уединённой частицы
    Vector dipolemoment0_;                  ///< Обезразмеренный
    Vector extra_electricfield;
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
    }

    /*!
     * Помещение частицы в контейнер
     */
    void setParticle(Particle* _particle);

    void setProperties();

    /*!
     * Получение частицы
     * \return Частицу
     */
    Particle* getParticle(int _nParticle);

    Vector getElectricField(int _nParticle);

    Vector getForce(int _nParticle);

    Energy getSelfEnergy(int _nParticle);
    Energy getInteractionEnergy(int _nParticle);
    Energy getInductionEnergy(int _nParticle);
};


#endif // PARTICLESYSTEM_H

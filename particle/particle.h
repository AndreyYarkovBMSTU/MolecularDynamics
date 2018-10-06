#ifndef PARTICLE_H
#define PARTICLE_H

#include "header.h"
#include "state.h"
#include "material/material.h"

/*!
 * \brief Класс частиц

 * С помощью этого класса создаётся частица.
 */
struct Particle
{
    Vector dipolemoment;        ///< Дипольный момент
    State* state;               ///< Состояние
    std::string name;           ///< Название
    Material* material;         ///< Материал

    /*!
     * Конструктор класса Particle
     */
    Particle(State* _state, Material* _material) :
        state(_state), material(_material)
    {
        name = "particle";
    }

    Vector getCoordinate();

    Vector getVelocity();

    virtual Vector getElectricField(Vector _r) = 0;
};

/*!
 * \brief Класс диполоидов
 */
struct Dipoloid : Particle
{
    /*!
     * Конструктор класса Dipoloid
     */
    Dipoloid(State* _state, Material* _material) :
        Particle(_state, _material)
    {
        name = "dipoloid";
    }

    Vector getElectricField(Vector _r);
};

#endif // PARTICLE_H

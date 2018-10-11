#ifndef PARTICLE_H
#define PARTICLE_H

#include "header.h"
#include "state.h"
#include "material/material.h"
#include "object.h"

/*!
 * \brief Класс частиц

 * С помощью этого класса создаётся частица.
 */
struct Particle
{
    double radius;
    double mass;                ///< Масса частицы
    Vector dipolemoment;        ///< Дипольный момент
    State* state;               ///< Состояние
    std::string name;           ///< Название
    Material* material;         ///< Материал

    /*!
     * Конструктор класса Particle
     */
    Particle(State* _state, Material* _material, Object* _obj) :
        state(_state), material(_material)
    {
        name = "particle";
        radius = _obj->radius;
        mass = material->density * 4/3 * M_PI * pow(radius, 3) * phys::kg();
    }

    Vector getCoordinate();
    Vector getVelocity();

    void setCoordinate(Vector _r);
    void setVelocity(Vector _v);

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
    Dipoloid(State* _state, Material* _material, Object* _obj) :
        Particle(_state, _material, _obj)
    {
        name = "dipoloid";
    }

    Vector getElectricField(Vector _r);
};

#endif // PARTICLE_H

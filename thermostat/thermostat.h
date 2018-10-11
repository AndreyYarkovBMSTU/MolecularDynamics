#ifndef THERMOSTAT_H
#define THERMOSTAT_H

#include "header.h"
#include "environment/environment.h"

/*!
 * \brief Класс термостатов

 * С помощью этого класса создаётся термостат, в который помещают систему (термостат включает в себя среду, в которой находится система).
 */
struct Thermostat
{
    Properties* prop;               ///< Свойства
    Environment* environment;       ///< Среда

    /*!
     * Конструктор класса Thermostat
     */
    Thermostat(Environment* _environment, Properties* _prop) :
        environment(_environment), prop(_prop)
    {

    }

    virtual Vector getForce() = 0;          ///< Расчёт силы, действующей со стороны среды
};

/*!
 * \brief Класс Броуновского термостата (описываемого уравнением Ланжевена)
 */
struct Langevin : Thermostat
{
    std::random_device rd;
    /*!
     * Конструктор класса Langevin
     */
    Langevin(Environment* _medium, Properties* _prop) :
        Thermostat(_medium, _prop)
    {

    }

    /*!
     * Рассчитывет силу, действующую на частицу, обусловленную броуновским движением и сопротивлением среды
     * \return Результирующую силу как сумму случайной броуновской силы и силы сопротивления среды
     */
    Vector getForce();
};

#endif // THERMOSTAT_H

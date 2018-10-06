#ifndef NUMERALEQUATIONS_H
#define NUMERALEQUATIONS_H

#include "header.h"

/*!
 * \brief Класс численных методов

 * С помощью этого класса рассчитывается состояние частицы в определённый момент времени.
 */
class NumeralEquations
{
public:
    double dt;          ///< Шаг по времени

    /*!
     * Конструктор класса NumeralEquations
     */
    NumeralEquations(double _dt) :
        dt(_dt)
    {

    }

    virtual Vector getCoordinates(Vector _r, Vector _r_previous, Vector _force, double _mass) = 0;      ///< Получение радиус-вектора частицы
    virtual Vector getVelocity(Vector _v, Vector _force_, Vector _force, double _mass) = 0;         ///< Получение вектора скорости частицы
};

/*!
 * \brief Класс численной схемы Верле
 */
class Verle : public NumeralEquations
{
public:

    /*!
     * Конструктор класса Verle
     */
    Verle(double _dt) :
        NumeralEquations(_dt)
    {

    }

    /*!
     * Получение радиус-вектора частицы
     * \return Радиус-вектор частицы, рассчитываемый по формуле: \n
       \f[
            r(t + dt) = 2 r(t) - r(t - dt) + \frac{F}{m} \cdot (dt)^{2}
       \f]
     */
    Vector getCoordinates(Vector _r, Vector _r_previous, Vector _force, double _mass);

    /*!
     * Получение вектора скорости частицы
     * \return Вектор скорости частицы, рассчитываемый по формуле: \n
       \f[
            v(t + dt) = \frac{r(t + dt)}{2 dt}
       \f]
     */
    Vector getVelocity(Vector _v, Vector _force_, Vector _force, double _mass);
};

#endif // NUMERALEQUATIONS_H

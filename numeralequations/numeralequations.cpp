#include "numeralequations.h"

Vector Verle::getCoordinates(Vector _r, Vector _r_previous, Vector _force, double _mass)
{
    return 2 * _r - _r_previous + _force * dt * dt / _mass;
}

/*!
 * Получение вектора скорости частицы
 * \return Вектор скорости частицы, рассчитываемый по формуле: \n
   \f[
        v(t + dt) = \frac{r(t + dt)}{2 dt}
   \f]
 */
Vector Verle::getVelocity(Vector _v, Vector _force_, Vector _force, double _mass)
{
    return _v + (_force_ + _force) * dt / (2 * _mass);
}

#include "numeralequations.h"

Vector Verle::getCoordinates(Vector _r, Vector _r_previous, Vector _a)
{
    return 2 * _r - _r_previous + _a * dt * dt;
}

/*!
 * Получение вектора скорости частицы
 * \return Вектор скорости частицы, рассчитываемый по формуле: \n
   \f[
        v(t + dt) = \frac{r(t + dt)}{2 dt}
   \f]
 */
Vector Verle::getVelocity(Vector _v, Vector _a_, Vector _a)
{
    return _v + (_a + _a_) * dt / 2;
}

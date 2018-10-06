#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "header.h"

/*!
 * \brief Класс потенциалов

 * С помощью этого класса задаётся потенциал взаимодействия между частицами и рассчитывается его градиент.
 */
class Potential
{
public:
    double scale;       ///< Величина размерности

    /*!
     * Конструктор класса Potential
     */

    Potential(double _scale)
    {
        scale = 1.0 / _scale; // т.к. это величина, а не единица размерности
    }

    virtual double getPotential(Vector _r_1, Vector _r_2) = 0;              ///< Получение потенциала
    virtual Vector getGradPotential(Vector _r_1, Vector _r_2) = 0;          ///< Получение градинта потенциала
};

/*!
 * \brief Класс IPL потенциалов
 */
class IPL : public Potential
{
public:

    int n;              ///< Порядок IPL потенциала

    /*!
     * Конструктор класса IPL
     */
    IPL(int _n, double _scale) :
        Potential(_scale)
    {
        n = _n;
    }
    ~IPL()
    {

    }

    /*!
     * Получение IPL потенциала взаимодействия
     * \return Потенциал взаимодействия, рассчитываемый по формуле: \n
       \f[
            \phi = \biggl(\frac{2}{|(r1 - r2)| \cdot scale}\biggr)^{n}
       \f]
     */
    double getPotential(Vector _r_1, Vector _r_2)
    {
        return pow(2.0 / ((_r_1 - _r_2).norm() * scale), n);
    }

    /*!
     * Получение градиента IPL потенциала взаимодействия
     * \return Градиент потенциала взаимодействия, рассчитываемый по формуле: \n
       \f[
            \nabla\phi = - n \cdot \phi \cdot \biggl(\frac{2}{|(r1 - r2)| \cdot scale}\biggr)^{2} \cdot (r1 - r2) \cdot scale
       \f]
     */
    Vector getGradPotential(Vector _r_1, Vector _r_2)
    {
        return - n * getPotential(_r_1, _r_2) * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 2) * (_r_1 - _r_2) * scale;
    }
};

/*!
 * \brief Класс потенциала Леннарда-Джонса
 */
class LJ : public Potential
{
public:

    /*!
     * Конструктор класса LJ
     */
    LJ(double _scale) :
        Potential(_scale)
    {

    }
    ~LJ()
    {

    }

    /*!
     * Получение потенциала взаимодействия Леннарда-Джонса
     * \return Потенциал взаимодействия, рассчитываемый по формуле: \n
       \f[
            \phi = \biggl(\frac{2}{|(r1 - r2)| \cdot scale}\biggr)^{12} - \biggl(\frac{2}{|(r1 - r2)| \cdot scale}\biggr)^{6}
       \f]
     */
    double getPotential(Vector _r_1, Vector _r_2)
    {
        if ((_r_1 - _r_2).norm() * scale < 2)
        {
            return pow(2.0 / ((_r_1 - _r_2).norm() * scale), 12);
        }
        else
        {
            return - pow(2.0 / ((_r_1 - _r_2).norm() * scale), 6);
            //return 0;
        }
    }

    /*!
     * Получение градиента потенциала взаимодействия Леннарда-Джонса
     * \return Градиент потенциала взаимодействия, рассчитываемый по формуле: \n
       \f[
            \nabla\phi = - 12 \cdot \phi \cdot \biggl(\frac{2}{|(r1 - r2)| \cdot scale}\biggr)^{2} \cdot (r1 - r2) \cdot scale + 6 \cdot \phi \cdot \biggl(\frac{2}{|(r1 - r2)| \cdot scale}\biggr)^{2} \cdot (r1 - r2) \cdot scale
       \f]
     */
    Vector getGradPotential(Vector _r_1, Vector _r_2)
    {
        if ((_r_1 - _r_2).norm() * scale < 2)
        {
            return - 12 * getPotential(_r_1, _r_2) * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 2) * (_r_1 - _r_2) * scale;
        }
        else
        {
            return - 6 * getPotential(_r_1, _r_2) * pow(2.0 / ((_r_1 - _r_2).norm() * scale), 2) * (_r_1 - _r_2) * scale;
        }
    }
};

#endif // POTENTIAL_H

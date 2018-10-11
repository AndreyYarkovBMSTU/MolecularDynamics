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
    virtual double getPotential_dipole(Vector _r_1, Vector _r_2, Vector _dipolemoment_1, Vector _dipolemoment_2) = 0;
    virtual Vector getGradPotential_dipole(Vector _r_1, Vector _r_2, Vector _dipolemoment_1, Vector _dipolemoment_2) = 0;
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
    double getPotential(Vector _r_1, Vector _r_2);

    /*!
     * Получение градиента IPL потенциала взаимодействия
     * \return Градиент потенциала взаимодействия, рассчитываемый по формуле: \n
       \f[
            \nabla\phi = - n \cdot \phi \cdot \biggl(\frac{2}{|(r1 - r2)| \cdot scale}\biggr)^{2} \cdot (r1 - r2) \cdot scale
       \f]
     */
    Vector getGradPotential(Vector _r_1, Vector _r_2);

    double getPotential_dipole(Vector _r_1, Vector _r_2, Vector _dipolemoment_1, Vector _dipolemoment_2);
    Vector getGradPotential_dipole(Vector _r_1, Vector _r_2, Vector _dipolemoment_1, Vector _dipolemoment_2);
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
    double getPotential(Vector _r_1, Vector _r_2);

    /*!
     * Получение градиента потенциала взаимодействия Леннарда-Джонса
     * \return Градиент потенциала взаимодействия, рассчитываемый по формуле: \n
       \f[
            \nabla\phi = - 12 \cdot \phi \cdot \biggl(\frac{2}{|(r1 - r2)| \cdot scale}\biggr)^{2} \cdot (r1 - r2) \cdot scale + 6 \cdot \phi \cdot \biggl(\frac{2}{|(r1 - r2)| \cdot scale}\biggr)^{2} \cdot (r1 - r2) \cdot scale
       \f]
     */
    Vector getGradPotential(Vector _r_1, Vector _r_2);

    double getPotential_dipole(Vector _r_1, Vector _r_2, Vector _dipolemoment_1, Vector _dipolemoment_2);
    Vector getGradPotential_dipole(Vector _r_1, Vector _r_2, Vector _dipolemoment_1, Vector _dipolemoment_2);
};

class StockMayer : public Potential
{
public:

    StockMayer(double _scale) :
        Potential(_scale)
    {

    }
    ~StockMayer()
    {

    }

    double getPotential(Vector _r_1, Vector _r_2);
    Vector getGradPotential(Vector _r_1, Vector _r_2);

    double getPotential_dipole(Vector _r_1, Vector _r_2, Vector _dipolemoment_1, Vector _dipolemoment_2);
    Vector getGradPotential_dipole(Vector _r_1, Vector _r_2, Vector _dipolemoment_1, Vector _dipolemoment_2);
};

#endif // POTENTIAL_H

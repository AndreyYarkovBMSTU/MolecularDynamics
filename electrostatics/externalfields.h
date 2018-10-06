#ifndef EXTERNALFIELDS_H
#define EXTERNALFIELDS_H

#include "header.h"
#include "properties.h"

/*!
 * \brief Класс внешнего электрического поля

 * С помощью этого класса задаётся внешнее электрическое поле.
 */
class ExternalFields
{
public:
    Vector electricfield;           ///< Напряжённость электрического поля
    double omega;                   ///< Угловая скорость вращения поля
    Properties* prop;

    /*!
     * Конструктор класса ExternalFields
     */
    ExternalFields(Properties* _prop, Vector _E, double _omega = 0.0) :
        prop(_prop), electricfield(_E), omega(_omega)
    {

    }
    virtual Vector getElectricField(int _numAngle = 0) = 0;
};

/*!
 * \brief Класс направленного электрического поля
 */
class DirectedField : public ExternalFields
{
public:   
    /*!
     * Конструктор класса DirectedField
     */
    DirectedField(Properties* _prop, Vector _E) :
        ExternalFields(_prop, _E)
    {

    }

    Vector getElectricField(int _numAngle);
};

class RotatingField : public ExternalFields
{
public:
    RotatingField(Properties* _prop, Vector _E, double _omega) :
        ExternalFields(_prop, _E, _omega)
    {

    }

    Vector getElectricField(int _numAngle);
};

#endif // EXTERNALFIELDS_H

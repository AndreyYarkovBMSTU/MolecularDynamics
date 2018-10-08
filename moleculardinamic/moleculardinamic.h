#ifndef MOLECULARDINAMIC_H
#define MOLECULARDINAMIC_H

#include "header.h"
#include "particle/particlesystem.h"
#include "numeralequations/numeralequations.h"
#include "thermostat/thermostat.h"
#include "moleculardinamic/potential.h"
#include "electrostatics/externalfields.h"
#include "electrostatics/interaction.h"
#include "electrostatics/methods.h"

/*!
 * \brief Класс метода молекулярной динамики

 * С помощью этого класса рассчитывается конфигурация системы в последующие моменты времени, а также, при необходимости, происходит запись в файл.
 */
class MolecularDinamic
{
private:
    int k;
    double rad;
    Vector f;
    Vector f_;
    Vector f_therm;
    Vector f_LenJon;
    Vector f_dipole;
    Matrix R_;
    vector<Vector> r_;            // Радиус-вектор частицы в настоящий момент времени
    vector<Vector> rp_;           // Радиус-вектор частицы в предыдущий момент времени
    std::string path;
    std::ofstream out;
    const char * c;
public:
    ParticleSystem* system;
    Thermostat* thermostat;             ///< Термостат
    NumeralEquations* numEq;            ///< Численная схема
    Potential* potential;               ///< Твёрдая стенка
    Methods* method;                    ///< Метод расчёта взаимодействия
    Interaction* interaction;
    Properties* prop;                   ///< Свойства

    /*!
     * Конструктор класса MolecularDinamic
     */
    MolecularDinamic(ParticleSystem* _system, Methods* _method, Properties* _prop, std::string _nameThermostat, std::string _nameNumEq, std::string _namePotential) :
        system(_system), method(_method), prop(_prop)
    {        
        interaction = new Interaction(_system, _method, _prop);

        if (_nameThermostat == "langevin")
        {
            thermostat = new Langevin(system->environment, _prop);
        }
        if (_nameNumEq == "verle")
        {
            numEq = new Verle(_prop->timestep);
        }
        if (_namePotential == "LJ")
        {
            potential = new LJ(prop->radius);
        }
    }

    /*!
     * Выводит состояния всех частиц
     */
    void dump();

    /*!
     * Расситывает и выводит на экран состояние определённой частицы на определённом кадре
     */
    void computer(int _nParticle, int _nFrame);

    /*!
     * Расcчитывает и записывает в файл конфигурацию системы в последующие моменты времени
     */
    void record();
};

#endif // MOLECULARDINAMIC_H

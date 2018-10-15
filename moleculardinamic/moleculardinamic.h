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
    double tau;                         ///< Коэффициент обезразмеривания по времени
    double koef_demping;
    double koef_LenJon;
    double koef_randForce;
    double koef_dipole;

    Vector v;
    Vector f;
    Vector f_;
    Vector a;
    Vector a_;
    Vector f_LenJon;
    Vector f_dipole;
    Matrix R_;
    vector<Vector> r_;            // Радиус-вектор частицы в настоящий момент времени
    vector<Vector> rp_;           // Радиус-вектор частицы в предыдущий момент времени
    Matrix _R;
    std::string path;
    std::ofstream out;
    std::string outputline;
    const char * c;
public:
    double t0;                          ///< Обезразмеренное время
    ParticleSystem* system;
    Thermostat* thermostat;             ///< Термостат
    NumeralEquations* numEq;            ///< Численная схема
    Potential* potential;               ///< Твёрдая стенка
    Methods* method;                    ///< Метод расчёта взаимодействия
    Interaction* interaction;
    Properties* prop;                   ///< Свойства
    std::string nameThermostat;

    /*!
     * Конструктор класса MolecularDinamic
     */
    MolecularDinamic(ParticleSystem* _system, Methods* _method, Properties* _prop, std::string _nameThermostat, std::string _nameNumEq, std::string _namePotential) :
        system(_system), method(_method), prop(_prop)
    {
        tau = 2 * system->particles[0]->radius / system->v_thermal;
        t0 = prop->timestep / tau;

        interaction = new Interaction(system, method, prop);   

        _R.resize(system->numParticles, 3);
        interaction = new Interaction(_system, _method, _prop);

        if (_nameThermostat == "langevin")
        {
            nameThermostat = _nameThermostat;
            thermostat = new Langevin(system->environment, prop);
        }
        else if (_nameThermostat == "brownian")
        {
            nameThermostat = _nameThermostat;
            thermostat = new Brownian(system->environment, prop);
        }

        if (_nameNumEq == "verle")
        {
            numEq = new Verle(t0);
        }
        if (_namePotential == "LJ")
        {
            potential = new LJ(1.0);
        }
    }

    /*!
     * Выводит состояния всех частиц
     */
    void dump();

    void computeKoef();
    /*!
     * Расситывает и выводит на экран состояние определённой частицы на определённом кадре
     */
    void computer(int _nParticle, int _nFrame);

    /*!
     * Расcчитывает и записывает в файл конфигурацию системы в последующие моменты времени
     */
    void record();

    void recordmove(int _numfile, int _nFrame);

    void setVelocityLangevin();
    void setCoordinateLangevin();
    void setVelocityBrownian();
    void setCoordinateBrownian();
};

#endif // MOLECULARDINAMIC_H

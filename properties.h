#ifndef PROPERTIES_H
#define PROPERTIES_H

#include "header.h"
#include "material/material.h"

class Properties
{
public:
    int numFrames = 10000;                   ///< Число кадров
    int numAngles = 40;                     ///< Число шагов при усреднении при вращении
    double timestep = 1e-2 * phys::us();           ///< Временной шаг
    double temperature = 293.15 * phys::K();      ///< Температура
    double k = 1.38e-23 * phys::J() / phys::K();        ///< Постоянная Больцмана
    double radius = 1.0 * phys::um();             ///< Радиус частицы
    double epsilon0 = 8.85 * 1e-12;             ///< Электрическая постоянная
    double ksi;                              ///< Фактор Клаузиуса-Моссотти
    double mass;                            ///< Масса частицы
    double prandtl;                         ///< Число Прандтля
    double prandtlEl;                        ///< Электрическое число Прандтля
    double schmidt;                         ///< Число Шмидта
    double knudsen;                         ///< Число Кнудсена
    double koef_LenJon;
    Material* particlematerial;             ///< Материал частиц
    Material* solvent;                      ///< Материал среды
    std::string path = "output/time/";       ///< Путь к файлу
    std::string filetype = ".txt";          ///< Тип файла

    /*!
     * Конструктор класса Properties
     * Расчёт массы частицы по формуле: \n
     \f[
        m = \rho \cdot \frac{4}{3} \pi r^{3}
     \f]
     */
    Properties(Material* _particlematerial, Material* _solvent) :
        particlematerial(_particlematerial), solvent(_solvent)
    {
        mass = _particlematerial->density * 4/3 * M_PI * pow(radius, 3) * phys::kg();

        ksi = phys::getClausiusMossotti(particlematerial->epsilon, solvent->epsilon);

        double friction = 6.0 * M_PI * radius * _solvent->viscosity;
        double diffusion = k * temperature / friction;

        schmidt = (friction / _solvent->density) / diffusion;
        double timedif = 1e-6;
        prandtl = 18 * pow((pow(diffusion * timedif, 0.5) / (2 * radius)), 2) * schmidt;
        knudsen = pow(timedif * diffusion, 0.5) / (2 * radius);

        timestep = timestep / timedif;
    }

    void computePrandtlEl(Eigen::Vector3d _dipolemoment0)
    {
        double viscosityEl = (1.0 / 3 * M_PI * 2 * radius) * (_dipolemoment0.dot(_dipolemoment0) / (pow(2 * radius, 2))) / (4 * M_PI * solvent->epsilon * epsilon0 * pow(2 * radius, 2));
        prandtlEl = solvent->viscosity / viscosityEl;
    }
};

#endif // PROPERTIES_H

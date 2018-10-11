#ifndef PROPERTIES_H
#define PROPERTIES_H

#include "header.h"
#include "material/material.h"

class Properties
{
public:
    int numFrames = 100000;                   ///< Число кадров
    int numAngles = 40;                     ///< Число шагов при усреднении при вращении
    double radius = 1e-6;
    double timestep = 1e-8;           ///< Временной шаг
    double temperature = 293.15 * phys::K();      ///< Температура
    double k = 1.38e-23 * phys::J() / phys::K();        ///< Постоянная Больцмана
    double epsilon0 = 8.85 * 1e-12;             ///< Электрическая постоянная
    double koef_LenJon = 1e-14; //1e-14
    double koef_Brown = 1e12;
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

    }
};

#endif // PROPERTIES_H

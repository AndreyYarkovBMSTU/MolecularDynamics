#ifndef PROPERTIES_H
#define PROPERTIES_H

#include "header.h"
#include "object.h"
#include "material/material.h"

class Properties
{
public:
    int numFrames = 500000;                   ///< Число кадров
    int numAngles = 40;                     ///< Число шагов при усреднении при вращении
    double radius = 1e-6;
    double timestep = 1e-4;           ///< Временной шаг
    double temperature = 293.15 * phys::K();      ///< Температура
    double k = 1.38e-23 * phys::J() / phys::K();        ///< Постоянная Больцмана
    double epsilon0 = 8.85 * 1e-12;             ///< Электрическая постоянная
    double koef_omega = 2 * 1e-4;//1.0e-3;//1e2;
    double koef_dipole = 75.0;//75;      // используется в particlesystem
    double koef_LenJon = 1e-19;
    double koef_Brown = 15;//15;
    Object* obj0;                           ///< стандартный объект
    Material* particlematerial;             ///< Материал частиц
    Material* solvent;                      ///< Материал среды
    std::string path = "output/";       ///< Путь к файлу
    std::string filetype = ".txt";          ///< Тип файла

    Properties(Material* _particlematerial, Material* _solvent) :
        particlematerial(_particlematerial), solvent(_solvent)
    {
        obj0 = new Object(radius);
    }
};

#endif // PROPERTIES_H

#ifndef MATERIAL_H
#define MATERIAL_H

#include "header.h"

/*!
 * \brief Класс материалов

 * С помощью этого класса вещества наделяются свойствами, присущими их материалам.
 */
struct Material {
    double sigma;
    double epsilon;             ///< Диэлектрическая пронецаемость
    double mu;
    double density;             ///< Плотность
    double viscosity;           ///< Вязкость
    std::string name;
    std::string formula;
    std::string electricType;

    /*!
     * Конструктор класса Material
     */
    Material(std::string _name) :
        name(_name)
    {
        if (name == "diamond") {
            formula = "C";
            electricType = "dielectric";
            epsilon = 1.0;
            mu = 1.0;
            viscosity = 8.9e-4;  // [Pa * s]
            density = 0.0;
        }
    }
    ~Material()
    {

    }
};

/*!
 * \brief Класс материала "Диоксид кремния"
 */
struct SiliconDioxide : Material
{
    /*!
     * Конструктор класса SiliconDioxide
     */
    SiliconDioxide(std::string _name) :
        Material(_name)
    {
        name = "silica";
        formula = "SiO2";
        electricType = "dielectric";
        setProperties();
    }

    /*!
     * Установление свойст материала
     */
    void setProperties()
    {
    //! epsilon = 4.0;    // 100%
        epsilon = 2.925;  // 75 %
        //epsilon = 2000.925;  // 75 %
        mu = 1.0;
        viscosity = 8.9e-4;  // [Pa * s]
        density = 2.648 * phys::g() / pow(phys::sm(),3);
    }
};

/*!
 * \brief Класс материала "Вода"
 */
struct Water : Material
{
    /*!
     * Конструктор класса Water
     */
    Water(std::string _name) :
        Material(_name)
    {
        name = "water";
        formula = "H2O";
        electricType = "dielectric";
        setProperties();
    }

    /*!
     * Установление свойст материала
     */
    void setProperties()
    {
        epsilon = 80.1;  // T = 20 C;
        mu = 1.0;
        viscosity = 8.9e-4;  // [Pa * s]
        density = 998; // [kg / m3]
    }
};

#endif // MATERIAL_H

#ifndef STATE_H
#define STATE_H

#include "header.h"

/*!
 * \brief Класс состояний частиц

 * С помощью этого класса задаётся состояние частицы.
 */
class State
{
public:
    int number;             ///< Номер частицы
    Vector r;               ///< Радиус-вектор частицы
    Vector v;               ///< Скорость частицы

    /*!
     * Конструктор класса State
     */
    State(int _number, Vector _r, Vector _v) :
        number(_number), r(_r), v(_v)
    {

    }

    /*!
     * Перегрузка оператора "<<" на вывод состояния частицы
     * \return Состояние частицы: порядковый номер, радиус-вектор, вектор скорости
     */
    friend std::ostream& operator <<(std::ostream& _flux, const State* _state)
    {
        std::ios::fmtflags mode = _flux.flags();
        _flux.setf(std::ios::fixed,std::ios::floatfield);
        auto prec = _flux.precision(9);

        _flux << _state->number << " " << _state->r(0)  << " " << _state->r(1) << " " << _state->r(2) << " " << _state->v(0)  << " " << _state->v(1) << " " << _state->v(2) << "\n";

        _flux.precision(prec);
        _flux.setf(mode,std::ios::floatfield);
        return _flux;
    }
};


#endif // STATE_H

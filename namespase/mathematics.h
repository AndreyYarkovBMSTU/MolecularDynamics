#ifndef MATHEMATICS_H
#define MATHEMATICS_H

#include "Eigen/Eigen"

typedef Eigen::MatrixXd Matrix;

#pragma once

namespace mathematics
{
    double getKronec(int i, int j);

    Matrix getUnitMatrix(int n);
}

#endif // MATHEMATICS_H

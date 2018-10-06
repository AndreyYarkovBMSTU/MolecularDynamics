#include "mathematics.h"

double mathematics::getKronec(int i, int j)
{
    if (i == j)
        return 1;
    else
        return 0;
}

Matrix mathematics::getUnitMatrix(int n)
{
    Matrix I(n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
                I(i, j) = 1;
            else
                I(i, j) = 0;
        }
    }

    return I;
}

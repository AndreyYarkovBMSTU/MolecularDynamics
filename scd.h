#ifndef SCD_H
#define SCD_H

#include "header.h"
#include "externalfields.h"
#include "particle/particlesystem.h"

class SCD
{
public:
    double d;                               ///< Диаметр частицы
    ParticleSystem* system;

    SCD(double _d, ParticleSystem* _system) :
        d(_d), system(_system)
    {

    }
//!  ERROR!!!!!
//! Matrix getBlock(Particle* _partilce1, Particle* _partilce2)
    Matrix getBlock(int _a, int _b, Vector _ra, Vector _rb)
    {
        Matrix B(3, 3);
        Vector n_ab = (_ra - _rb) / (_ra - _rb).norm();
        double kronec_ab = mathematics::getKronec(_a,_b);

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                double kronec_ij = mathematics::getKronec(i, j);

                if (kronec_ab == 1)
                {
                    B(i, j) = 0;
                }
                else
                {
                    B(i, j) = (3 * n_ab(i) * n_ab(j) - kronec_ij) * (1 - kronec_ab) / pow((_ra - _rb).norm(), 3);
                }
            }
        }

        return B;
    }

    Matrix getInteractionsMatrix()
    {
        int a, b;
        Matrix B[system->numParticles * system->numParticles];
        Matrix block_Matrix(3 * system->numParticles, 3 * system->numParticles);

        Matrix I = mathematics::getUnitMatrix(3 * system->numParticles);

        for (int t = 0; t < system->numParticles * system->numParticles; t++)
        {
            a = t / system->numParticles;
            b = t % system->numParticles;
            B[t] = getBlock(a, b, system->particles[a]->state->r, system->particles[b]->state->r);
        }

        for (int i = 0; i < 3 * system->numParticles; i++)
        {
            for (int j  = 0; j < 3 * system->numParticles; j++)
            {
                int k = system->numParticles * (i  / 3) + j / 3;
                block_Matrix(i, j) = B[k](i % 3, j % 3);
            }
        }
//!   ERROR!!!!
        // chi = CM  d^3 / 8;
        return ((8.0 / (phys::getClausiusMossotti(system->particles[0]->material->epsilon, system->environment->material->epsilon) * pow(d, 3)))) * I - block_Matrix;
    }

    void setDipoleMoment()
    {
        Matrix E(3 * system->numParticles, 1);
        Matrix dipolemoment(E.rows(), 1);

        for (int i = 0; i < 3 * system->numParticles; i++)
        {
            E(i, 0) = system->environment->externalfield->electricfield(i % 3);
        }
        dipolemoment = getInteractionsMatrix().ldlt().solve(E);

        int k = 0;
        for (int i = 0; i < system->numParticles; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                system->particles[i]->dipolemoment(j) = dipolemoment(k);
                k++;
            }
        }
    }
};

#endif // SCD_H

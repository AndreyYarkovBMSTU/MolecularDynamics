#include "methods.h"

void NonIteractingDipoles::setDipoleMoment()
{
    dipolemoment = system->ksi * pow(2 * system->particles[0]->radius, 3) * system->environment->externalfield->electricfield / 8.0;
    for (int i = 0; i < system->numParticles; i++)
    {
        system->particles[i]->dipolemoment = dipolemoment;
    }
}

Matrix SelfConsistentDipoles::getBlock(Particle* particle1, Particle* particle2)
{
    Block.resize(3, 3);
    n_ab = (particle1->getCoordinate() - particle2->getCoordinate()) / (particle1->getCoordinate() - particle2->getCoordinate()).norm();
    kronec_ab = mathematics::getKronec(particle1->state->number, particle2->state->number);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            kronec_ij = mathematics::getKronec(i, j);

            if (kronec_ab == 1)
            {
                Block(i, j) = 0;
            }
            else
            {
                Block(i, j) = (3 * n_ab(i) * n_ab(j) - kronec_ij) * (1 - kronec_ab) / pow((particle1->getCoordinate() - particle2->getCoordinate()).norm(), 3);
            }
        }
    }

    return Block;
}

Matrix SelfConsistentDipoles::getInteractionsMatrix()
{
    B.resize(system->numParticles * system->numParticles);
    Block_Matrix.resize(3 * system->numParticles, 3 * system->numParticles);

    for (int t = 0; t < system->numParticles * system->numParticles; t++)
    {
        a = t / system->numParticles;
        b = t % system->numParticles;
        B[t] = getBlock(system->particles[a], system->particles[b]);
    }

    for (int i = 0; i < 3 * system->numParticles; i++)
    {
        for (int j  = 0; j < 3 * system->numParticles; j++)
        {
            k = system->numParticles * (i  / 3) + j / 3;
            Block_Matrix(i, j) = B[k](i % 3, j % 3);
        }
    }

    return ((8.0 / (system->ksi * pow(system->particles[0]->radius * 2, 3)))) * mathematics::getUnitMatrix(3 * system->numParticles) - Block_Matrix;
}

void SelfConsistentDipoles::setDipoleMoment()
{
    E.resize(3 * system->numParticles, 1);
    Dipolemoment.resize(E.rows(), 1);

    for (int i = 0; i < 3 * system->numParticles; i++)
    {
        E(i, 0) = system->environment->externalfield->electricfield(i % 3);
    }
    Dipolemoment = getInteractionsMatrix().ldlt().solve(E);

    k = 0;
    for (int i = 0; i < system->numParticles; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            system->particles[i]->dipolemoment(j) = Dipolemoment(k);
            k++;
        }
    }
}

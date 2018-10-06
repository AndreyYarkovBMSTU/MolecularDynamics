#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include "header.h"
#include "material/material.h"
#include "electrostatics/externalfields.h"

struct Environment
{
    Material* material;
    ExternalFields* externalfield;

    Environment(Material* _material, ExternalFields* _externalfield) :
        material(_material), externalfield(_externalfield)
    {

    }
};

#endif // ENVIRONMENT_H

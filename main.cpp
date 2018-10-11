#include <iostream>
#include "header.h"
#include "electrostatics/externalfields.h"
#include "electrostatics/methods.h"
#include "electrostatics/interaction.h"
#include "environment/environment.h"
#include "material/material.h"
#include "particle/particlesystem.h"
#include "thermostat/thermostat.h"
#include "moleculardinamic/moleculardinamic.h"
#include "properties.h"
#include "object.h"


int main()
{
    // Создание системы частиц.
    Material* particlematerial = new SiliconDioxide("silica");
    Material* solvent = new Water("water");

    Properties* prop = new Properties(particlematerial, solvent);

    Object* obj = new Object(prop->radius);

    State* state1 = new State(1, Vector(-1.1, 0.0, 0.0), Vector(0.0, 0.0, 0.0));
    State* state2 = new State(2, Vector(1.1, 0.0, 0.0), Vector(0.0, 0.0, 0.0));
    State* state3 = new State(3, Vector(2.0, 4.0, 0.0), Vector(0.0, 0.0, 0.0));
    State* state4 = new State(2, Vector(-4.0, 3.0, 0.0), Vector(0.0, 0.0, 0.0));
    State* state5 = new State(3, Vector(-4.0, 1.0, 0.0), Vector(0.0, 0.0, 0.0));

    Particle* particle1 = new Dipoloid(state1, particlematerial, obj);
    Particle* particle2 = new Dipoloid(state2, particlematerial, obj);
    Particle* particle3 = new Dipoloid(state3, particlematerial, obj);
    Particle* particle4 = new Dipoloid(state4, particlematerial, obj);
    Particle* particle5 = new Dipoloid(state5, particlematerial, obj);

    //ExternalFields* externalfield = new DirectedField(prop, Vector(1e3, 0.0, 0.0));
    ExternalFields* externalfield = new RotatingField(prop, Vector(1e3, 0.0, 0.0), 30 * 1e3);
    Environment* environment = new Environment(solvent, externalfield);

    ParticleSystem* system = new ParticleSystem(environment, prop);
    system->setParticle(particle1);
    system->setParticle(particle2);
    system->setParticle(particle3);
    system->setParticle(particle4);
    system->setParticle(particle5);
    system->setProperties();

    // Расчёт

    Methods* method = new SelfConsistentDipoles(system, prop);

    method->setDipoleMoment();

    MolecularDinamic* moleculardinamic = new MolecularDinamic(system,
                                                              method,
                                                              prop,
                                                              "langevin",
                                                              "verle",
                                                              "LJ");

    moleculardinamic->record();

    std::cout << "Reinolds: " << system->reinolds << std::endl;
    std::cout << "t0: " << moleculardinamic->t0 << std::endl;

//    moleculardinamic->computer(2, 100000);

//    for (auto &iParticle : particlesystem->particles)
//    {
//        std::cout << iParticle->state;
//    }


//    std::vector<int> hist;
//    hist.resize(1030);
//    std::string path = "output/gauss.txt";
//    std::ofstream out;
//    out.open(path, std::ios::out | std::ios::binary);
//    for (int iter = 0; iter < 1000000; iter++)
//    {
//        std::random_device rd;
//        std::mt19937 gen(rd());
//        std::normal_distribution<> d(0, 1);
//        double f = d(gen);
//        for (int i = 0; i < hist.size(); i++)
//        {
//            if (f >= -3 + i * 6.0 / hist.size() && f <= -3 + (1 + i) * 6.0 / hist.size())
//            {
//                hist[i] += 1;
//            }
//        }
//    }
//    double delta = 3.0 / hist.size();
//    std::string o = std::to_string(-3 + delta) + "   " + std::to_string(hist[0]) + "\n";
//    for (int i = 1; i < hist.size(); i++)
//    {
//        delta = (2 * i + 1) * 3.0 / hist.size();
//        o += std::to_string(-3 + delta) + "   " + std::to_string(hist[i]) + "\n";
//    }

//    out << o;
//    out.close();

    return 0;
}

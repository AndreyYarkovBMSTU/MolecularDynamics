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


int main()
{
    // Создание системы частиц.
    Material* particlematerial = new SiliconDioxide("silica");
    Material* solvent = new Water("water");

    Properties* prop = new Properties(particlematerial, solvent);

    State* state1 = new State(1, Vector(-2.0, 1.0, 0.0) * phys::um() * 10.0, Vector(0.0, 0.0, 0.0) * phys::um() / phys::s());
    State* state2 = new State(2, Vector(2.0, -1.0, 0.0) * phys::um(), Vector(0.0, 0.0, 0.0) * phys::um() / phys::s());
//    State* state1 = new State(1, Vector(-4.0, -2.0, 0.0) * phys::um(), Vector(0.0, 0.0, 0.0) * phys::um() / phys::s());
//    State* state2 = new State(2, Vector(1.0, -4.0, 0.0) * phys::um(), Vector(0.0, 0.0, 0.0) * phys::um() / phys::s());
//    State* state3 = new State(3, Vector(2.0, 4.0, 0.0) * phys::um(), Vector(0.0, 0.0, 0.0) * phys::um() / phys::s());
//    State* state4 = new State(2, Vector(-4.0, 3.0, 0.0) * phys::um(), Vector(0.0, 0.0, 0.0) * phys::um() / phys::s());
//    State* state5 = new State(3, Vector(-4.0, 1.0, 0.0) * phys::um(), Vector(0.0, 0.0, 0.0) * phys::um() / phys::s());

    Particle* particle1 = new Dipoloid(state1, particlematerial);
    Particle* particle2 = new Dipoloid(state2, particlematerial);
//    Particle* particle3 = new Dipoloid(state3, particlematerial);
//    Particle* particle4 = new Dipoloid(state4, particlematerial);
//    Particle* particle5 = new Dipoloid(state5, particlematerial);

    //ExternalFields* externalfield = new DirectedField(prop, Vector(1e3, 0.0, 0.0));
    ExternalFields* externalfield = new RotatingField(prop, Vector(1e3, 0.0, 0.0), 30 * 1e3);
    Environment* environment = new Environment(solvent, externalfield);

    ParticleSystem* system = new ParticleSystem(environment, prop);
    system->setParticle(particle1);
    system->setParticle(particle2);
//    system->setParticle(particle3);
//    system->setParticle(particle4);
//    system->setParticle(particle5);

    // Расчёт

    Methods* method = new SelfConsistentDipoles(system, prop);

    method->setDipoleMoment();
    std::cout << system->particles[0]->dipolemoment;

    MolecularDinamic* moleculardinamic = new MolecularDinamic(system,
                                                              method,
                                                              prop,
                                                              "langevin",
                                                              "verle",
                                                              "LJ");

//    moleculardinamic->record();

    moleculardinamic->computer(2, 1000);

//    for (auto &iParticle : particlesystem->particles)
//    {
//        std::cout << iParticle->state;
//    }


//    std::string path = "output/gauss.txt";
//    std::ofstream out;
//    out.open(path, std::ios::out | std::ios::binary);
//    for (int i = 0; i < 1000; i++)
//    {
//        double t = -300 + std::rand() % 600;
//        double f = phys::getGaussian(t / 100, 0.0 , 1.0);
//        std::string o = std::to_string(t / 100) + "   " + std::to_string(f) + "\n";

//        out << o;
//    }
//    out.close();


    return 0;
}

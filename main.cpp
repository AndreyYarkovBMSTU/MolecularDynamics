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

//    State* state1 = new State(1, Vector(0.0, 0.0, 0.0), Vector(0.0, 0.0, 0.0));
//    State* state2 = new State(2, Vector(2.0, 0.0, 0.0), Vector(0.0, 0.0, 0.0));
    State* state1 = new State(1, Vector(-2.0, 0.0, 0.0), Vector(0.0, 0.0, 0.0));
    State* state2 = new State(2, Vector(0.0, 2.0, 0.0), Vector(0.0, 0.0, 0.0));
    State* state3 = new State(3, Vector(2.0, 0.0, 0.0), Vector(0.0, 0.0, 0.0));
//    State* state4 = new State(4, Vector(0.0, -2.0, 0.0), Vector(0.0, 0.0, 0.0));
//    State* state5 = new State(5, Vector(-6.0, 1.0, 0.0), Vector(0.0, 0.0, 0.0));
//    State* state6 = new State(6, Vector(-1.1, 4.0, 0.0), Vector(0.0, 0.0, 0.0));
//    State* state7 = new State(7, Vector(1.1, 2.0, 0.0), Vector(0.0, 0.0, 0.0));
//    State* state8 = new State(8, Vector(4.0, 4.0, 0.0), Vector(0.0, 0.0, 0.0));
//    State* state9 = new State(9, Vector(-4.0, 1.0, 0.0), Vector(0.0, 0.0, 0.0));
//    State* state10 = new State(10, Vector(-4.0, 6.0, 0.0), Vector(0.0, 0.0, 0.0));

    Particle* particle1 = new Dipoloid(state1, particlematerial, obj);
    Particle* particle2 = new Dipoloid(state2, particlematerial, obj);
    Particle* particle3 = new Dipoloid(state3, particlematerial, obj);
//    Particle* particle4 = new Dipoloid(state4, particlematerial, obj);
//    Particle* particle5 = new Dipoloid(state5, particlematerial, obj);
//    Particle* particle6 = new Dipoloid(state6, particlematerial, obj);
//    Particle* particle7 = new Dipoloid(state7, particlematerial, obj);
//    Particle* particle8 = new Dipoloid(state8, particlematerial, obj);
//    Particle* particle9 = new Dipoloid(state9, particlematerial, obj);
//    Particle* particle10 = new Dipoloid(state10, particlematerial, obj);

    //ExternalFields* externalfield = new DirectedField(prop, Vector(1.0, 0.0, 0.0));
    ExternalFields* externalfield = new RotatingField(prop, Vector(1.0, 0.0, 0.0), 30 * 1e3);
    Environment* environment = new Environment(solvent, externalfield);

    ParticleSystem* system = new ParticleSystem(environment, prop);
    system->setParticle(particle1);
    system->setParticle(particle2);
    system->setParticle(particle3);
//    system->setParticle(particle4);
//    system->setParticle(particle5);
//    system->setParticle(particle6);
//    system->setParticle(particle7);
//    system->setParticle(particle8);
//    system->setParticle(particle9);
//    system->setParticle(particle10);
    system->setProperties();

    // Расчёт

    Methods* method = new SelfConsistentDipoles(system, prop);

    method->setDipoleMoment();
    std::cout << system->particles[0]->dipolemoment << std::endl;

    std::cout << "T / timestep: " << (2 * M_PI / system->environment->externalfield->omega) / prop->timestep << std::endl;

    MolecularDinamic* moleculardinamic = new MolecularDinamic(system,
                                                              method,
                                                              prop,
                                                              "brownian",//"langevin"
                                                              "verle",
                                                              "LJ",
                                                              "average"); //"average"

    moleculardinamic->record(prop->path + "MD/MDtripleAver");

//    Interaction* interaction = new Interaction(system, method, prop, "average");
//    interaction->recordPotentials(1000, "induction");
//    interaction->recordPotentials(1000, "interaction");
//    interaction->recordPotentials(1000, "self");
//    interaction->recordPotentials(1000, "ipl3");

//    moleculardinamic->recordVMD(prop->path + "vmd/vmdTriple");
//    moleculardinamic->recordtest(prop->path + "MD/MDtest");

    std::cout << "First modeling finished" << std::endl;

    externalfield->omega = externalfield->omega / (moleculardinamic->t0 * prop->koef_omega); // возвращение значения

//    system->particles[0]->setCoordinate(Vector(-2, 0.0, 0.0));
//    system->particles[0]->setVelocity(Vector(0.0, 0.0, 0.0));
//    system->particles[1]->setCoordinate(Vector(0.0, 2.0, 0.0));
//    system->particles[1]->setVelocity(Vector(0.0, 0.0, 0.0));
//    system->particles[2]->setCoordinate(Vector(2.0, 0.0, 0.0));
//    system->particles[2]->setVelocity(Vector(0.0, 0.0, 0.0));
//    system->particles[3]->setCoordinate(Vector(0.0, -2.0, 0.0));
//    system->particles[3]->setVelocity(Vector(0.0, 0.0, 0.0));

    method->setDipoleMoment();

//    MolecularDinamic* moleculardinamic1 = new MolecularDinamic(system,
//                                                              method,
//                                                              prop,
//                                                              "brownian",//"langevin"
//                                                              "verle",
//                                                              "LJ",
//                                                              "average"); //"exact"
//    moleculardinamic1->record(prop->path + "time1/");

//    std::cout << "Second modeling finished" << std::endl;

//    prop->koef_dipole = 1e1;
//    moleculardinamic->record(prop->path + "time2/");

//    system->particles[0]->setCoordinate(Vector(-1.5, 0.0, 0.0));
//    system->particles[0]->setVelocity(Vector(0.0, 0.0, 0.0));
//    system->particles[1]->setCoordinate(Vector(1.5, 0.0, 0.0));
//    system->particles[1]->setVelocity(Vector(0.0, 0.0, 0.0));
//    system->particles[2]->setCoordinate(Vector(2.0, 4.0, 0.0));
//    system->particles[2]->setVelocity(Vector(0.0, 0.0, 0.0));

//    prop->koef_dipole = 1e2;
//    moleculardinamic->record(prop->path + "time3/");

    std::cout << "Reinolds: " << system->reinolds << std::endl;
    std::cout << "v_thermal: " << system->v_thermal << std::endl;
    std::cout << "tau: " << moleculardinamic->tau << std::endl;
    std::cout << "t0: " << moleculardinamic->t0 << std::endl;
    std::cout << "m/gamma: " << system->particles[0]->mass / system->friction << std::endl;
    std::cout << "D: " << system->diffusion << std::endl;
    std::cout << "dipolemoment0: " << system->dipolemoment0.norm() << std::endl;
    std::cout << "dipolemoment0_: " << system->dipolemoment0_.norm() << std::endl;
    std::cout << "ksi: " << system->ksi << std::endl;

//    moleculardinamic->computer(0, 10000);

//    for (int numfile = 0; numfile < 100; numfile++)
//    {
//        moleculardinamic->recordmove(numfile, 100000);
//        system->particles[0]->setCoordinate(Vector(-1.1, 0.0, 0.0));
//        system->particles[0]->setVelocity(Vector(0.0, 0.0, 0.0));
//    }

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

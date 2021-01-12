//
//  main.cpp
//  PSO
//
//  Created by Simon Peter on 05/01/21.
//

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//PSO parameters
//Particles are like pigeons / entities
//Iterations - the time available to arrive at the solution
#define MAX_PARTICLES           100
#define MAX_ITERATIONS          100
#define INITIAL_FITNESS         -750

#define INERTIA_WEIGHT            0.529f;
#define LOCAL_WEIGHT            1.49445f; /* cognitive/local weight  */
#define GLOBAL_WEIGHT            1.49445f; /* social/global weight    */


#define POS_COORD               2
#define VEL_COORD               2

typedef struct _pos {
    double p[POS_COORD];
} Position;

typedef struct _vel {
    double p[VEL_COORD];
} Velocity;

typedef struct _particle {
    Velocity velocity;               //Particle velocity
    Position position;                 //current position
    Position bestPosition;           //best position history
    double   fitness;                //current fitness
    double   bestFitness;            //best fitness history
    int index;                       //for debugging - particle id
} Particle;

typedef struct _swarm {
    Particle *members;               //The members of the Swarm
    double   bestGlobalFitness;      //Best fit particle in the swarm
    Position bestGlobalPosition;     //Position of the above particle
    int      population;             //Swarm population
} Swarm;


Particle *InitParticles(std::string FileName, int numberParticles,
            double *bestGlobalFitness, Position *bestGlobalPosition);
void ProcessPSO(std::string FileName, Swarm *pSwarm, int numberIterations);
void DisplayParam(Swarm *pSwarm, std::string qualifier);
void DisplayResults(Swarm *pSwarm);

//Some constants.
double minX = 1.0;
double maxX = 100.0;
double maxV = 5, minV = 1;


/*
 * randomNum
 *   Generates a random number (- (-1, 1)
 */
static double randomNum(void)
{
    double val;
    double v1, v2;
    
    v1 =  rand();
    v2 = 1.1 + RAND_MAX;
    val = -1.0 + (v1 / v2) * 2.0 ;
    return val;
}

/*
 * ObjectiveFunction
 *   Model of the problem statement.. Write your own according to the
 * problem domain.
 */
static double ObjectiveFunction(std::string FileName, Position *pPos)
{
    // f(x) = 3 + x^2 + y^2
    int x;
    int y;
    //double acc = 0;
    
    x = (int) pPos->p[0];
    y = (int) pPos->p[1];
    
    //if (x < 1)
     //   x = 1;

    //acc = donn(FileName,x);
    //return acc;

#if 1
    if (x > 50)
        return 1000-10*(x-50);
    else
        return 1000-10*(50-x);
#else
    return 50-x + 2*x;
#endif
}

void pso (std::string FileName)
{
    Swarm ThisSwarm;
    Position guessBest;

    int numberParticles  = MAX_PARTICLES;
    int numberIterations = MAX_ITERATIONS;

    //Initialize the RNG..
    srand((unsigned)time(NULL));

    printf("\n --Begin PSO demo v1.0-----------------------------\n");
    printf("\nObjective function is f(x) = 3 + (x0^2 + x1^2)");
    printf("\nParticles = %d, Iterations = %d", numberParticles, numberIterations);
    printf("\n -------------------------------------------------\n");
    guessBest.p[0] = 100.0;

    ThisSwarm.bestGlobalPosition = guessBest;
    ThisSwarm.bestGlobalFitness = INITIAL_FITNESS;
    ThisSwarm.population = numberParticles;
    ThisSwarm.members = InitParticles(FileName, numberParticles, &(ThisSwarm.bestGlobalFitness), &(ThisSwarm.bestGlobalPosition));

    if (ThisSwarm.members == NULL) {
        printf("\n Unable to create the swarm. Exiting..");
        exit(-1);
    }

    DisplayParam(&ThisSwarm, "Initial");
    
    // Main processing loop
    ProcessPSO(FileName, &ThisSwarm, numberIterations);
    
    // Display results
    DisplayResults(&ThisSwarm);

    //getchar();
}

void DisplayParam(Swarm *pSwarm, std::string qualifier)
{
    int i;

    printf("\n %s best fitness = %1.4f", qualifier.c_str(), pSwarm->bestGlobalFitness);
    printf("\n %s Best position/solution:", qualifier.c_str());
    for (i = 0; i < POS_COORD; ++i) {
        printf("\n x[%d] : %2.4f", i, pSwarm->bestGlobalPosition.p[i]);
    }
}

void DisplayResults(Swarm *pSwarm)
{
    printf("\nProcessing complete");
    printf("\n------------------------------------------------------");
    DisplayParam(pSwarm, "Final");
    printf("\n------------------------------------------------------");
    printf("\n ######### End PSO demonstration ###############\n");
}

Particle *InitParticles(std::string FileName, int numberParticles,
            double *bestGlobalFitness, Position *bestGlobalPosition)
{
    Particle *member, *pPart;
    Position randomPosition;
    Velocity randomVelocity;
    double lo, hi, fitness, minV, maxV;
    int i, j;

    minV = 1; maxV = 5;
    
    //Allocated the number of particles..
    member = (Particle *) malloc (sizeof(Particle) * numberParticles);
    if (NULL == member) {
        //Failed to allocate members..
        return NULL;
    }
    
    for (i = 0; i < numberParticles; ++i) {
        pPart = &member[i];
        pPart->index = i;
        lo = minX * 1.0;
        hi = maxX * 1.0;
        for (j = 0; j < POS_COORD; ++j) {
            randomPosition.p[j] = ((hi - lo) * randomNum()*0.8) ;
            if (randomPosition.p[j] < minX)
                randomPosition.p[j] = minX;
            else if (randomPosition.p[j] > maxX)
                randomPosition.p[j] = maxX;
        }

        fitness = ObjectiveFunction(FileName, &randomPosition);
        pPart->position = pPart->bestPosition = randomPosition;
        pPart->fitness = pPart->bestFitness = fitness;
        
        if (fitness > *bestGlobalFitness) {
            *bestGlobalFitness = fitness;
            *bestGlobalPosition = randomPosition;
        }
    
        for (j = 0; j < VEL_COORD; ++j) {
            hi = (maxX > minX) ? (maxX - minX) : (minX - maxX);
            lo = -hi * 1.0;
            randomVelocity.p[j] = (hi - lo) * randomNum() + lo;

            if (randomVelocity.p[j] < minV)
                randomVelocity.p[j] = minV;
            else if (randomVelocity.p[j] > maxV)
                randomVelocity.p[j] = maxV;

        }
        pPart->velocity = randomVelocity;
    }
    
    return member;
    
}

/*
 * ProcessPSO : Apply the Particle Swarm Optimization algo.
 * see http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=00870279
 */
void ProcessPSO(std::string FileName, Swarm *pSwarm, int numberIterations)
{
    Particle *member;
    Velocity newVelocity;
    Position newPosition;
    Position bestGlobalPos;

    double w, c1, c2;
    double newFitness, bestGlobalFitness;
    double r1, r2; // cognitive and social randomizations
    int i, j, iteration, total_particles;

    w = INERTIA_WEIGHT;
    c1 = LOCAL_WEIGHT;
    c2 = GLOBAL_WEIGHT;
    printf("\nEntering main PSO processing loop");
    total_particles = pSwarm->population;
    bestGlobalPos = pSwarm->bestGlobalPosition;
    bestGlobalFitness = pSwarm->bestGlobalFitness;
    member = pSwarm->members;
    for (iteration = 0; iteration < numberIterations; ++iteration) {
        // each Particle
        printf("\nIteration: %d ",iteration+1);
        for (i = 0; i < total_particles; ++i)  {
            //member[i] = member[i];

            // each x value of the velocity
            for (j = 0; j < VEL_COORD; ++j) {
                r1 = randomNum();
                r2 = randomNum();

                newVelocity.p[j] = (double)(w * member[i].velocity.p[j]) +
                (c1 * r1 * (member[i].bestPosition.p[j] - member[i].position.p[j])) +
                (c2 * r2 * (bestGlobalPos.p[j] - member[i].position.p[j]));

                if (newVelocity.p[j] < minV)
                    newVelocity.p[j] = minV;
                else if (newVelocity.p[j] > maxV)
                    newVelocity.p[j] = maxV;
            }

            member[i].velocity = newVelocity;
            
            //For all position coordinates..
            for (j = 0; j < POS_COORD; ++j) {
                double dir = randomNum();
                if (dir > 0)
                    newPosition.p[j] = member[i].position.p[j] + newVelocity.p[j];
                else
                    newPosition.p[j] = member[i].position.p[j] - newVelocity.p[j];
                //Normalize the coordinates..
                if (newPosition.p[j] < minX)
                    newPosition.p[j] = minX;
                else if (newPosition.p[j] > maxX)
                    newPosition.p[j] = maxX;
            }

            member[i].position = newPosition;
            
            newFitness = ObjectiveFunction(FileName, &newPosition);
            member[i].fitness = newFitness;
            if (newFitness > member[i].bestFitness) {
                member[i].bestPosition = newPosition;
                printf("\n\tParticle:%d got new fitness: %f > %f",i+1,newFitness,member[i].bestFitness);
                member[i].bestFitness = newFitness;
            }

            if (newFitness > bestGlobalFitness) {
                bestGlobalPos = newPosition;
                printf("\n!! Global fitness: %f > %f !!",newFitness, bestGlobalFitness);
                bestGlobalFitness = newFitness;
            }
            
        } // each Particle

    } // interations..
    pSwarm->bestGlobalPosition = bestGlobalPos;
    pSwarm->bestGlobalFitness = bestGlobalFitness;
}


int main(int argc, const char * argv[]) {
    pso("annpso.csv");
    return 0;
}

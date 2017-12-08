//============================================================================
//  Created on : 2014-08-19
//      Author : Mohammad Nazrul Ishlam Patoary
// Version     :
// Copyright   : Your copyright notice
// Description : CaWave in C++, Ansi-style
//============================================================================
//# include "mpi.h"
#include "Main.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <time.h>
#include <fstream>
#include <functional>
#include <iostream>
#include <list>
#include <queue>
#include <sstream>
#include <string>
#include <vector>
#include <queue>
#include "Heap.h"
#include "SV.h"
#include "species.h"

using namespace std;

const double CONVERSION_FACTOR = 602214.129;
const double CONVERSION_FACTOR_VOLUME = CONVERSION_FACTOR / 64.0;

vector<double> getTotalStates(int numberSpecies);
void REACTION(SV* topSV, double uniRand);
bool LT(long double simTime, long double nextSimTime);
// void mpiInitialization();
// void mpiFinalization();

double simulationTime = 0;
struct timeval start_time, end_time;
struct rusage cluster_start, cluster_end, cluster_mid;
long double cl_start_time = 0;
long double cl_end_time = 0;
long double sim_time = 0;
long double sget_time = 0;
unsigned long simulation_time = 0;  // to find out the simulation time .
int testCounter = 0;

int SV_size;
int CaWaveFlag;  //*********
double* diffusionContant;
string speciesNames[50];
double diffusionConstants[50];
double uniRand1;

vector<reaction*> reactions;
vector<species*> spcs;
int numberOfSpecies;
int numberOfReactions;
int globalReactionCounter = 0;
reaction* globalReaction[50];

vector<vector<int> > listOfNeighbours;  // list of neighbor lists
vector<vector<int> > spatialLocation;   // dx, dy dz
vector<SV*> sv;
// vector<SV*> sv_init;

double p = 0;
SV* topSV = 0;

// int node;
int dataIterator = 0;
int diffusionCounter = 0;
int reactionCounter = 0;
double* states;
// int lastID;

int initialized, finalized, node;

// ofstream myfile;

struct Compare_pqTau {
    bool operator()(SV* p1, SV* p2) {
        return p1->pqTau > p2->pqTau;
    }
};

void MAKE_HEAP(vector<SV*> sv) {
  std::make_heap(const_cast<SV**>(&sv[0]),
                 const_cast<SV**>(&sv[0]) + sv.size(),
                 Compare_pqTau());
}

void REACTION(SV* topSV, double uniRand) {
  int index = topSV->getMiu(uniRand);
  int reactionDirection = topSV->reactionDirection;
  // for reactant . .
  bool flag = true;  //

  for (int i = 0;
       i < (int)(topSV->Reactions[index]->reactantMetricValue).size(); i++) {
    if (topSV->Species[topSV->Reactions[index]->reactantSpeciesID[i]]
            ->numberOfMolecules <= 0) {
      flag = false;
    }
  }

  if (flag) {
    // for reactant ..
    for (int i = 0;
         i < (int)(topSV->Reactions[index]->reactantMetricValue).size(); i++) {
      int n = topSV->Reactions[index]->reactantMetricValue[i];
      int m = topSV->Reactions[index]->reactantSpeciesID[i];
      // cout << n << " n; m: " << m << endl;
      topSV->Species[m]->numberOfMolecules -= n * reactionDirection;
      // states[(m * SV_size) + topSV->ID] -= n * reactionDirection;
      //##################################
      // cout<<topSV->Species[m]->numberOfMolecules<<" =
      // "<<states[(m*SV_size)+topSV->ID]<<endl;
    }
    // for Product . .
    for (int i = 0; i < (int)topSV->Reactions[index]->productMetricValue.size();
         i++) {
      int n = topSV->Reactions[index]->productMetricValue[i];
      int m = topSV->Reactions[index]->productSpeciesID[i];
      topSV->Species[m]->numberOfMolecules += n * reactionDirection;
      // states[(m * SV_size) + topSV->ID] += n * reactionDirection;
      // ##################################
      // cout<<topSV->Species[m]->numberOfMolecules<<" =
      // "<<states[(m*SV_size)+topSV->ID]<<endl;
    }
    topSV->setTau(simulationTime);  // setTau();//Tau is updated of the SV
    reactionCounter++;
  } else {
    topSV->setTau(simulationTime);
  }
}

bool LT(long double simTime, long double nextSimTime) {
  double epsilon = .0000000001;
  if (simTime < 0 || nextSimTime < 0) {
    cout << "simTime: " << simTime << endl;
    cout << "nextSimTime: " << nextSimTime << endl;
    cout << "[ERROR] less than 0" << endl;
  }

  return (nextSimTime - simTime) >
         ((fabs(simTime) < fabs(nextSimTime) ? fabs(nextSimTime)
                                             : fabs(simTime)) *
          epsilon);
}

// Binary min heap. . . .
Heap* pq = new Heap(SV_size);

// NTW function 1
/*
void mpiInitialization(){
        //MPI initializaiton !
        //int initialized, finalized,node;

        MPI_Initialized(&initialized);
        if (!initialized){
                MPI_Init(NULL, NULL);
        }
        //srand((unsigned)time(0));
        //int node;
        //srand(10);
        MPI_Comm_rank(MPI_COMM_WORLD, &node);

}
*/

void setRandomSeed(int s) { srand(s); }

void getInformationsFromPython(int sv_size, vector<vector<int> > lon,
                               double* st, int numOfSpecies, double* d,
                               int waveFlag) {
  // Getting neighbour information !
  listOfNeighbours = lon;
  // Getting total number of sub-volumes !
  SV_size = sv_size;
  cout << "Total SVs: " << SV_size << endl;
  Heap* pq1 = new Heap(SV_size);
  pq = pq1;
  states = st;
  diffusionContant = d;
  numberOfSpecies = numOfSpecies;
  CaWaveFlag = waveFlag;

  // cout<<CaWaveFlag<<"  diffusion constant  "<<diffusionContant[0]<<"
  // "<<diffusionContant[1]<<" "<<diffusionContant[2]<<endl;  exit(-1);
}

void createReactions(vector<double> putreactantMetricValue,
                     vector<int> putreactantSpeciesID,
                     vector<double> putproductMetricValue,
                     vector<int> putproductSpeciesID,
                     double putreactionConstant, double dx_,
                     double (*getReactionRate_C)(species* Species[])) {
  // cout<<"Reach here . . . "<<endl;
  // func();
  globalReaction[globalReactionCounter] = new reaction(globalReactionCounter);
  globalReaction[globalReactionCounter]->getReactionRate_C = getReactionRate_C;
  // assigning reaction values to reaction object
  globalReaction[globalReactionCounter]->reactantMetricValue =
      putreactantMetricValue;
  globalReaction[globalReactionCounter]->reactantSpeciesID =
      putreactantSpeciesID;

  if (CaWaveFlag == 0) {
    globalReaction[globalReactionCounter]->reactionConstant =
        (double)(putreactionConstant);
  } else {
    // cout<<"REAched "<<CaWaveFlag<<endl;exit(-1);
    globalReaction[globalReactionCounter]->reactionConstant =
        (double)(putreactionConstant);
  }
  // cout<<factor<<"   "<<globalReactionCounter<<endl;

  globalReaction[globalReactionCounter]->productMetricValue =
      putproductMetricValue;
  globalReaction[globalReactionCounter]->productSpeciesID = putproductSpeciesID;
  reactions.push_back(globalReaction[globalReactionCounter]);  //
  globalReactionCounter++;
}

void initSVs(double dx_) {
  // myfile.open ("./Data/Stochastic/Time.txt");
  // necessary to resize
  sv.resize(SV_size);
  numberOfReactions = globalReactionCounter;

  getrusage(RUSAGE_SELF, &cluster_start);
  // gettimeofday(&start_time, NULL);
  for (int i = 0; i < numberOfSpecies; i++) {
    diffusionConstants[i] = (double)diffusionContant[i];
  }

  // 5. CREATING subvolumes and pushing to priority quequ, pq
  SV* Sv;
  double dx = dx_;  //.250;//as I get from NEURON

  for (int i = 0; i < SV_size; i++) {
    // bool flag = true;
    // 5.1 create SV object and put into sv[] vector
    Sv = new SV(i, numberOfSpecies, numberOfReactions);  // id ta 1 theke suru
    sv[i] = Sv;

    sv[i]->setCaWaveFlag(CaWaveFlag);

    // 5.2 create of list of species . . . and assign to sv[i]
    for (int s = 0; s < numberOfSpecies; s++) {
      species* spc = new species(s);

      spc->diffusionConstant = diffusionConstants[s];
      spcs.push_back(spc);
    }
    sv[i]->Species = spcs;  // providing species
    spcs.clear();

    // 5.3 assign reactions to sv[i]
    for (int r = 0; r < (int)reactions.size(); r++) {
      reaction* rect = new reaction(reactions[r]->reactionID);
      rect->reactantMetricValue = reactions[r]->reactantMetricValue;
      rect->reactantSpeciesID = reactions[r]->reactantSpeciesID;
      rect->reactionConstant = reactions[r]->reactionConstant;
      rect->productMetricValue = reactions[r]->productMetricValue;
      rect->productSpeciesID = reactions[r]->productSpeciesID;
      rect->getReactionRate_C = reactions[r]->getReactionRate_C;
      sv[i]->Reactions.push_back(rect);
    }

    /*
    if(i==0){
    //To display your reaction set !
    cout<<"\n\tTHE REACTION SET "<<endl;
    for(int ii=0;ii<sv[i]->Reactions.size();ii++){
     cout<<"R["<<ii<< "],  ";

     for(int j=0;j<sv[i]->Reactions[ii]->reactantMetricValue.size();j++){
                    cout<<""<<sv[i]->Reactions[ii]->reactantMetricValue[j]<<"*"<<sv[i]->Reactions[ii]->reactantSpeciesID[j]<<"
    + ";
     }
     cout<<" ====>  ";
     for(int k=0;k<(int)sv[i]->Reactions[ii]->productMetricValue.size();k++){
                    cout<<""<<sv[i]->Reactions[ii]->productMetricValue[k]<<"*"<<sv[i]->Reactions[ii]->productSpeciesID[k]<<"
    + ";
     }
     cout<<"\t, Rate: "<<sv[i]->Reactions[ii]->reactionConstant;
     cout<<endl;
    }
    }//
    */

    if (CaWaveFlag == 0) {
      // NEURON's initial data is assigned here !
      for (int s = 0; s < numberOfSpecies; s++) {
        // cout<<states[((s%numberOfSpecies)*SV_size)+i]<<endl;
        // states[((s % numberOfSpecies) * SV_size) + i] =
        //     (long)(CONVERSION_FACTOR_VOLUME *
        //            states[((s % numberOfSpecies) * SV_size) + i] +
        //            0.00000001);
        // cout << "states" << states[((s % numberOfSpecies) * SV_size) + i] <<
        // endl;
        double tmp = CONVERSION_FACTOR_VOLUME *
                     states[((s % numberOfSpecies) * SV_size) + i];
        sv[i]->Species[s]->numberOfMolecules = (long)(tmp);

        double rand_0_1 = (double)rand() / (double)RAND_MAX;
        if (rand_0_1 < tmp - sv[i]->Species[s]->numberOfMolecules) {
          sv[i]->Species[s]->numberOfMolecules += 1;
        }
        // cout << rand_0_1 << " " << tmp << " " <<
        // sv[i]->Species[s]->numberOfMolecules << endl;

        states[((s % numberOfSpecies) * SV_size) + i] =
            (double)(sv[i]->Species[s]->numberOfMolecules);
        // cout << states[((s % numberOfSpecies) * SV_size) + i] << endl;
      }
    } else {
      for (int s = 0; s < numberOfSpecies; s++) {
        // cout<<states[((s%numberOfSpecies)*SV_size)+i]<<endl;
        sv[i]->Species[s]->numberOfMolecules =
            (long)(states[((s % numberOfSpecies) * SV_size) + i]);
      }
    }
    // 5.5 provide neighbour information to sv[i]
    sv[i]->Neighbours = listOfNeighbours[i];

    // 5.6 provide surfaceArea, volume, area and distance information
    // sv[i]->Location = spatialLocation[i];
    sv[i]->surfaceArea = dx * dx;    // surfaceArea, from NEURON
    sv[i]->volume = dx * dx * dx;    // volume, from NEURON
    sv[i]->area = dx * dx;           // neighbor area, from NEURON
    sv[i]->distance = (double)(dx);  // /2.175;

    // cout << dx << endl; exit(-1);
    // 5.7 compute initial Tau and others for sv[i]
    sv[i]->setTau(simulationTime);
    // 5.8 push sv[i] to priority queue pq
    pq->push(sv[i]);
  }
  // pq->Display();
}

void mainLoop(double nextTimeStep) {
  // cout<<"  Simulation is going on . . . ."<<endl;
  // cout<<"  Time step is : "<<nextTimeStep<<endl;
  // if(nextTimeStep==.002){
  // cout<<"First Fired SV: "<<pq.top()->ID<<endl;
  //}
  // if(node==0){

  // get simulation time t from NEURON . . .
  double getNextSimulationTime = nextTimeStep;
  // cout<<"First Fired SV: "<<pq.top()->ID<<endl;
  while (simulationTime < getNextSimulationTime) {
    topSV = pq->pop();
    simulationTime = topSV->Gt;
    // cout<<simulationTime<<endl;
    //    myfile<< simulationTime<<"\n";
    // simulationTime = topSV->Lt;
    // simulationTime = topSV->Tau + simulationTime;
    // pq->Display();

    uniRand1 = (double)rand() / (double)RAND_MAX;
    // cout << topSV->r << "; " << topSV->s << endl;
    // exit(-1);
    double tmp = (double)(topSV->r / ((double)(topSV->r + topSV->s)));

    if (topSV->s > 0) {
      if (uniRand1 < tmp) {  // NSM basic condition !
        REACTION(topSV, uniRand1);
        pq->push(topSV);
        reactionCounter++;

      } else {
        // double tmp1 = topSV->r / (topSV->r + topSV->s);
        // double tmp2 = topSV->s / (topSV->r + topSV->s);
        // double uniRand3 = (uniRand1 - tmp1) / tmp2;
        // int index = topSV->getdMiu(uniRand3);
        // Just re-use uniRand1 since this parameter is not used in getdMiu
        int index = topSV->getdMiu(uniRand1);
        int randomFace = (int)(rand() % topSV->Neighbours.size());

        // Diffusion
        int targetID = topSV->Neighbours[randomFace];

        if (topSV->ID == targetID) {
          cout << "[ERROR] Connectivity" << endl;
        }
        if (topSV->Species[index]->numberOfMolecules == 0) {
          cout << "[ERROR] numberOfMolecules is 0" << endl;
        }

        topSV->Species[index]->numberOfMolecules -= 1;
        // states[index * SV_size + topSV->ID] -= 1;
        topSV->setTau(simulationTime);
        pq->push(topSV);

        sv[targetID]->Species[index]->numberOfMolecules += 1;
        // states[index * SV_size + targetID] += 1;
        sv[targetID]->setTau(simulationTime);
        pq->update(targetID);

        diffusionCounter++;
      }
    } else {
      REACTION(topSV, uniRand1);
      pq->push(topSV);
      reactionCounter++;
    }
  }

  for (int i = 0; i < SV_size; i++) {
    for (int s = 0; s < numberOfSpecies; s++) {
      states[s * SV_size + i] = sv[i]->Species[s]->numberOfMolecules;
    }
  }

  // myfile.close();
  // cout << "nextTimeStep: " << nextTimeStep << endl;
  // exit(-1);
  /*
  int token=99;
  MPI::COMM_WORLD.Send(&token, 1, MPI::INT, 1, 1);
  MPI::COMM_WORLD.Recv(&token, 1, MPI::INT, 1, 1);
  cout << " P0 receiving Token value :" << token << endl;
  */

  //}//if node = 0

  // on test
  // if(node==1)  {
  /*
          int p1_rc_token = 0;
          MPI::COMM_WORLD.Recv(&p1_rc_token, 1, MPI::INT, 0, 1);
          cout << " P1 receiving Token value :" << p1_rc_token << endl;
          p1_rc_token--;
          MPI::COMM_WORLD.Send(&p1_rc_token, 1, MPI::INT, 0, 1);
  */

  //}

  // To print all SVs states in a file
  /*  ofstream myFile;
          myFile.open("all_SV_states.txt");
          for(int i=0; i<SV_size;i++){
                  //sum = sum + sv[i]->Species[0]->numbeOfMolecules;
                  myFile<<sv[i]->Species[0]->numberOfMolecules<<"\n";
          }
          myFile.close();
  */
}

void statistics() {
  // myfile.close();

  getrusage(RUSAGE_SELF, &cluster_end);
  gettimeofday(&end_time, NULL);
  simulation_time = time(0) - simulation_time;  // get simulation interval
  cl_start_time = (long double)cluster_start.ru_stime.tv_sec +
                  (long double)cluster_start.ru_stime.tv_usec / 1000000 +
                  (long double)cluster_start.ru_utime.tv_sec +
                  (long double)cluster_start.ru_utime.tv_usec / 1000000;
  cl_end_time = (long double)cluster_end.ru_stime.tv_sec +
                (long double)cluster_end.ru_stime.tv_usec / 1000000 +
                (long double)cluster_end.ru_utime.tv_sec +
                (long double)cluster_end.ru_utime.tv_usec / 1000000;
  sim_time = (cl_end_time - cl_start_time);
  sget_time =
      ((long double)end_time.tv_sec - (long double)start_time.tv_sec) +
      ((long double)end_time.tv_usec - (long double)start_time.tv_usec) /
          1000000;
  cout << "Simulation time :" << sim_time << " sec." << endl;
  cout << "Consumed memory :" << (cluster_end.ru_maxrss) << " kilobytes"
       << endl;
  cout << "Reaction :" << reactionCounter << " Diffusion :" << diffusionCounter
       << endl;
}

/*
void mpiFinalization(){
        //MPI finalization !
        MPI_Finalized(&finalized);

        //if (!finalized){
        //  MPI_Finalize();
        //}
        cout<<"Finalized "<<endl;

}
*/

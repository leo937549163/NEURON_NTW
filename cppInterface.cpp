/*
 * SV.cpp
 *
 *  Created on: 2014-08-19
 *      Author: Mohammad Nazrul Ishlam Patoary
 *
 */

#include "cppInterface.h"
#include <functional>
#include <iostream>
#include <vector>
#include "Main.h"
#include "species.h"

using namespace std;

vector<vector<int> > listFromNeuron;
std::vector<int> Neighbour;

myClass::myClass() {}
myClass::~myClass() {}

void myClass::cpp_function() { std::cout << "Simulator is called here  \n"; }

void *get_myClass() {
  myClass *out = new myClass();
  return ((void *)out);
}

void Send_to_cppInterface(void *objectPtr, int sv_size, long (*neighbours)[6],
                          double *states, int numberOfSpecies, double *d,
                          int CaWaveFlag) {
  for (int i = 0; i < sv_size; i++) {
    for (int j = 0; j < 6; j++) {
      if (neighbours[i][j] >= 0) {
        Neighbour.push_back((int)neighbours[i][j]);
      }
    }
    listFromNeuron.push_back(Neighbour);
    Neighbour.clear();
  }

  //   mpiInitialization();//***1. for without mpi

  // this funciton will be implemented in NTW to get the following informations
  // from NEURON
  getInformationsFromPython(sv_size, listFromNeuron, states, numberOfSpecies, d,
                            CaWaveFlag);  // NTW function
  // 1. sv_size = Number of SVs
  // 2. listFromNeuron = neighbor list
  // 3. status = is the pointer of all species in a same list = rxd.node._states
  // 4. numberOfSpecies = total number of species
  // 5. d = pointer of list of diffusion constants
}

void createReaction_ccpInterface(
    int numberOfReactant, double *reactantMetricValue, long *reactantSpeciesID,
    int numberOfProduct, double *productMetricValue, long *productSpeciesID,
    double *reactionConstant, double *dx,
    double (*getReactionRate_C)(species *Species[])) {
  vector<double> putreactantMetricValue;
  vector<int> putreactantSpeciesID;
  double putreactionConstant = (double)reactionConstant[0];
  vector<double> putproductMetricValue;
  vector<int> putproductSpeciesID;

  for (int i = 0; i < numberOfReactant; i++) {
    // cout<<reactantMetricValue[i]<<endl;
    putreactantMetricValue.push_back((double)reactantMetricValue[i]);
  }
  for (int i = 0; i < numberOfReactant; i++) {
    // cout<<reactantSpeciesID[i]<<endl;
    putreactantSpeciesID.push_back((int)reactantSpeciesID[i]);
  }

  for (int i = 0; i < numberOfProduct; i++) {
    // cout<<productMetricValue[i]<<endl;
    putproductMetricValue.push_back((double)productMetricValue[i]);
  }
  for (int i = 0; i < numberOfProduct; i++) {
    // cout<<productSpeciesID[i]<<endl;
    putproductSpeciesID.push_back((int)productSpeciesID[i]);
  }
  // cout<<"REACH here . . ."<<*(dx)<<endl;
  // NTW function to pass the following values:
  createReactions(putreactantMetricValue, putreactantSpeciesID,
                  putproductMetricValue, putproductSpeciesID,
                  putreactionConstant, (double)dx[0], getReactionRate_C);
  // 1. putreactantMetricValue, putreactantSpeciesID, putproductMetricValue,
  // putproductSpeciesID, putreactionConstant
}  //

void initSVs_cppInterface(double *dx) {
  // NTW function to get:
  initSVs((double)dx[0]);
  // 1. dx = length of SVs length, height and width
}

void runMainLoop_ccpInterface(double *nextTimeStep) {
  // NTW function :
  mainLoop((double)(nextTimeStep[0]));
  // 1. nextTimeStep = how long simulation will run
}

void finalization_cppInterface() {
  // NTW functions
  statistics();

  //  mpiFinalization();//***2
}

void setRandomSeed_cppInterface(int s) { setRandomSeed(s); }

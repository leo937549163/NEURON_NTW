/*
 * SV.cpp
 *
 *  Created on: 2014-08-19
 *      Author: Mohammad Nazrul Ishlam Patoary
 *
 */
#include "SV.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>

using namespace std;
SV::SV() {}
SV::~SV() {}

SV::SV(int id, int nSpecies, int nReaction) {
  numberOfOpen = 0;
  r = 0.0;   // reaction rate
  s = 0.0;   // diffusion rate
  Tau = 0;   // putative time
  Lt = 0.0;  // local time
  Gt = 0;
  pqTau = 0;
  ID = id;    // sv id
  Miu = -1;   // selected reaction
  dMiu = -1;  // selected species to diffuse
  numberOfReactions = nReaction;
  numberOfSpecies = nSpecies;  // lage na
  CaWaveFlag = 0;
  numberOfEvents = 0;
  reactionDirection = 1;
}

void SV::setReactionRates() {
  r = 0;
  // cout << "numberOfReactions: " << numberOfReactions << endl;
  for (int i = 0; i < numberOfReactions; i++) {
    if (CaWaveFlag == 0) {
      r = r + Reactions[i]->getReactionRate(Species);  // CA
    } else {
      // cout<<CaWaveFlag<<endl;
      r = r + Reactions[i]->getReactionRate_wave(Species);  // CA
    }
  }
  // rate which reaction happen
}  //

long double SV::getReactionRates() { return r; }  //

void SV::setDiffusionRates() {
  double directions = Neighbours.size();
  // number of directions in which the molecules can diffuse
  s = 0;
  for (int i = 0; i < numberOfSpecies; i++) {
    s += Species[i]->getDiffusionRate(directions, distance, surfaceArea, volume,
                                      CaWaveFlag);  // I fixed it to 6 aug 14
  }
}

long double SV::getDiffusionRates() { return s; }

void SV::updateRates() {
  setReactionRates();
  if (r < 0) {
  	r = -r; // let Reaction Rate >= 0
    reactionDirection = -1;
  } else {
    reactionDirection = 1;
  }
  setDiffusionRates();
}

void SV::setTau(double simTime) {
  updateRates();
  long double Total_Event_Rate = r + s;

  if (Total_Event_Rate <= 0) {
    cout << "[ERROR] Total Event Rate is 0" << endl;
  }

  if (Total_Event_Rate <= 0) {
    pqTau = (long double)RAND_MAX;
    return;
  }

  // double Urandom1 = (double)rand() / (double)RAND_MAX;
  // Tau = (long double)((long double)(1 / Total_Event_Rate) *
  //                     (long double)log(1 / (Urandom1)));
  Tau = (long double)((long double)(1 / Total_Event_Rate) *
                      (long double)log((double)RAND_MAX / (double)rand()));
  // cout<<"tau "<<Tau<<endl;
  Gt = simTime + Tau;
  // Lt = Lt + Tau;
  Lt = Gt;
  pqTau = (long double)Lt;
}

long double SV::getTau() { return Tau; }

// ##################################################### WORK here !
void SV::setMiu_new(double uniRand) {
  if (numberOfReactions == 1) {
    Miu = 0;  // when there is only one reaction
    return;
  }
  // Otherwise
  double Urandom1 = (double)rand() / (double)RAND_MAX;
  long double tmp = Urandom1 * r;
  long double sum = 0;
  for (int i = 0; i < numberOfReactions; i++) {
    if (LT(0.0, Reactions[i]->reactionRate)) {
      // cout<<" r"<<i<<" :"<<Reactions[i]->h<<endl;
      if (CaWaveFlag == 0) {
        // sum = sum + Reactions[i]->getReactionRate(Species);
        sum = sum + Reactions[i]->getReactionRate(Species);
      } else {
        sum = sum + Reactions[i]->getReactionRate_wave(Species);
      }
      if (LT(tmp, sum)) {  // tmp<sum
        Miu = i;
        // cout<<"Miu :"<<Miu<<endl;
        break;
      }
    }
  }  // end for loop
}

int SV::getMiu(double uniRand) {
  setMiu_new(uniRand);
  return Miu;
}

void SV::setdMiu(double uniRand) {
  if (numberOfSpecies == 1) {
    dMiu = 0;  // when there is only one specie
    return;
  }
  double Urandom1 = (double)rand() / (double)RAND_MAX;
  double tmp = Urandom1 * s;  //
  double directions = Neighbours.size();

  double sum = 0;
  for (int i = 0; i < numberOfSpecies; i++) {
    // if((double)(Species[i]->diffusionRate)>0.0){
    if (LT(0.0, Species[i]->diffusionRate)) {
      sum = sum + Species[i]->getDiffusionRate(directions, distance,
                                               surfaceArea, volume, CaWaveFlag);
      if (LT(tmp, sum)) {
        dMiu = i;
        // cout<<"SV::setdMiu "<<dMiu<<endl;
        break;
      }
    }  //
  }
}  // end of setdMiu

int SV::getdMiu(double uniRand) {
  setdMiu(uniRand);
  return dMiu;
}

void SV::setCaWaveFlag(int waveFlag) { CaWaveFlag = waveFlag; }

// less than LT
bool SV::LT(long double a, long double b) {
  double epsilon = .0000000001;  // 12 ta precision
  return (b - a) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

long SV::double2long(long double d) {
  long l = static_cast<long>((static_cast<long double>(d * 10000000)));
  // 8 ta 0
  return l;
}

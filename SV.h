/*
 * SV.h
 *
 *  Created on: 2012-11-20
 *      Author: mpatoa
 */

#include <vector>
#include "reaction.h"
#include "species.h"
using namespace std;

#ifndef SV_H_
#define SV_H_

class SV {
 public:
  int numberOfOpen;
  int numberOfEvents;  // any events RD
  int CaWaveFlag;
  void setCaWaveFlag(int waveFlag);

  long double surfaceArea;  // the area of the boundary of the neuron lying
                            // withing the cube
  long double
      volume;  // the volume of the intersection of the cube with the neuron
  long double area;  // the area of the shared face between two adjacent cubes
  long double distance;  // how far apart the centers are
  long double segment;   // mapping to 1D electrophys
  int ID;

  vector<int> Neighbours;  // list of neighbors
  vector<int> Location;    // list of Location[0]=dx, Location[1]=dy
                         // Location[2]=dz

  int numberOfSpecies;
  vector<species *> Species;  /////////////

  int numberOfReactions;  //
  vector<reaction *> Reactions;

  long double pqTau;  //***
  long double Tau;
  long double Lt;  // local time
  long double Gt;
  long double r;  // this is a0
  long double s;

  long double testTau;

  int Miu;
  int dMiu;

  int reactionDirection;

  SV();
  virtual ~SV();
  SV(int id, int nSpecies, int nReactions);

  void updateRates();

  void setTau(double simTime);
  long double getTau();

  void setMiu(double uniRand);
  void setMiu_new(double uniRand);
  int getMiu(double uniRand);

  void setdMiu(double uniRand);
  int getdMiu(double uniRand);

  void setReactionRates();
  long double getReactionRates();

  void setDiffusionRates();
  long double getDiffusionRates();
  bool LT(long double simTime, long double nextSimTime);
  long double2long(long double d);
};

#endif /* SV_H_ */

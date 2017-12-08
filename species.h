/*
 * species.h
 *
 *  Created on: 2013-07-15
 *      Author: mpatoa
 */
#ifdef __cplusplus  // only actually define the class if this is C++
#include <string>
using namespace std;
#ifndef SPECIES_H_
#define SPECIES_H_

class species {
 public:
  string speciesName;
  int CaWaveFlag;
  long numberOfMolecules;
  // double h;
  double diffusionConstant;
  long double diffusionRate;
  int speciesID;
  species();
  virtual ~species();
  species(int id);
  long double getDiffusionRate(double directions, double distance,
                               double surfaceArea, double volume,
                               int caWaveFlag);
  bool LT(long double simTime, long double nextSimTime);
};
#endif /* SPECIES_H_ */
#else  // C doesn't know about classes, just say it's a struct
typedef struct species species;
#endif

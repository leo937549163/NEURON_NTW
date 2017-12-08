#include <math.h>
#include <stdio.h>
#include "species.h"

// rate = 0.01 * pow(getNumberOfMolecules(Species[0]), 1) *
//            getNumberOfMolecules(Species[1]) -
//        getNumberOfMolecules(Species[2]) * 0.01;
double C_getNumberOfMolecules(species *Species);
double C_getNumberOfMoleculesPerVol(species *Species);

double reaction_0(species *Species[]) {
  // double rate = pow(C_getNumberOfMolecules(Species[0]),2) *
  // C_getNumberOfMolecules(Species[1]);
  double rate = C_getNumberOfMoleculesPerVol(Species[0]) *
                C_getNumberOfMoleculesPerVol(Species[1]);
  return rate;
}

double reaction_1(species *Species[]) {
  double rate = C_getNumberOfMoleculesPerVol(Species[2]);
  return rate;
}

// double reaction_2(species *Species[]) {
//   double rate = (0.000106274491  - C_getNumberOfMoleculesPerVol(Species[0])) *
//                 (0.000637646946 - C_getNumberOfMoleculesPerVol(Species[0])) *
//                 (0.00191294084 - C_getNumberOfMoleculesPerVol(Species[0])) *
//                 1000;
//   return rate;
// }

double reaction_2(species *Species[]) {
  double rate = (0.0001  - C_getNumberOfMoleculesPerVol(Species[0])) *
                (0.0006 - C_getNumberOfMoleculesPerVol(Species[0])) *
                (0.002 - C_getNumberOfMoleculesPerVol(Species[0])) *
                1000;
  return rate;
}
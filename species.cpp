/*
 * species.cpp
 *
 *  Created on: 2014-08-19
 *      Author: Mohammad Nazrul Ishlam Patoary
 */
#include "species.h"
#include <math.h>
#include <iostream>
// using namespace std;

species::species() {}
species::~species() {}

species::species(int id) {
  speciesID = id;
  speciesName = "";
  numberOfMolecules = 0;
  diffusionConstant = 0.0;
  CaWaveFlag = 0;
}

long double species::getDiffusionRate(double directions, double distance,
                                      double surfaceArea, double volume,
                                      int caWaveFlag) {
  double n = directions;
  // number of directions in which the molecules can diffuse
  if (caWaveFlag > 0) {
    n = 2;
  }
  // cout << "diffusionConstant: " << diffusionConstant << endl;
  // cout << "surfaceArea: " << surfaceArea << endl;
  // cout << "distance: " << distance << endl;
  // cout << "volume: " << volume << endl;

  long double dj = (long double)((double)(diffusionConstant * surfaceArea) /
                                 (double)(distance * volume));

  diffusionRate = (long double)(n * dj * numberOfMolecules);
  
  return diffusionRate;
}

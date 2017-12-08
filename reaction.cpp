/*
 * reaction.cpp
 *
 *  Created on: 2014-08-19
 *      Author: Mohammad Nazrul Ishlam Patoary
 */
#include "reaction.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "species.h"

reaction::reaction() {
  // TODO Auto-generated constructor stub
}

reaction::~reaction() {
  // TODO Auto-generated destructor stub
}

reaction::reaction(int id) {
  reactionConstant = 0.0;
  reactionID = id;
  sampleCounter = 0;
  FC = .83;  //.83
  FE = .17;
  VIP3R = .0002;  //.0002
  KIP3 = .0013;   //.0013
  KAct = .0004;   //.0004
  kInh = .0019;   //.0019
  tau = 400;
  VLeak = 0.00003;               // 0.00003
  VSERCA = (double)(.00003249);  // Ager value:.00003249 komale say:000003249
                                 // korle pump o kome fole, Ca_C level beshi
                                 // thake.
  KSERCA = .0001;                // .0001
  getReactionRate_C = NULL;
}

double C_getNumberOfMolecules(species *Species) {
  return Species->numberOfMolecules;
}

// converting from mM um^3 to molecules
// = 6.02214129e23 * 1000. / 1.e18 / 1000
// = avogadro * (L / m^3) * (m^3 / um^3) * (mM / M)
// value for avogardro's constant from NIST webpage, accessed 25 April 2012:
// http://physics.nist.gov/cgi-bin/cuu/Value?na
const double CONVERSION_FACTOR = 602214.129;
const double CONVERSION_FACTOR_VOLUME = CONVERSION_FACTOR / 64.0;

double C_getNumberOfMoleculesPerVol(species *Species) {
  return (double)Species->numberOfMolecules / CONVERSION_FACTOR_VOLUME;
}

long double reaction::getReactionRate(vector<species *> Species) {
  long double rc = reactionConstant;

  rc = (long double)(rc * getReactionRate_C(&Species[0]) *
                     CONVERSION_FACTOR_VOLUME);

  reactionRate = rc;  // reaction rate is updated here . . .
  return reactionRate;
}

long double reaction::getReactionRate_wave(vector<species *> Species) {
  double r = 0;  //
  sampleCounter++;

  Ca_C = (double)((Species[0]->numberOfMolecules));
  Ca_E = (double)((Species[1]->numberOfMolecules));
  IP3 = (double)((Species[2]->numberOfMolecules));
  inhibit = (double)((Species[3]->numberOfMolecules));
  not_inhibit = (double)((Species[4]->numberOfMolecules));
  not_active = (double)((Species[5]->numberOfMolecules));
  active = (double)((Species[6]->numberOfMolecules));

  Ca_C = (double)(Ca_C /
                  (6.023 * .25 * .25 * .25 * 100000));  // 2/9411 = 0.000212 mM
  Ca_E =
      (double)(Ca_E / (6.023 * .25 * .25 * .25 * 100000));  // 200/9411 = .02125
  IP3 = (double)(IP3 /
                 (6.023 * .25 * .25 * .25 * 100000));  // middle svs, .02125 mM
  inhibit = (double)(inhibit / (6.023 * .25 * .25 * .25 *
                                100000));  // 1/9411 = 0.000106 mM
  not_inhibit = (double)(not_inhibit / (6.023 * .25 * .25 * .25 * 100000));
  active = (double)(active / (6.023 * .25 * .25 * .25 * 100000));
  not_active = (double)(not_active / (6.023 * .25 * .25 * .25 * 100000));

  hinf = (double)(kInh / (kInh + Ca_C));
  h = (double)(not_inhibit / (inhibit + not_inhibit));
  ninf = (double)(Ca_C / (Ca_C + KAct));

  n = (double)(active / (active + not_active));

  m = (double)(IP3 / (IP3 + KIP3));

  JIP3R = (double)(VIP3R * n * m * h * (Ca_E - Ca_C));

  // JIP3R = (double)(VIP3R*n*m*h*(Ca_E));
  // cout<<JIP3R<<" -------- "<<Ca_E-Ca_C<<endl;

  JSERCA = (double)(VSERCA *
                    ((double)((1.0 * Ca_C * Ca_C) /
                              (double)((KSERCA * KSERCA) + (Ca_C * Ca_C)))));
  // Ca_C upore chilo 22 04 16
  JLeak = (double)(VLeak * (Ca_E - Ca_C));
  // JLeak = (double) (VLeak*(Ca_E));

  if (reactionID == 0) {
    r = reactionConstant * Species[3]->numberOfMolecules * hinf;
    // cout<<"R : "<<reactionID<<"  rate: "<<r<<endl;
    // cout<<reactionID<<"  "<<reactionID<<"  "<<reactionID<<"
    // "<<reactionID<<endl;  r=0;

  } else if (reactionID == 1) {
    // cout<<reactionConstant*Species[4]->numberOfMolecules<<endl;
    r = reactionConstant * Species[4]->numberOfMolecules * (1 - hinf);
    // cout<<"R : "<<reactionID<<"  rate: "<<r<<endl;
    // cout<<reactionID<<"  "<<reactionID<<"  "<<reactionID<<"
    // "<<reactionID<<endl;  r=0;

  } else if (reactionID == 2) {
    r = reactionConstant * Species[5]->numberOfMolecules * (1 - ninf);
    // cout<<"R : "<<reactionID<<"  rate: "<<r<<endl;
    // cout<<reactionID<<"  "<<reactionID<<"  "<<reactionID<<"
    // "<<reactionID<<endl;  r=0;

  } else if (reactionID == 3) {
    r = reactionConstant * Species[6]->numberOfMolecules * (ninf);
    // cout<<"R : "<<reactionID<<"  rate: "<<r<<endl;
    // cout<<reactionID<<"  "<<reactionID<<"  "<<reactionID<<"
    // "<<reactionID<<endl;  r=0;

  }

  else if (reactionID == 4) {
    if (JIP3R < 0) {
      JIP3R = -JIP3R;
      // cout<<"reach here !"<<endl;
    }
    r = reactionConstant * Species[1]->numberOfMolecules * (JIP3R * 9411);
    // cout<<"R : "<<reactionID<<"  rate: "<<r<<endl;
    // cout<<reactionID<<"  "<<reactionID<<"  "<<reactionID<<"
    // "<<reactionID<<"*****"<<endl;

  } else if ((reactionID == 5)) {
    r = reactionConstant * (JSERCA * 9411);
    // cout<<"R : "<<reactionID<<"  rate: "<<r<<endl;
    // cout<<reactionID<<"  "<<reactionID<<"  "<<reactionID<<"
    // "<<reactionID<<"##############"<<endl;

  } else {
    if (JLeak < 0) {
      JLeak = -JLeak;
    }
    r = reactionConstant * (JLeak * 9411);
    // cout<<"R : "<<reactionID<<"  rate: "<<r<<endl;
    // cout<<reactionID<<"  "<<reactionID<<"  "<<reactionID<<"
    // "<<reactionID<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
  }

  reactionRate = r * 650;  // 750

  return reactionRate;
}  // end of funciton

bool reaction::LT(long double a, long double b) {
  double epsilon = .00000001;
  return (b - a) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

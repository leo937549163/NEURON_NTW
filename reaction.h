/*
 * reaction.h
 *
 *  Created on: 2013-07-18
 *      Author: mpatoa
 */

#ifndef REACTION_H_
#define REACTION_H_
#include <string>
#include <vector>
#include "species.h"

using namespace std;
class reaction {
 public:
 public:
  double FC;  //.83
  double FE;
  double VIP3R;  // 8.5
  double KIP3;
  double KAct;
  double Ca_C;
  double Ca_E;
  double IP3;
  double inhibit;      // IPR_C;
  double not_inhibit;  // IPR_O;
  double active;
  double not_active;
  double kInh;
  double tau;
  double hinf;
  double ninf;
  double VLeak;
  double VSERCA;
  double KSERCA;
  double h;
  double n;
  double m;
  double JIP3R;
  double JLeak;
  double JSERCA;

  char ca_c_rate[256];
  char ca_e_rate[256];
  int sampleCounter;
  double reactionConstant;
  int reactionID;
  int reactionType;

  long double reactionRate;
  // long double reactionRatePerVol;

  vector<double> reactantMetricValue;
  vector<int> reactantSpeciesID;
  vector<double> productMetricValue;
  vector<int> productSpeciesID;
  reaction();
  virtual ~reaction();
  reaction(int id);
  reaction(int id, double (*getReactionRate_C)(species *Species[]));
  long double getReactionRate(vector<species *> Species);
  // long double getReactionRatePerVol(vector<species *> Species);

  double (*getReactionRate_C)(species *Species[]);

  long double getReactionRate_wave(vector<species *> Species);
  // long double getReactionRate_ca_r(vector<species *> Species);
  bool LT(long double simTime, long double nextSimTime);
};

#endif /* REACTION_H_ */

#ifdef __cplusplus
extern "C" {
#endif
double C_getNumberOfMolecules(species *Species);
double C_getNumberOfMoleculesPerVol(species *Species);
#ifdef __cplusplus
}
#endif
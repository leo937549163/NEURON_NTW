#include <vector>
#include "reaction.h"
#include "species.h"
using namespace std;

void createReactions(
    vector<double> putreactantMetricValue, vector<int> putreactantSpeciesID,
    vector<double> putproductMetricValue, vector<int> putproductSpeciesID,
    double putreactionConstant, double dx_,
    double (*getReactionRate_C)(species *Species[]));
void mainLoop(double nextTimeStep);  // need to define here !

void initSVs(double dx_);
void getInformationsFromPython(int sv_size, vector<vector<int> > lon,
                               double *st, int numOfSpecies, double *d,
                               int waveFlag);
void statistics();
void mpiInitialization();
void mpiFinalization();
void setRandomSeed(int x);

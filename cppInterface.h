#include "species.h"

#ifdef __cplusplus
class myClass {
 public:
  myClass();
  ~myClass();
  void cpp_function();
};
#endif

#ifdef __cplusplus
extern "C" {
#endif
void *get_myClass();
void Send_to_cppInterface(void *objectPtr, int tt, long (*neighbours)[6],
                          double *states, int numberOfSpecies, double *d,
                          int CaWaveFlag);
void createReaction_ccpInterface(
    int numberOfReactant, double *reactantMetricValue, long *reactantSpeciesID,
    int numberOfProduct, double *productMetricValue, long *productSpeciesID,
    double *reactionConstant, double *dx,
    double (*getReactionRate_C)(species *Species[]));
void initSVs_cppInterface(double *dx);
void runMainLoop_ccpInterface(double *nextTimeStep);
void finalization_cppInterface();
void setRandomSeed_cppInterface(int s);
#ifdef __cplusplus
}
#endif

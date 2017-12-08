/*
 * SV.cpp
 *
 *  Created on: 2014-08-19
 *      Author: Mohammad Nazrul Ishlam Patoary
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "cppInterface.h"
#include "species.h"

long length(long *flat_offsets, long list_id) {
  /* return the length of an item from a flattened list of lists */
  return (flat_offsets[list_id + 1] - flat_offsets[list_id]);
}

long item(long *flat_items, long *flat_offsets, long list_id, long item_id) {
  /* return a selected item from a flattened list of lists */
  return (flat_items[flat_offsets[list_id] + item_id]);
}

// sending: n = number of SVs, states = staeList, numberOfSpecies, d =
// diffusionList. . .
void write_connectivity(int n, long *flat_items, long *flat_offsets,
                        double *states, int numberOfSpecies, double *d,
                        int CaWaveFlag) {
  long list_id, len;
  long neighbours[n][6];
  int i = 0, j = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < 6; j++) {
      neighbours[i][j] = -1;
    }
  }

  for (list_id = 0; list_id < n; list_id++) {
    len = length(flat_offsets, list_id);
    for (i = 0; i < len; i++) {
      neighbours[list_id][i] = (long)item(flat_items, flat_offsets, list_id, i);
    }
  }

  void *temp_obj = get_myClass();
  Send_to_cppInterface(
      temp_obj, n, neighbours, states, numberOfSpecies, d,
      CaWaveFlag);  // datastructure is sending to c++ interface . .
}  // end of write connectivity

void initSVs(double *dx) { initSVs_cppInterface(dx); }

void createReaction(int numberOfReactant, double *reactantMetricValue,
                    long *reactantSpeciesID, int numberOfProduct,
                    double *productMetricValue, long *productSpeciesID,
                    double *reactionConstant, double *dx,
                    double (*getReactionRate_C)(species *Species[])) {
  createReaction_ccpInterface(numberOfReactant, reactantMetricValue,
                              reactantSpeciesID, numberOfProduct,
                              productMetricValue, productSpeciesID,
                              reactionConstant, dx, getReactionRate_C);
}

void runMainLoop(double *nextTimeStep) {
  runMainLoop_ccpInterface(nextTimeStep);
}

void finalization() { finalization_cppInterface(); }

void setRandomSeed(int s) { setRandomSeed_cppInterface(s); }

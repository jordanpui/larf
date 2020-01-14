#pragma once
#include "global.h"

class TdmNet;

void partition(vector<vector<int>> &outputClusters, int k);
void getTdmNets(vector<int> &instClusterIndex, vector<TdmNet *> &tdmNets);
void formPlSubproblem(vector<vector<int>> &outputClusters, vector<TdmNet *> &nets);
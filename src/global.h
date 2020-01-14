#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include <thread>
#include <mutex>
#include <iostream>
#include <fstream>
#include <sstream>

#include <string>
#include <array>
#include <vector>
#include <stack>
#include <queue>
#include <list>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <algorithm>

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <cfloat>
#include <cmath>
#include <cassert>

#include <boost/functional/hash.hpp>

#include "utils/log.h"
#include "utils/misc.h"
#include "utils/geo.h"

using namespace std;

class Setting {
public:
    enum AlgoFlow { Flow_Tdm_Part, Flow_Tdm_Place, Flow_Tdm_Time };

    enum ContMethod { Tdm_ILP, Tdm_LP, Tdm_Lag, Tdm_Iter, Tdm_None };

    enum LgMethod { Lg_Disp, Lg_SR, Lg_None, Lg_MaxDisp };

    string io_out;
    string io_aux;
    string io_nodes;
    string io_nets;
    string io_pl;
    string io_wts;
    string io_scl;
    string io_lib;
    int nPartition;

    AlgoFlow flow;
    ContMethod cont;
    LgMethod lg;
    int nThreads;
    int lagIter;
    bool doRefine;
    bool computeDual;

    Setting() {
        cont = Tdm_Lag;
        lg = Lg_MaxDisp;
        flow = Flow_Tdm_Time;
        nThreads = 8;
        lagIter = 1000;
        doRefine = false;
        computeDual = false;
    }
};

extern Setting setting;

#endif

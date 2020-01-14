#include "db/db.h"
using namespace db;
#include "gp/gp.h"
#include "tdm/tdm_db.h"
#include "tdm/tdm_part.h"
#include "tdm/tdm_solve_lp.h"
#include "tdm/tdm_solve_lag.h"
#include "tdm/tdm_leg.h"
#include "tdm/tdm_refine_lp.h"
#include "tdm/tdm_refine_greedy.h"

Setting setting;

bool get_args(int argc, char **argv);

int main(int argc, char **argv) {
    init_log(LOG_NORMAL);

    log() << "---------------------------------------------------------------------" << endl;
    log() << "  Timing-Division-Multiplexing Optimization for Multi-FPGAs Systems  " << endl;
    log() << "  Author: Chak-Wa Pui  " << endl;
    log() << "---------------------------------------------------------------------" << endl;

    if (!get_args(argc, argv)) return 1;
    database.readAux(setting.io_aux);
    database.readLib(setting.io_lib);
    database.readNodes(setting.io_nodes);
    database.readScl(setting.io_scl);
    database.readPl(setting.io_pl);
    database.readNets(setting.io_nets);

    database.setup();
    database.print();

    gpSetting.init();

    if (setting.flow == Setting::Flow_Tdm_Part) {
        vector<vector<int>> clusters;
        partition(clusters, setting.nPartition);

        vector<int> instClusterIndex(database.instances.size());
        for (unsigned i = 0; i < clusters.size(); i++) {
            for (auto instId : clusters[i]) {
                instClusterIndex[instId] = i;
            }
        }

        vector<TdmNet *> tdmNets;
        getTdmNets(instClusterIndex, tdmNets);
        formPlSubproblem(clusters, tdmNets);
        log() << "finish generating sub problems" << endl;
    } else if (setting.flow == Setting::Flow_Tdm_Time) {
        log() << "------------------------------------------------------" << endl;
        log() << "                begin TDM optimization                " << endl;
        log() << "------------------------------------------------------" << endl;

        vector<Group> groups(database.instances.size());
        for (unsigned int i = 0; i < groups.size(); i++) {
            groups[i].instances.push_back(database.instances[i]);
            groups[i].id = i;
        }
        ifstream instPlFile("instance.pos");
        unordered_map<string, int> nameToIndx;
        for (unsigned int i = 0; i < groups.size(); i++) nameToIndx[groups[i].instances[0]->name] = i;

        for (unsigned int i = 0; i < groups.size(); i++) {
            string instName;
            double x, y;
            instPlFile >> instName;
            instPlFile >> x >> y;
            int idx = nameToIndx[instName];
            groups[idx].x = x;
            groups[idx].y = y;
        }
        instPlFile.close();

        tdmDatabase.init(setting.nPartition, &groups);
        tdmDatabase.getOptXdrVars();

        if (setting.cont == Setting::Tdm_LP) {
            TdmLpSolver tdmLpSolver(true);
            tdmLpSolver.solve();
        } else if (setting.cont == Setting::Tdm_ILP) {
            TdmLpSolver tdmLpSolver(false);
            tdmLpSolver.solve();
        } else if (setting.cont == Setting::Tdm_Lag) {
            TdmLagSolver tdmLagSolver;
            tdmLagSolver.solve();
        } else if (setting.cont == Setting::Tdm_None) {
            tdmDatabase.readSol(database.bmName + "_cont.tdm");
        }
        tdmDatabase.writeSol(database.bmName + "_cont.tdm");

        if (setting.doRefine) {
            TdmRefineLP contRefiner(true);
            contRefiner.solve();
            tdmDatabase.writeSol(database.bmName + "_ref.tdm");
            tdmDatabase.readSol(database.bmName + "_ref.tdm");
        }

        if (setting.lg != Setting::Lg_None) {
            TdmLegalize legalizer;
            legalizer.solve();
        }

        // TdmRefineLP discRefiner(false);
        // discRefiner.solve();

        TdmRefine greedyRefiner;
        greedyRefiner.solve();

        tdmDatabase.reportSol();
        tdmDatabase.writeSol(setting.io_out);

        log() << "finish tdm optimization" << endl;
    } else if (setting.flow == Setting::Flow_Tdm_Place) {
        vector<Group> groups(database.instances.size());
        for (unsigned int i = 0; i < groups.size(); i++) {
            groups[i].instances.push_back(database.instances[i]);
            groups[i].id = i;
        }

        gpSetting.set2();
        gplace(groups);
        gpSetting.set3();
        gplace(groups);

        ofstream fs(setting.io_out);
        for (auto &group : groups) {
            fs << group.instances[0]->name << " " << group.x << " " << group.y << endl;
        }
        fs.close();

        log() << "finish placement for sub problem " << database.bmName << endl;
    }

    log() << "-----------------------------------" << endl;
    log() << "           terminating...          " << endl;
    log() << "-----------------------------------" << endl;
    return 0;
}

bool get_args(int argc, char **argv) {
    bool valid = true;
    for (int a = 1; a < argc; a++) {
        if (strcmp(argv[a], "-out") == 0) {
            setting.io_out.assign(argv[++a]);
            unsigned pos = setting.io_out.find_last_of('.');
            database.bmName = setting.io_out.substr(0, pos);
        } else if (strcmp(argv[a], "-aux") == 0) {
            setting.io_aux.assign(argv[++a]);
        } else if (strcmp(argv[a], "-flow") == 0) {
            string flowname(argv[++a]);
            if (flowname == "tdm_part") {
                setting.flow = Setting::Flow_Tdm_Part;
            } else if (flowname == "tdm_place") {
                setting.flow = Setting::Flow_Tdm_Place;
            } else if (flowname == "tdm_time") {
                setting.flow = Setting::Flow_Tdm_Time;
            } else {
                cerr << "unknown flow: " << flowname << endl;
                valid = false;
            }
        } else if (strcmp(argv[a], "-cont") == 0) {
            string methodname(argv[++a]);
            if (methodname == "ILP") {
                setting.cont = Setting::Tdm_ILP;
            } else if (methodname == "LP") {
                setting.cont = Setting::Tdm_LP;
            } else if (methodname == "Lag") {
                setting.cont = Setting::Tdm_Lag;
            } else if (methodname == "None") {
                setting.cont = Setting::Tdm_None;
            } else {
                cerr << "unknown method: " << methodname << endl;
                valid = false;
            }
        } else if (strcmp(argv[a], "-lg") == 0) {
            string methodname(argv[++a]);
            if (methodname == "Disp") {
                setting.lg = Setting::Lg_Disp;
            } else if (methodname == "None") {
                setting.lg = Setting::Lg_None;
            } else if (methodname == "MaxDisp") {
                setting.lg = Setting::Lg_MaxDisp;
            } else {
                cerr << "unknown method: " << methodname << endl;
                valid = false;
            }
        } else if (strcmp(argv[a], "-partition") == 0) {
            setting.nPartition = atoi(string(argv[++a]).c_str());
        } else if (strcmp(argv[a], "-thread") == 0) {
            setting.nThreads = atoi(string(argv[++a]).c_str());
        } else if (strcmp(argv[a], "-refine") == 0) {
            setting.doRefine = true;
        } else if (strcmp(argv[a], "-lagIter") == 0) {
            setting.lagIter = atoi(string(argv[++a]).c_str());
        } else if (strcmp(argv[a], "-computeDual") == 0) {
            setting.computeDual = true;
        } else {
            cerr << "unknown parameter: " << argv[a] << endl;
            valid = false;
        }
    }
    if (valid && setting.io_aux == "") valid = false;
    if (!valid) log() << "usage: " << argv[0] << " -aux <aux_file> [-out <output_file>]" << endl;
    return valid;
}

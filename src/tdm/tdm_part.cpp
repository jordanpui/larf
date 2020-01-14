#include "tdm_part.h"
#include "alg/patoh/patoh.h"
#include "db/db.h"
using namespace db;

#include "gp/gp.h"
#include "tdm_net.h"

void partition(vector<vector<int>> &clusters, int k) {
    vector<int> inputCluster(database.instances.size());
    for (unsigned i = 0; i < inputCluster.size(); ++i) inputCluster[i] = i;

    if (k == 1) {
        clusters.push_back(inputCluster);
        return;
    }

    clusters.resize(k);

    // Preprocessing
    PaToH_Parameters args;
    // cell/instance
    int c = 0;  // num of cells
    vector<int> inst2cell(database.instances.size(), -1);
    for (auto instId : inputCluster) inst2cell[instId] = c++;
    int *cwghts = new int[c];
    for (int i = 0; i < c; ++i) cwghts[i] = 1;
    // net/pin
    int n = 0;  // num of nets
    vector<bool> isCsdNet(database.nets.size(), false);
    int p = 0;  // num of pins
    for (auto net : database.nets) {
        bool exist = false;
        for (auto pin : net->pins)
            if (inst2cell[pin->instance->id] != -1) {
                ++p;
                exist = true;
            }
        if (exist) {
            ++n;
            isCsdNet[net->id] = true;
        }
    }
    // net to cell
    int *xpins = new int[n + 1];
    xpins[0] = 0;
    int *pins = new int[p];
    int in = 0, ip = 0;
    for (auto net : database.nets)
        if (isCsdNet[net->id]) {
            for (auto pin : net->pins)
                if (inst2cell[pin->instance->id] != -1) pins[ip++] = inst2cell[pin->instance->id];
            xpins[++in] = ip;
        }
    assert(ip == p && in == n);

    // Init
    PaToH_Initialize_Parameters(&args, PATOH_CONPART, PATOH_SUGPARAM_QUALITY);
    args.seed = 0;
    args._k = k;
    args.final_imbal = 0.05;
    int *partvec = new int[c];
    int *partweights = new int[args._k];
    int cut;
    PaToH_Check_User_Parameters(&args, true);
    PaToH_Alloc(&args, c, n, 1, cwghts, NULL, xpins, pins);

    // Processing
    PaToH_Part(&args, c, n, 1, 0, cwghts, NULL, xpins, pins, NULL, partvec, partweights, &cut);

    // Postprocessing
    for (int i = 0; i < c; ++i) clusters[partvec[i]].push_back(inputCluster[i]);

    printlog(LOG_INFO, "cut size: %d, n: %d, c: %d", cut, n, c);

    for (int i = 0; i < k; i++) printlog(LOG_INFO, "%d: weight=%d", i, partweights[i]);

    delete[] cwghts;
    delete[] xpins;
    delete[] pins;
    delete[] partweights;
    delete[] partvec;
    PaToH_Free();
}

void getTdmNets(vector<int> &instDeviceIdx, vector<TdmNet *> &tdmNets) {
    vector<set<int>> netDeviceIdxs(database.nets.size());
    for (unsigned i = 0; i < database.nets.size(); i++) {
        for (auto pin : database.nets[i]->pins) {
            netDeviceIdxs[i].insert(instDeviceIdx[pin->instance->id]);
        }
    }

    map<int, vector<int>> crossDevCnt;
    for (unsigned i = 0; i < database.nets.size(); i++) {
        int devCnt = netDeviceIdxs[i].size();
        crossDevCnt[devCnt].push_back(i);
    }
    log() << "cross device net cnt" << endl;
    log() << "#device\t#net" << endl;
    for (auto &pair : crossDevCnt) log() << pair.first << "\t\t" << pair.second.size() << endl;

    for (unsigned i = 0; i < database.nets.size(); i++) {
        if (database.nets[i]->isClk) continue;

        int devCnt = netDeviceIdxs[i].size();
        if (devCnt == 1) {
            TdmNet *tdmNet = new TdmNet(database.nets[i]);
            tdmNet->setIntraNet(true);
            for (auto pin : database.nets[i]->pins) tdmNet->addPin(pin);

            tdmNet->setFromDevice(*netDeviceIdxs[i].begin());
            tdmNet->setToDevice(*netDeviceIdxs[i].begin());
            tdmNets.push_back(tdmNet);
        } else {
            // decomposite into 2-pin net
            vector<TdmNet *> subnets(devCnt, NULL);

            unordered_map<int, int> deviceIdxToSubnets;
            int cnt = 0;
            for (auto idx : netDeviceIdxs[i]) {
                subnets[cnt] = new TdmNet(database.nets[i]);
                subnets[cnt]->setToDevice(idx);
                deviceIdxToSubnets[idx] = cnt;
                cnt++;
            }

            int driverIdx = -1;
            int driverDeviceIdx = -1;
            for (auto pin : database.nets[i]->pins) {
                int idx = deviceIdxToSubnets[instDeviceIdx[pin->instance->id]];
                if (pin->type->type == 'o') {
                    driverIdx = idx;
                    driverDeviceIdx = instDeviceIdx[pin->instance->id];
                }
                subnets[idx]->addPin(pin);
            }

            Pin *driverPin = subnets[driverIdx]->getPin(subnets[driverIdx]->getDriverIdx());
            for (int j = 0; j < devCnt; j++) {
                subnets[j]->setIntraNet(j == driverIdx);
                if (j != driverIdx) {
                    subnets[j]->addPin(driverPin);
                }
                subnets[j]->setFromDevice(driverDeviceIdx);
            }

            // ignore the intranet if it only contain the driver pin
            for (auto tdmNet : subnets)
                if (tdmNet->getNumPins() != 1) tdmNets.push_back(tdmNet);
        }
    }

    int cnt = 0;
    for (auto net : tdmNets)
        if (net->isInterNet()) cnt++;
    log() << "#tdmNet=" << tdmNets.size() << ", #inter-Net=" << cnt << endl;
}

void formPlSubproblem(vector<vector<int>> &clusters, vector<TdmNet *> &tdmNets) {
    int nClusters = clusters.size();
    string parentFolder = database.bmName;
    string rmCommand = "rm -rf ./";
    string mkdirCommand = "mkdir -p ./";
    string cpCommand = "cp ";
    int code;
    code = system((rmCommand + parentFolder + "/*").c_str());
    code = system((mkdirCommand + parentFolder).c_str());

    vector<vector<TdmNet *>> clusterToNet(nClusters);
    for (auto tdmNet : tdmNets) {
        if (tdmNet->isIntraNet()) {
            clusterToNet[tdmNet->getFromDevice()].push_back(tdmNet);
        }
    }

    int intraNetCnt = 0;
    for (int c = 0; c < nClusters; c++) {
        code = system((mkdirCommand + parentFolder + "/" + to_string(c)).c_str());

        code = system((cpCommand + setting.io_aux + " ./" + parentFolder + "/" + to_string(c)).c_str());
        code = system((cpCommand + setting.io_wts + " ./" + parentFolder + "/" + to_string(c)).c_str());
        code = system((cpCommand + setting.io_lib + " ./" + parentFolder + "/" + to_string(c)).c_str());
        code = system((cpCommand + setting.io_scl + " ./" + parentFolder + "/" + to_string(c)).c_str());

        ofstream plFile("./" + parentFolder + "/" + to_string(c) + "/design.pl");
        ofstream nodeFile("./" + parentFolder + "/" + to_string(c) + "/design.nodes");
        for (auto idx : clusters[c]) {
            nodeFile << database.instances[idx]->name << " "
                     << Master::NameEnum2String(database.instances[idx]->master->name) << endl;
            if (database.instances[idx]->fixed) {
                int x = database.instances[idx]->pack->site->x;
                int y = database.instances[idx]->pack->site->y;
                int slot = database.instances[idx]->slot;
                plFile << database.instances[idx]->name << " " << x << " " << y << " " << slot << " FIXED" << endl;
            }
        }
        plFile.close();
        nodeFile.close();

        ofstream intraNetFile("./" + parentFolder + "/" + to_string(c) + "/design.nets");
        for (auto tdmNet : clusterToNet[c]) {
            intraNetFile << "net intranet_" << to_string(intraNetCnt) << " " << tdmNet->getNumPins() << endl;
            for (auto pin : tdmNet->getPins()) {
                intraNetFile << "\t" << pin->instance->name << " " << pin->type->name << endl;
            }
            intraNetFile << "endnet" << endl;
            intraNetCnt++;
        }

        intraNetFile.close();
    }

    ofstream interNetFile("./" + parentFolder + "/design_inter.nets");
    int interNetCnt = 0;
    for (auto tdmNet : tdmNets) {
        if (tdmNet->isInterNet()) {
            interNetFile << "net internet_" << to_string(interNetCnt) << " " << tdmNet->getNumPins() << endl;
            for (auto pin : tdmNet->getPins()) {
                interNetFile << "\t" << pin->instance->name << " " << pin->type->name << endl;
            }
            interNetFile << "endnet" << endl;
            interNetCnt++;
        }
    }
    interNetFile.close();

    vector<int> instClusterIndex(database.instances.size());
    for (unsigned i = 0; i < clusters.size(); i++) {
        for (auto instId : clusters[i]) {
            instClusterIndex[instId] = i;
        }
    }
    ofstream instDeviceFile("./" + parentFolder + "/instance.device");
    for (unsigned i = 0; i < database.instances.size(); i++) {
        instDeviceFile << database.instances[i]->name << " " << instClusterIndex[i] << endl;
    }
    instDeviceFile.close();
}
#include "tdm_solve_lag.h"
#include "tdm_db.h"
#include "tdm_net.h"
#include "timing_graph.h"
#include "gurobi_c++.h"

void TdmLagMultiplierUpdater::removeAccIssue(Node *node) {
    auto *timingGraph = _tdmLagData._timingGraph;

    double driverSum = 0, fanoutSum = 0;
    if (node == timingGraph->getSink())
        fanoutSum = 1;
    else if (node == timingGraph->getSource())
        return;

    for (auto driver : node->_drivers) driverSum += _tdmLagData.getMuVal(driver);
    for (auto fanout : node->_fanouts) fanoutSum += _tdmLagData.getMuVal(fanout);

    if (driverSum != fanoutSum) {
        double diff = fanoutSum - driverSum;
        double maxDriverMu = -1;
        Edge *maxDriver = NULL;
        for (auto driver : node->_drivers) {
            double val = abs(_tdmLagData.getMuVal(driver));
            if (val > maxDriverMu) {
                maxDriverMu = val;
                maxDriver = driver;
            }
        }
        double &mu = _tdmLagData.getMu(maxDriver);
        mu += diff;
        if (mu < 0) mu = 0;
    }
}

void TdmLagMultiplierUpdater::critFlow(Node *node, double driverSum, double fanoutSum) {
    int driverSize = node->_drivers.size();

    auto sortedDrivers = node->_drivers;
    sort(sortedDrivers.begin(), sortedDrivers.end(), [&](Edge *e1, Edge *e2) {
        return _tdmLagData.getMuVal(e1) < _tdmLagData.getMuVal(e2);
    });

    double diff = fanoutSum - driverSum;
    for (int d = 0; d < driverSize; d++) {
        double &mu = _tdmLagData.getMu(sortedDrivers[d]);
        if (driverSum != 0) {
            mu = mu + diff * (mu / driverSum);
        } else {
            mu = mu + diff / (driverSize - d);
            diff -= diff / (driverSize - d);
            if (mu < 0) {
                diff += mu;
                mu = 0;
            }
        }
    }
}

void TdmLagMultiplierUpdater::decreaseFlow(
    double driverSum, double fanoutSum, vector<pair<Edge *, double>> &gradients, int startIdx, Node *node) {
    int driverSize = gradients.size();

    double diff = fanoutSum - driverSum;

    double muSum = 0;
    for (int d = startIdx; d < driverSize; d++) muSum += _tdmLagData.getMuVal(gradients[d].first);
    if (muSum < driverSum - fanoutSum) {
        printlog(LOG_ERROR, "muSum is smaller than diff");
        cout << node->_id << endl;
        cout << "driver: ";
        for (auto driver : node->_drivers) cout << _tdmLagData.getMuVal(driver) << "," << driver->_driver->_id << " ";
        cout << endl;
        cout << "fanout: ";
        for (auto fanout : node->_fanouts) cout << _tdmLagData.getMuVal(fanout) << "," << fanout->_fanout->_id << " ";
        cout << endl;
        getchar();
        return;
    }

    sort(gradients.begin() + startIdx, gradients.end(), [&](pair<Edge *, double> p1, pair<Edge *, double> p2) {
        return p1.second < p2.second;
    });
    double maxAbsGradient = abs(gradients[startIdx].second);
    int lastIdx = driverSize;
    for (int d = driverSize - 1; d >= startIdx; d--) {
        if (_tdmLagData.getMuVal(gradients[d].first) * gradients[d].second != 0)
            break;
        else
            lastIdx = d;
    }

    double sum = 0;
    for (int d = startIdx; d < driverSize; d++)
        sum += _tdmLagData.getMuVal(gradients[d].first) * (abs(gradients[d].second) / maxAbsGradient);

    if (lastIdx == startIdx) {
        for (int d = startIdx; d < driverSize; d++) {
            double &mu = _tdmLagData.getMu(gradients[d].first);

            mu += diff * (mu / muSum);
            mu = max(mu, 0.0);
        }
        return;
    }

    double ratio = abs(diff) / sum;
    int curIdx = startIdx;
    while (ratio > 1) {
        double &mu = _tdmLagData.getMu(gradients[curIdx].first);
        diff += mu;
        sum -= mu;
        muSum -= mu;
        mu = 0;
        curIdx++;

        if (muSum == 0) return;

        if (curIdx == lastIdx) {
            for (int d = curIdx; d < driverSize; d++) {
                double &mu = _tdmLagData.getMu(gradients[d].first);
                mu += diff * (mu / muSum);
                mu = max(mu, 0.0);
            }
            return;
        }

        sum *= maxAbsGradient;
        maxAbsGradient = abs(gradients[curIdx].second);
        sum /= maxAbsGradient;

        ratio = abs(diff) / sum;
    }

    for (int d = curIdx; d < driverSize; d++) {
        double &mu = _tdmLagData.getMu(gradients[d].first);
        mu += mu * ratio * (gradients[d].second / maxAbsGradient);
        mu = max(mu, 0.0);
    }
}

void TdmLagMultiplierUpdater::increaseFlow(Node *node, double driverSum, double fanoutSum) {
    double diff = fanoutSum - driverSum;
    double critSum = 0;
    int numCritMu = 0;

    sort(node->_drivers.begin(), node->_drivers.end(), [&](Edge *edge1, Edge *edge2) {
        return edge1->getArrivalTimeAlongEdge() > edge2->getArrivalTimeAlongEdge();
    });

    const double maxDiffRatio = 0.05;
    double threshold = (*node->_drivers.begin())->getArrivalTimeAlongEdge() * (1 - maxDiffRatio);

    for (auto driver : node->_drivers) {
        if (driver->getArrivalTimeAlongEdge() >= threshold) {
            critSum += _tdmLagData.getMuVal(driver);
            numCritMu++;
        } else {
            break;
        }
    }

    for (auto driver : node->_drivers) {
        if (driver->getArrivalTimeAlongEdge() >= threshold) {
            double &mu = _tdmLagData.getMu(driver);
            if (critSum != 0) {
                mu = mu + diff * (mu / critSum);
            } else {
                mu = mu + diff / numCritMu;
            }
        } else {
            break;
        }
    }
}

void TdmLagMultiplierUpdater::sinkFlow(double driverSum, double fanoutSum, vector<pair<Edge *, double>> &gradients) {
    const double maxDiffRatio = 0.05;
    const double maxNumRatio = 0.01;
    const double maxMuIncr2SumRatio = 0.002;

    double gradientSum = 0;
    for (auto &pair : gradients) gradientSum += pair.second;

    // set stepsize
    sort(gradients.begin(), gradients.end(), [&](pair<Edge *, double> p1, pair<Edge *, double> p2) {
        Edge *e1 = p1.first;
        Edge *e2 = p2.first;
        return e1->getArrivalTimeAlongEdge() > e2->getArrivalTimeAlongEdge();
    });
    Edge *critEdge = gradients.front().first;
    double threshold = critEdge->getArrivalTimeAlongEdge() * (1 - maxDiffRatio);
    int maxNum = max(1.0, gradients.size() * maxNumRatio);
    int lastIncrIdx = -1;
    double curMuSum = 0;
    double curDeltaSum = 0;
    for (int d = 0; d < maxNum && gradients[d].first->getArrivalTimeAlongEdge() >= threshold; d++) {
        double mu = _tdmLagData.getMuVal(gradients[d].first);

        if (curDeltaSum + mu * _ratio > (driverSum - curMuSum - mu) * maxMuIncr2SumRatio) {
            if (d == 0) {
                _ratio = (driverSum - mu) * maxMuIncr2SumRatio / mu;
            } else {
                break;
            }
        }

        lastIncrIdx = d;
        curMuSum += mu;
        curDeltaSum += mu * _ratio;
    }

    // update
    int nIncr = 0;
    for (int d = 0; d <= lastIncrIdx; d++) {
        Edge *e = gradients[d].first;
        double &mu = _tdmLagData.getMu(e);
        double delta = mu * _ratio;

        if (delta > 0) nIncr++;

        mu = mu + delta;
        driverSum += delta;
        gradientSum -= gradients[d].second;
    }

    decreaseFlow(driverSum, fanoutSum, gradients, lastIncrIdx + 1);
}

void TdmLagMultiplierUpdater::updateMu() {
    auto *timingGraph = _tdmLagData._timingGraph;

    vector<vector<Node *>> &revLevels = timingGraph->getRevLevels();
    for (int l = 0, sz = revLevels.size(); l < sz; l++) {
        vector<Node *> &level = revLevels[l];

        int vIndex = 0;
        std::mutex idxMutex;
        std::mutex counterMutex;

        const int batchSize = 20;

        auto propagateFlow = [&]() {
            while (true) {
                idxMutex.lock();
                int idx = vIndex;
                vIndex += batchSize;
                idxMutex.unlock();

                if (idx >= (int)level.size()) break;

                for (int i = idx; i < min((int)level.size(), idx + batchSize); i++) {
                    Node *node = level[i];

                    int driverSize = node->_drivers.size();
                    vector<pair<Edge *, double>> gradients;
                    for (int d = 0; d < driverSize; d++)
                        gradients.emplace_back(node->_drivers[d], _tdmLagData.getMuGrad(node->_drivers[d]));

                    int numCritMu = 0;
                    for (auto &driver : node->_drivers) numCritMu += driver->isCritical();

                    double fanoutSum = 0;
                    if (node == timingGraph->getSink()) {
                        fanoutSum = 1;
                    } else {
                        for (auto fanout : node->_fanouts) fanoutSum += _tdmLagData.getMuVal(fanout);
                    }

                    double driverSum = 0;
                    for (auto driver : node->_drivers) driverSum += _tdmLagData.getMuVal(driver);

                    if (node == timingGraph->getSink()) {
                        if (numCritMu == driverSize) continue;
                        sinkFlow(driverSum, fanoutSum, gradients);
                    } else if (numCritMu == driverSize) {
                        critFlow(node, driverSum, fanoutSum);
                    } else if (driverSum > fanoutSum) {
                        decreaseFlow(driverSum, fanoutSum, gradients, 0, node);
                    } else if (driverSum <= fanoutSum) {
                        increaseFlow(node, driverSum, fanoutSum);
                    } else {
                        printlog(LOG_ERROR,
                                 "uncatch case in updateMu, %d: %d %f %f %d %d",
                                 l,
                                 node->_id,
                                 driverSum,
                                 fanoutSum,
                                 numCritMu,
                                 driverSize);
                        getchar();
                    }

                    removeAccIssue(node);
                }
            }
        };

        int nThreads = max(1, setting.nThreads);
        std::thread threads[nThreads];
        for (int i = 0; i < nThreads; i++) threads[i] = std::thread(propagateFlow);
        for (int i = 0; i < nThreads; i++) threads[i].join();
    }
}

void TdmLagMultiplierUpdater::updateLambda() {
    // update lambda such that the xdr use up all the resources
    for (auto &pair : _tdmLagData._tronconToXdrVar) {
        vector<double> muVec;
        Troncon *troncon = pair.first;
        double sum = 0;
        for (auto var : pair.second) {
            double tmpSum = 0;
            for (auto edge : _tdmLagData._timingGraph->getEdges(_tdmLagData._optXdrVars[var]))
                tmpSum += edge->_tdmCoef * _tdmLagData.getMuVal(edge);

            muVec.push_back(sqrt(tmpSum));
            sum += muVec.back();
        }

        sort(muVec.begin(), muVec.end());
        if (muVec.back() / muVec.front() > _tdmLagData._maxChoice) {
            // int sz = muVec.size();
            // printlog(LOG_ERROR,
            //          "sz=%lu, lambda/u=%f/%f %f %f %f %f %f %f %f %f %f",
            //          muVec.size(),
            //          _tdmLagData.getLambdaVal(troncon),
            //          muVec[0],
            //          muVec[1],
            //          muVec[sz * 0.01],
            //          muVec[sz * 0.1],
            //          muVec[sz * 0.2],
            //          muVec[sz * 0.4],
            //          muVec[sz * 0.6],
            //          muVec[sz * 0.8],
            //          muVec[sz * 0.9],
            //          muVec[sz - 1]);
            break;
        }

        double &lambda = _tdmLagData.getLambda(troncon);
        lambda = pow(sum / troncon->_limit, 2);
        // cout << troncon->_id << " " << lambda << " " << muVec.back() << endl;
        lambda = max(lambda, muVec.back());
    }
}

void TdmLagMultiplierUpdater::getRatio(int iter) { _ratio = _baseRate * pow(0.5, _changeRate * iter); }

double TdmLagMultiplierUpdater::getMuStepSize() {
    auto *timingGraph = _tdmLagData._timingGraph;

    double stepSize = DBL_MAX;
    double maxStepSize = DBL_MIN, minStepSize = DBL_MAX;

    vector<double> stepSizes;

    for (int e = 0, sz = timingGraph->getNumEdges(); e < sz; e++) {
        Edge *edge = timingGraph->getEdge(e);
        double mu = _tdmLagData.getMuVal(edge);

        double gradient = _tdmLagData.getMuGrad(edge);

        if (abs(gradient) > _tdmLagData._epsilon && mu > 0) {
            stepSizes.push_back(_ratio * mu / abs(gradient));
            minStepSize = min(minStepSize, stepSizes.back());
            maxStepSize = max(maxStepSize, stepSizes.back());
        }
    }

    stepSize = (maxStepSize + minStepSize) / 2;

    printlog(LOG_INFO, "mu: min=%f, max=%f, ratio=%.3f, step=%f", minStepSize, maxStepSize, _ratio, stepSize);

    return stepSize;
}

double TdmLagMultiplierUpdater::getLambdaStepSize() {
    double stepSize = DBL_MAX;
    double maxPosStepSize = DBL_MIN, minPosStepSize = DBL_MAX;
    double maxNegStepSize = DBL_MIN, minNegStepSize = DBL_MAX;
    bool hasNeg = false, hasPos = false;

    for (auto &pair : _tdmLagData._tronconToXdrVar) {
        Troncon *troncon = pair.first;
        double lambda = _tdmLagData.getLambdaVal(troncon);
        double gradient = _tdmLagData.getLambdaGrad(troncon);

        if (abs(gradient) > _tdmLagData._epsilon) {
            if (gradient > 0) {
                minPosStepSize = min(minPosStepSize, _ratio * lambda / abs(gradient));
                maxPosStepSize = max(maxPosStepSize, _ratio * lambda / abs(gradient));
                hasPos = true;
            } else if (gradient < 0) {
                minNegStepSize = min(minNegStepSize, _ratio * lambda / abs(gradient));
                maxNegStepSize = max(maxNegStepSize, _ratio * lambda / abs(gradient));
                hasNeg = true;
            }
        }
    }

    if (hasPos && hasNeg) {
        stepSize = min(minNegStepSize, minPosStepSize);
        printlog(LOG_INFO,
                 "lambda(pos/neg): min=%f/%f, max=%f/%f, ratio=%.3f, step=%f",
                 minPosStepSize,
                 minNegStepSize,
                 maxPosStepSize,
                 maxNegStepSize,
                 _ratio,
                 stepSize);
    } else if (hasPos) {
        stepSize = minPosStepSize;
        printlog(LOG_INFO,
                 "lambda(pos/neg): min=%f/na, max=%f/na, ratio=%.3f, step=%f",
                 minPosStepSize,
                 maxPosStepSize,
                 _ratio,
                 stepSize);
    } else if (hasNeg) {
        stepSize = minNegStepSize;
        printlog(LOG_INFO,
                 "lambda(pos/neg): min=na/%f, max=na/%f, ratio=%.3f, step=%f",
                 minNegStepSize,
                 maxNegStepSize,
                 _ratio,
                 stepSize);
    } else {
        stepSize = 0;
        printlog(LOG_INFO, "lambda: ratio=%.3f, step=%f", _ratio, stepSize);
    }

    return stepSize;
}

TdmLagMultiplierUpdater::TdmLagMultiplierUpdater(TdmLagData &tdmLagData) : _tdmLagData(tdmLagData) {
    _preserve.assign(_tdmLagData._timingGraph->getNumEdges(), true);
}

void TdmLagMultiplierUpdater::run(int iter) {
    getRatio(iter);

    updateMu();
    updateLambda();
}

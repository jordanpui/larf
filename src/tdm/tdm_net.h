#pragma once

#include "global.h"

namespace db {
class Pin;
class Net;
}  // namespace db
class XdrVar;
class Troncon;

class TdmNet {
public:
    TdmNet(db::Net* parentNet);

    void setIntraNet(bool flag);
    void setInterNet(bool flag);
    bool isInterNet();
    bool isIntraNet();

    void setFromDevice(int deviceIdx) { _fromDeviceIdx = deviceIdx; }
    void setToDevice(int deviceIdx) { _toDeviceIdx = deviceIdx; }
    int getFromDevice() { return _fromDeviceIdx; }
    int getToDevice() { return _toDeviceIdx; }
    int getNumDevices() { return _fromDeviceIdx == _toDeviceIdx ? 1 : 0; }

    void addPins(vector<db::Pin*>& pins);
    void addPin(db::Pin* pin);
    vector<db::Pin*>& getPins();
    db::Pin* getPin(int i);
    int getDriverIdx();
    int getNumPins();

    string getName();
    XdrVar* getXdrVar() { return _xdrVar; }

    XdrVar* _xdrVar;
    Troncon* _troncon;

private:
    db::Net* _parentNet;
    vector<db::Pin*> _pins;
    bool _intraNet;
    bool _interNet;
    int _driverIdx;
    int _fromDeviceIdx;
    int _toDeviceIdx;
};
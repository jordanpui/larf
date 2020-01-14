#include "tdm_net.h"
#include "db/net.h"

TdmNet::TdmNet(db::Net* parentNet) {
    _parentNet = parentNet;
    _intraNet = false;
    _interNet = false;
    _driverIdx = -1;
    _xdrVar = NULL;
    _troncon = NULL;
}

void TdmNet::setIntraNet(bool flag) {
    _intraNet = flag;
    _interNet = !_intraNet;
}
void TdmNet::setInterNet(bool flag) {
    _interNet = flag;
    _intraNet = !_interNet;
}
bool TdmNet::isInterNet() { return _interNet; }

bool TdmNet::isIntraNet() { return _intraNet; }

void TdmNet::addPins(vector<db::Pin*>& pins) { _pins = pins; }

void TdmNet::addPin(db::Pin* pin) {
    if (pin->type->type == 'o') {
        if (_driverIdx != -1) cout << "error when adding pins" << endl;
        _driverIdx = _pins.size();
    }
    _pins.push_back(pin);
}
vector<db::Pin*>& TdmNet::getPins() { return _pins; }

db::Pin* TdmNet::getPin(int i) { return _pins[i]; }

int TdmNet::getDriverIdx() { return _driverIdx; }

int TdmNet::getNumPins() { return _pins.size(); }

string TdmNet::getName() { return _parentNet->name; }

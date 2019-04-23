#ifndef HPRACOUSTICTENSOR_H
#define HPRACOUSTICTENSOR_H

#include "hprwavespeeds.h"
using namespace std;

class HPRAcousticTensor
{
public:
    HPRAcousticTensor();

    static vector<vector<double> > computeAcousticTensorComponent1(HPRStateVector stateVector, HPRMaterialParameters materialParameters, int direction);
    static vector<vector<double> > computeAcousticTensorComponent2(HPRStateVector stateVector, HPRMaterialParameters materialParameters, int direction);

    static vector<vector<double> > computeAcousticTensor(HPRStateVector stateVector, HPRMaterialParameters materialParameters, int direction);
    static double computeMaximumWaveSpeed(HPRStateVector stateVector, HPRMaterialParameters materialParameters, int direction);
};

#endif // HPRACOUSTICTENSOR_H

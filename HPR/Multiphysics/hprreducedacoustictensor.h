#ifndef HPRREDUCEDACOUSTICTENSOR_H
#define HPRREDUCEDACOUSTICTENSOR_H

#include "HPR/hprwavespeeds.h"
#include "hprreducedstatevector.h"
using namespace std;

class HPRReducedAcousticTensor
{
public:
    HPRReducedAcousticTensor();

    static vector<vector<double> > computeMaterial1AcousticTensorComponent1(HPRReducedStateVector stateVector, HPRMaterialParameters material1Parameters, int direction);
    static vector<vector<double> > computeMaterial2AcousticTensorComponent1(HPRReducedStateVector stateVector, HPRMaterialParameters material2Parameters, int direction);

    static vector<vector<double> > computeMaterial1AcousticTensorComponent2(HPRReducedStateVector stateVector, HPRMaterialParameters material1Parameters, int direction);
    static vector<vector<double> > computeMaterial2AcousticTensorComponent2(HPRReducedStateVector stateVector, HPRMaterialParameters material2Parameters, int direction);

    static vector<vector<double> > computeMaterial1AcousticTensor(HPRReducedStateVector stateVector, HPRMaterialParameters material1Parameters, int direction);
    static vector<vector<double> > computeMaterial2AcousticTensor(HPRReducedStateVector stateVector, HPRMaterialParameters material2Parameters, int direction);

    static double computeMaximumWaveSpeed(HPRReducedStateVector stateVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters, int direction);
};

#endif // HPRREDUCEDACOUSTICTENSOR_H

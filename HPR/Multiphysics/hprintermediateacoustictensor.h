#ifndef HPRINTERMEDIATEACOUSTICTENSOR_H
#define HPRINTERMEDIATEACOUSTICTENSOR_H

#include "HPR/hprwavespeeds.h"
#include "hprintermediatestatevector.h"
using namespace std;

class HPRIntermediateAcousticTensor
{
public:
    HPRIntermediateAcousticTensor();

    static vector<vector<double> > computeMaterial1AcousticTensorComponent1(HPRIntermediateStateVector stateVector, HPRMaterialParameters material1Parameters, int direction);
    static vector<vector<double> > computeMaterial2AcousticTensorComponent1(HPRIntermediateStateVector stateVector, HPRMaterialParameters material2Parameters, int direction);

    static vector<vector<double> > computeMaterial1AcousticTensorComponent2(HPRIntermediateStateVector stateVector, HPRMaterialParameters material1Parameters, int direction);
    static vector<vector<double> > computeMaterial2AcousticTensorComponent2(HPRIntermediateStateVector stateVector, HPRMaterialParameters material2Parameters, int direction);

    static vector<vector<double> > computeMaterial1AcousticTensor(HPRIntermediateStateVector stateVector, HPRMaterialParameters material1Parameters, int direction);
    static vector<vector<double> > computeMaterial2AcousticTensor(HPRIntermediateStateVector stateVector, HPRMaterialParameters material2Parameters, int direction);

    static double computeMaximumWaveSpeed(HPRIntermediateStateVector stateVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters, int direction);
};

#endif // HPRINTERMEDIATEACOUSTICTENSOR_H

#ifndef HPRMULTIPHYSICSACOUSTICTENSOR_H
#define HPRMULTIPHYSICSACOUSTICTENSOR_H

#include "HPR/hprwavespeeds.h"
#include "hprmultiphysicsstatevector.h"
using namespace std;

class HPRMultiphysicsAcousticTensor
{
public:
    HPRMultiphysicsAcousticTensor();

    static vector<vector<double> > computeMaterial1AcousticTensorComponent1(HPRMultiphysicsStateVector stateVector, HPRMaterialParameters material1Parameters, int direction);
    static vector<vector<double> > computeMaterial2AcousticTensorComponent1(HPRMultiphysicsStateVector stateVector, HPRMaterialParameters material2Parameters, int direction);

    static vector<vector<double> > computeMaterial1AcousticTensorComponent2(HPRMultiphysicsStateVector stateVector, HPRMaterialParameters material1Parameters, int direction);
    static vector<vector<double> > computeMaterial2AcousticTensorComponent2(HPRMultiphysicsStateVector stateVector, HPRMaterialParameters material2Parameters, int direction);

    static vector<vector<double> > computeMaterial1AcousticTensor(HPRMultiphysicsStateVector stateVector, HPRMaterialParameters material1Parameters, int direction);
    static vector<vector<double> > computeMaterial2AcousticTensor(HPRMultiphysicsStateVector stateVector, HPRMaterialParameters material2Parameters, int direction);

    static double computeMaximumWaveSpeed(HPRMultiphysicsStateVector stateVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters, int direction);
};

#endif // HPRMULTIPHYSICSACOUSTICTENSOR_H

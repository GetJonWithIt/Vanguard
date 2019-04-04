#ifndef ELASTICACOUSTICTENSOR_H
#define ELASTICACOUSTICTENSOR_H

#include "elasticstatevector.h"
using namespace std;

class ElasticAcousticTensor
{
public:
    ElasticAcousticTensor();

    static vector<vector<double> > computeStressTensorCentredDifference(vector<vector<double> > distortionTensor, double entropy, HyperelasticMaterialParameters materialParameters,
                                                                         double epsilon, int direction);

    static vector<vector<double> > computeAcousticTensor(vector<vector<double> > distortionTensor, double entropy, HyperelasticMaterialParameters materialParameters, double epsilon,
                                                         int direction);
    static double computeMaximumWaveSpeed(vector<vector<double> > distortionTensor, double entropy, HyperelasticMaterialParameters materialParameters, int direction);
};

#endif // ELASTICACOUSTICTENSOR_H

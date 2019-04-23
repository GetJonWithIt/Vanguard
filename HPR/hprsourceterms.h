#ifndef HPRSOURCETERMS_H
#define HPRSOURCETERMS_H

#include "hprequationofstate.h"
using namespace std;

class HPRSourceTerms
{
public:
    HPRSourceTerms();

    static double computeTheta1Reciprocal(double density, vector<vector<double> > distortionTensor, HPRMaterialParameters materialParameters);
    static double computeTheta2Reciprocal(double density, double temperature, HPRMaterialParameters materialParameters);
};

#endif // HPRSOURCETERMS_H

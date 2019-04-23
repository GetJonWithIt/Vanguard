#ifndef HPRDERIVATIVES_H
#define HPRDERIVATIVES_H

#include "hprequationofstate.h"
using namespace std;

class HPRDerivatives
{
public:
    HPRDerivatives();

    static double computeTotalEnergyDerivativeDensity(double density, double pressure, HPRMaterialParameters materialParameters);
    static double computeTotalEnergyDerivativePressure(double density, HPRMaterialParameters materialParameters);

    static vector<vector<double> > computeTotalEnergyDerivativeDistortionTensor(vector<vector<double> > distortionTensor, HPRMaterialParameters materialParameters);
    static vector<double> computeTotalEnergyDerivativeThermalImpulse(double xThermalImpulse, double yThermalImpulse, double zThermalImpulse, HPRMaterialParameters materialParameters);

    static double computeTemperatureDerivativeDensity(double density, double pressure, HPRMaterialParameters materialParameters);
    static double computeTemperatureDerivativePressure(double density, HPRMaterialParameters materialParameters);
};

#endif // HPRDERIVATIVES_H

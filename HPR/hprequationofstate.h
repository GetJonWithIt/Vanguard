#ifndef HPREQUATIONOFSTATE_H
#define HPREQUATIONOFSTATE_H

#include "Mathematics/tensoralgebra.h"
#include "hprmiegruneisen.h"
#include "hprderivatives.h"
using namespace std;

class HPREquationOfState
{
public:
    HPREquationOfState();

    static double computeMicroscaleEnergy(double density, double pressure, HPRMaterialParameters materialParameters);

    static double computeMesoscaleDistortionEnergy(vector<vector<double> > distortionTensor, HPRMaterialParameters materialParameters);
    static double computeMesoscaleThermalEnergy(double xThermalImpulse, double yThermalImpulse, double zThermalImpulse, HPRMaterialParameters materialParameters);

    static double computeMacroscaleEnergy(double xVelocity, double yVelocity, double zVelocity);

    static double computeTotalEnergy(double density, double pressure, double xVelocity, double yVelocity, double zVelocity, vector<vector<double> > distortionTensor, double xThermalImpulse,
                                     double yThermalImpulse, double zThermalImpulse, HPRMaterialParameters materialParameters);

    static double computePressure(double density, double totalEnergy, double xVelocity, double yVelocity, double zVelocity, vector<vector<double> > distortionTensor, double xThermalImpulse,
                                  double yThermalImpulse, double zThermalImpulse, HPRMaterialParameters materialParameters);

    static double computeTemperature(double density, double pressure, HPRMaterialParameters materialParameters);
    static vector<double> computeHeatFluxVector(double temperature, double xThermalImpulse, double yThermalImpulse, double zThermalImpulse, HPRMaterialParameters materialParameters);

    static vector<vector<double> > computeShearStressTensor(double density, vector<vector<double> > distortionTensor, HPRMaterialParameters materialParameters);

    static vector<vector<vector<vector<double> > > > computeShearStressTensorDerivativeDistortionTensor(double density, vector<vector<double> > distortionTensor,
                                                                                                        HPRMaterialParameters materialParameters);
    static vector<vector<double> > computeShearStressTensorDerivativeDensity(vector<vector<double> > distortionTensor, HPRMaterialParameters materialParameters);
};

#endif // HPREQUATIONOFSTATE_H

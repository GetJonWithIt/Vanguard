#ifndef HPRSTATEVECTOR_H
#define HPRSTATEVECTOR_H

#include "hprequationofstate.h"
#include "hprsourceterms.h"
#include "Elasticity/elasticequationofstate.h"
using namespace std;

class HPRStateVector
{
public:
    HPRStateVector();
    HPRStateVector(double newDensity, double newPressure, vector<vector<double> > newDistortionTensor, double newXVelocity, double newYVelocity, double newZVelocity, double newXThermalImpulse,
                   double newYThermalImpulse, double zThermalImpulse);
    HPRStateVector(double newXVelocity, double newYVelocity, double newZVelocity, vector<vector<double> > newDistortionTensor, double newEntropy,
                   HyperelasticMaterialParameters hyperelasticMaterialParameters, HPRMaterialParameters materialParameters);

    void setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector, HPRMaterialParameters materialParameters);
    void setConservedVariableVector(vector<double> newConservedVariableVector, HPRMaterialParameters materialParameters);

    vector<double> computePrimitiveVariableVector(HPRMaterialParameters materialParameters);
    vector<double> computeConservedVariableVector(HPRMaterialParameters materialParameters);

    static vector<double> computeXFluxVector(vector<double> conservedVariableVector, HPRMaterialParameters materialParameters);
    vector<double> computeXFluxVector(HPRMaterialParameters materialParameters);

    static vector<double> computeYFluxVector(vector<double> conservedVariableVector, HPRMaterialParameters materialParameters);
    vector<double> computeYFluxVector(HPRMaterialParameters materialParameters);

    static vector<double> computeSourceTermVector(vector<double> conservedVariableVector, HPRMaterialParameters materialParameters);
    vector<double> computeSourceTermVector(HPRMaterialParameters materialParameters);

    double computeTotalEnergy(HPRMaterialParameters materialParameters);
    double computeTemperature(HPRMaterialParameters materialParameters);
    vector<double> computeHeatFluxVector(HPRMaterialParameters materialParameters);

    double computeTotalEnergyDerivativeDensity(HPRMaterialParameters materialParameters);
    double computeTotalEnergyDerivativePressure(HPRMaterialParameters materialParameters);
    vector<vector<double> > computeTotalEnergyDerivativeDistortionTensor(HPRMaterialParameters materialParameters);
    vector<double> computeTotalEnergyDerivativeThermalImpulse(HPRMaterialParameters materialParameters);

    double computeTemperatureDerivativeDensity(HPRMaterialParameters materialParameters);
    double computeTemperatureDerivativePressure(HPRMaterialParameters materialParameters);

    double computeTheta1Reciprocal(HPRMaterialParameters materialParameters);
    double computeTheta2Reciprocal(HPRMaterialParameters materialParameters);

    vector<vector<double> > computeShearStressTensor(HPRMaterialParameters materialParameters);
    vector<vector<double> > computeShearStressTensorDerivativeDensity(HPRMaterialParameters materialParameters);
    vector<vector<vector<vector<double> > > > computeShearStressTensorDerivativeDistortionTensor(HPRMaterialParameters materialParameters);

    void setDensity(double newDensity);
    void setPressure(double newPressure);
    void setDistortionTensor(vector<vector<double> > distortionTensor);

    void setXVelocity(double newXVelocity);
    void setYVelocity(double newYVelocity);
    void setZVelocity(double newZVelocity);

    void setXThermalImpulse(double newXThermalImpulse);
    void setYThermalImpulse(double newYThermalImpulse);
    void setZThermalImpulse(double newZThermalImpulse);

    double getDensity();
    double getPressure();
    vector<vector<double> > getDistortionTensor();

    double getXVelocity();
    double getYVelocity();
    double getZVelocity();

    double getXThermalImpulse();
    double getYThermalImpulse();
    double getZThermalImpulse();

private:
    double density;
    double pressure;
    vector<vector<double> > distortionTensor;

    double xVelocity;
    double yVelocity;
    double zVelocity;

    double xThermalImpulse;
    double yThermalImpulse;
    double zThermalImpulse;
};

#endif // HPRSTATEVECTOR_H

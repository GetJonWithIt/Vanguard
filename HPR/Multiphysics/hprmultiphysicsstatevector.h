#ifndef HPRMULTIPHYSICSSTATEVECTOR_H
#define HPRMULTIPHYSICSSTATEVECTOR_H

#include "HPR/hprequationofstate.h"
#include "HPR/hprsourceterms.h"
#include "Elasticity/elasticequationofstate.h"
using namespace std;

class HPRMultiphysicsStateVector
{
public:
    HPRMultiphysicsStateVector();
    HPRMultiphysicsStateVector(double newMaterial1VolumeFraction, double newInterfaceXVelocity, double newInterfaceYVelocity, double newInterfaceZVelocity, double newMaterial1Density,
                               double newMaterial1Pressure, vector<vector<double> > newMaterial1DistortionTensor, double newMaterial1XThermalImpulse, double newMaterial1YThermalImpulse,
                               double newMaterial1ZThermalImpulse, double newMaterial2Density, double newMaterial2Pressure, vector<vector<double> > newMaterial2DistortionTensor,
                               double newMaterial2XThermalImpulse, double newMaterial2YThermalImpulse, double newMaterial2ZThermalImpulse);
    HPRMultiphysicsStateVector(double newMaterial1VolumeFraction, double newInterfaceXVelocity, double newInterfaceYVelocity, double newInterfaceZVelocity,
                               vector<vector<double> > newMaterial1DistortionTensor, double newMaterial1Entropy, vector<vector<double> > newMaterial2DistortionTensor, double newMaterial2Entropy,
                               HyperelasticMaterialParameters hyperelasticMaterial1Parameters, HyperelasticMaterialParameters hyperelasticMaterial2Parameters,
                               HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    void setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    void setConservedVariableVector(vector<double> newConservedVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    vector<double> computePrimitiveVariableVector(HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    vector<double> computeConservedVariableVector(HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static vector<double> computeXFluxVector(vector<double> conservedVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    vector<double> computeXFluxVector(HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static vector<double> computeSourceTermVector(vector<double> conservedVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    vector<double> computeSourceTermVector(HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    double computeMaterial1TotalEnergy(HPRMaterialParameters material1Parameters);
    double computeMaterial1Temperature(HPRMaterialParameters material1Parameters);
    vector<double> computeMaterial1HeatFluxVector(HPRMaterialParameters material1Parameters);

    double computeMaterial1TotalEnergyDerivativeDensity(HPRMaterialParameters material1Parameters);
    double computeMaterial1TotalEnergyDerivativePressure(HPRMaterialParameters material1Parameters);

    vector<vector<double> > computeMaterial1TotalEnergyDerivativeDistortionTensor(HPRMaterialParameters material1Parameters);
    vector<double> computeMaterial1TotalEnergyDerivativeThermalImpulse(HPRMaterialParameters material1Parameters);

    double computeMaterial1TemperatureDerivativeDensity(HPRMaterialParameters material1Parameters);
    double computeMaterial1TemperatureDerivativePressure(HPRMaterialParameters material1Parameters);

    double computeMaterial1Theta1Reciprocal(HPRMaterialParameters material1Parameters);
    double computeMaterial1Theta2Reciprocal(HPRMaterialParameters material1Parameters);

    vector<vector<double> > computeMaterial1ShearStressTensor(HPRMaterialParameters material1Parameters);
    vector<vector<double> > computeMaterial1ShearStressTensorDerivativeDensity(HPRMaterialParameters material1Parameters);
    vector<vector<vector<vector<double> > > > computeMaterial1ShearStressTensorDerivativeDistortionTensor(HPRMaterialParameters material1Parameters);

    double computeMaterial2TotalEnergy(HPRMaterialParameters material2Parameters);
    double computeMaterial2Temperature(HPRMaterialParameters material2Parameters);
    vector<double> computeMaterial2HeatFluxVector(HPRMaterialParameters material2Parameters);

    double computeMaterial2TotalEnergyDerivativeDensity(HPRMaterialParameters material2Parameters);
    double computeMaterial2TotalEnergyDerivativePressure(HPRMaterialParameters material2Parameters);

    vector<vector<double> > computeMaterial2TotalEnergyDerivativeDistortionTensor(HPRMaterialParameters material2Parameters);
    vector<double> computeMaterial2TotalEnergyDerivativeThermalImpulse(HPRMaterialParameters material2Parameters);

    double computeMaterial2TemperatureDerivativeDensity(HPRMaterialParameters material2Parameters);
    double computeMaterial2TemperatureDerivativePressure(HPRMaterialParameters material2Parameters);

    double computeMaterial2Theta1Reciprocal(HPRMaterialParameters material2Parameters);
    double computeMaterial2Theta2Reciprocal(HPRMaterialParameters material2Parameters);

    vector<vector<double> > computeMaterial2ShearStressTensor(HPRMaterialParameters material2Parameters);
    vector<vector<double> > computeMaterial2ShearStressTensorDerivativeDensity(HPRMaterialParameters material2Parameters);
    vector<vector<vector<vector<double> > > > computeMaterial2ShearStressTensorDerivativeDistortionTensor(HPRMaterialParameters material2Parameters);

    double computeTotalDensity();
    double computeTotalPressure();
    vector<vector<double> > computeTotalDistortionTensor();

    double computeTotalXThermalImpulse();
    double computeTotalYThermalImpulse();
    double computeTotalZThermalImpulse();

    void relaxTotalDensity();
    void relaxTotalPressure();
    void relaxTotalDistortionTensor();

    void relaxTotalXThermalImpulse();
    void relaxTotalYThermalImpulse();
    void relaxTotalZThermalImpulse();

    void setMaterial1VolumeFraction(double newMaterial1VolumeFraction);

    void setInterfaceXVelocity(double newInterfaceXVelocity);
    void setInterfaceYVelocity(double newInterfaceYVelocity);
    void setInterfaceZVelocity(double newInterfaceZVelocity);

    void setMaterial1Density(double newMaterial1Density);
    void setMaterial1Pressure(double newMaterial1Pressure);
    void setMaterial1DistortionTensor(vector<vector<double> > newMaterial1DistortionTensor);

    void setMaterial1XThermalImpulse(double newMaterial1XThermalImpulse);
    void setMaterial1YThermalImpulse(double newMaterial1YThermalImpulse);
    void setMaterial1ZThermalImpulse(double newMaterial1ZThermalImpulse);

    void setMaterial2Density(double newMaterial2Density);
    void setMaterial2Pressure(double newMaterial2Pressure);
    void setMaterial2DistortionTensor(vector<vector<double> > newMaterial2DistortionTensor);

    void setMaterial2XThermalImpulse(double newMaterial2XThermalImpulse);
    void setMaterial2YThermalImpulse(double newMaterial2YThermalImpulse);
    void setMaterial2ZThermalImpulse(double newMaterial2ZThermalImpulse);

    double getMaterial1VolumeFraction();

    double getInterfaceXVelocity();
    double getInterfaceYVelocity();
    double getInterfaceZVelocity();

    double getMaterial1Density();
    double getMaterial1Pressure();
    vector<vector<double> > getMaterial1DistortionTensor();

    double getMaterial1XThermalImpulse();
    double getMaterial1YThermalImpulse();
    double getMaterial1ZThermalImpulse();

    double getMaterial2Density();
    double getMaterial2Pressure();
    vector<vector<double> > getMaterial2DistortionTensor();

    double getMaterial2XThermalImpulse();
    double getMaterial2YThermalImpulse();
    double getMaterial2ZThermalImpulse();

private:
    double material1VolumeFraction;

    double interfaceXVelocity;
    double interfaceYVelocity;
    double interfaceZVelocity;

    double material1Density;
    double material1Pressure;
    vector<vector<double> > material1DistortionTensor;

    double material1XThermalImpulse;
    double material1YThermalImpulse;
    double material1ZThermalImpulse;

    double material2Density;
    double material2Pressure;
    vector<vector<double> > material2DistortionTensor;

    double material2XThermalImpulse;
    double material2YThermalImpulse;
    double material2ZThermalImpulse;
};

#endif // HPRMULTIPHYSICSSTATEVECTOR_H

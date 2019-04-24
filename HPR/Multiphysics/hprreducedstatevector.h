#ifndef HPRREDUCEDSTATEVECTOR_H
#define HPRREDUCEDSTATEVECTOR_H

#include "HPR/hprequationofstate.h"
#include "HPR/hprsourceterms.h"
#include "Elasticity/elasticequationofstate.h"
using namespace std;

class HPRReducedStateVector
{
public:
    HPRReducedStateVector();
    HPRReducedStateVector(double newMaterial1VolumeFraction, double newInterfacePressure, vector<vector<double> > newInterfaceDistortionTensor, double newInterfaceXVelocity,
                          double newInterfaceYVelocity, double newInterfaceZVelocity, double newInterfaceXThermalImpulse, double newInterfaceYThermalImpulse, double newInterfaceZThermalImpulse,
                          double newMaterial1Density, double newMaterial2Density);
    HPRReducedStateVector(double newMaterial1VolumeFraction, double newInterfaceXVelocity, double newInterfaceYVelocity, double newInterfaceZVelocity,
                          vector<vector<double> > newMaterial1DistorionTensor, double newMaterial1Entropy, vector<vector<double> > newMaterial2DistortionTensor, double newMaterial2Entropy,
                          HyperelasticMaterialParameters hyperelasticMaterial1Parameters, HyperelasticMaterialParameters hyperelasticMaterial2Parameters, HPRMaterialParameters material1Parameters,
                          HPRMaterialParameters material2Parameters);

    void setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    void setConservedVariableVector(vector<double> newConservedVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    vector<double> computePrimitiveVariableVector(HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    vector<double> computeConservedVariableVector(HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static vector<double> computeXFluxVector(vector<double> conservedVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    vector<double> computeXFluxVector(HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static vector<double> computeYFluxVector(vector<double> conservedVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    vector<double> computeYFluxVector(HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static vector<double> computeSourceTermVector(vector<double> conservedVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    vector<double> computeSourceTermVector(HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    double computeMaterial1TotalEnergy(HPRMaterialParameters material1Parameters);
    double computeMaterial1Temperature(HPRMaterialParameters material1Parameters);
    vector<double> computeMaterial1HeatFluxVector(HPRMaterialParameters material1Parameters);

    double computeMaterial1TotalEnergyDerivativeDensity(HPRMaterialParameters material1Parameters);
    double computeMaterial1TotalEnergyDerivativePressure(HPRMaterialParameters material1Parameters);

    double computeMaterial1TemperatureDerivativeDensity(HPRMaterialParameters material1Parameters);
    double computeMaterial1TemperatureDerivativePressure(HPRMaterialParameters material1Parameters);

    double computeMaterial1Theta1Reciprocal(HPRMaterialParameters material1Parameters);
    double computeMaterial1Theta2Reciprocal(HPRMaterialParameters material1Parameters);

    vector<vector<double> > computeMaterial1ShearStressTensor(HPRMaterialParameters material1Parameters);
    vector<vector<vector<vector<double> > > > computeMaterial1ShearStressTensorDerivativeDistortionTensor(HPRMaterialParameters material1Parameters);

    double computeMaterial2TotalEnergy(HPRMaterialParameters material2Parameters);
    double computeMaterial2Temperature(HPRMaterialParameters material2Parameters);
    vector<double> computeMaterial2HeatFluxVector(HPRMaterialParameters material2Parameters);

    double computeMaterial2TotalEnergyDerivativeDensity(HPRMaterialParameters material2Parameters);
    double computeMaterial2TotalEnergyDerivativePressure(HPRMaterialParameters material2Parameters);

    double computeMaterial2TemperatureDerivativeDensity(HPRMaterialParameters material2Parameters);
    double computeMaterial2TemperatureDerivativePressure(HPRMaterialParameters material2Parameters);

    double computeMaterial2Theta1Reciprocal(HPRMaterialParameters material2Parameters);
    double computeMaterial2Theta2Reciprocal(HPRMaterialParameters material2Parameters);

    vector<vector<double> > computeMaterial2ShearStressTensor(HPRMaterialParameters material2Parameters);
    vector<vector<vector<vector<double> > > > computeMaterial2ShearStressTensorDerivativeDistortionTensor(HPRMaterialParameters material2Parameters);

    vector<vector<double> > computeTotalEnergyDerivativeDistortionTensor(HPRMaterialParameters materialParameters);
    vector<double> computeTotalEnergyDerivativeThermalImpulse(HPRMaterialParameters materialParameters);

    vector<vector<double> > computeShearStressTensorDerivativeDensity(HPRMaterialParameters materialParameters);

    double computeTotalDensity();
    void relaxTotalDensity();

    void setMaterial1VolumeFraction(double newMaterial1VolumeFraction);
    void setInterfacePressure(double newInterfacePressure);
    void setInterfaceDistortionTensor(vector<vector<double> > newInterfaceDistortionTensor);

    void setInterfaceXVelocity(double newInterfaceXVelocity);
    void setInterfaceYVelocity(double newInterfaceYVelocity);
    void setInterfaceZVelocity(double newInterfaceZVelocity);

    void setInterfaceXThermalImpulse(double newInterfaceXThermalImpulse);
    void setInterfaceYThermalImpulse(double newInterfaceYThermalImpulse);
    void setInterfaceZThermalImpulse(double newInterfaceZThermalImpulse);

    void setMaterial1Density(double newMaterial1Density);
    void setMaterial2Density(double newMaterial2Density);

    double getMaterial1VolumeFraction();
    double getInterfacePressure();
    vector<vector<double> > getInterfaceDistortionTensor();

    double getInterfaceXVelocity();
    double getInterfaceYVelocity();
    double getInterfaceZVelocity();

    double getInterfaceXThermalImpulse();
    double getInterfaceYThermalImpulse();
    double getInterfaceZThermalImpulse();

    double getMaterial1Density();
    double getMaterial2Density();

private:
    double material1VolumeFraction;
    double interfacePressure;
    vector<vector<double> > interfaceDistortionTensor;

    double interfaceXVelocity;
    double interfaceYVelocity;
    double interfaceZVelocity;

    double interfaceXThermalImpulse;
    double interfaceYThermalImpulse;
    double interfaceZThermalImpulse;

    double material1Density;
    double material2Density;
};

#endif // HPRREDUCEDSTATEVECTOR_H

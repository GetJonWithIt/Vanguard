#ifndef EULERMULTIPHYSICSSTATEVECTOR_H
#define EULERMULTIPHYSICSSTATEVECTOR_H

#include "Euler/eulerequationofstate.h"
#include <vector>
using namespace std;

class EulerMultiphysicsStateVector
{
public:
    EulerMultiphysicsStateVector();
    EulerMultiphysicsStateVector(double newMaterial1VolumeFraction, double newInterfaceXVelocity, double newInterfaceYVelocity, double newInterfaceZVelocity, double newMaterial1Density,
                                 double newMaterial1Pressure, double newMaterial2Density, double newMaterial2Pressure);

    void setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector);
    void setConservedVariableVector(vector<double> newConservedVariableVector, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    vector<double> computePrimitiveVariableVector();
    vector<double> computeConservedVariableVector(EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> computeXFluxVector(vector<double> conservedVariableVector, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    vector<double> computeXFluxVector(EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> computeYFluxVector(vector<double> conservedVariableVector, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    vector<double> computeYFluxVector(EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    double computeMaterial1SpecificInternalEnergy(EulerMaterialParameters material1Parameters);
    double computeMaterial1TotalEnergy(EulerMaterialParameters material1Parameters);
    double computeMaterial1SoundSpeed(EulerMaterialParameters material1Parameters);
    double computeMaterial1Entropy(EulerMaterialParameters material1Parameters);

    double computeMaterial2SpecificInternalEnergy(EulerMaterialParameters material2Parameters);
    double computeMaterial2TotalEnergy(EulerMaterialParameters material2Parameters);
    double computeMaterial2SoundSpeed(EulerMaterialParameters material2Parameters);
    double computeMaterial2Entropy(EulerMaterialParameters material2Parameters);

    double computeTotalDensity();
    double computeTotalPressure();

    void relaxTotalDensity();
    void relaxTotalPressure(EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    void setMaterial1VolumeFraction(double newMaterial1VolumeFraction);
    void setInterfaceXVelocity(double newInterfaceXVelocity);
    void setInterfaceYVelocity(double newInterfaceYVelocity);
    void setInterfaceZVelocity(double newInterfaceZVelocity);

    void setMaterial1Density(double newMaterial1Density);
    void setMaterial1Pressure(double newMaterial1Pressure);

    void setMaterial2Density(double newMaterial2Density);
    void setMaterial2Pressure(double newMaterial2Pressure);

    double getMaterial1VolumeFraction();
    double getInterfaceXVelocity();
    double getInterfaceYVelocity();
    double getInterfaceZVelocity();

    double getMaterial1Density();
    double getMaterial1Pressure();

    double getMaterial2Density();
    double getMaterial2Pressure();

private:
    double material1VolumeFraction;
    double interfaceXVelocity;
    double interfaceYVelocity;
    double interfaceZVelocity;

    double material1Density;
    double material1Pressure;

    double material2Density;
    double material2Pressure;
};

#endif // EULERMULTIPHYSICSSTATEVECTOR_H

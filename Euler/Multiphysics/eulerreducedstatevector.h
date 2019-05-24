#ifndef EULERREDUCEDSTATEVECTOR_H
#define EULERREDUCEDSTATEVECTOR_H

#include "Euler/eulerequationofstate.h"
#include <vector>
using namespace std;

class EulerReducedStateVector
{
public:
    EulerReducedStateVector();
    EulerReducedStateVector(double newMaterial1VolumeFraction, double newInterfaceXVelocity, double newInterfaceYVelocity, double newInterfaceZVelocity, double newMaterial1Density,
                            double newMaterial2Density, double newInterfacePressure);

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

    void relaxTotalDensity();

    void setMaterial1VolumeFraction(double newMaterial1VolumeFraction);
    void setInterfaceXVelocity(double newInterfaceXVelocity);
    void setInterfaceYVelocity(double newInterfaceYVelocity);
    void setInterfaceZVelocity(double newInterfaceZVelocity);

    void setMaterial1Density(double newMaterial1Density);
    void setMaterial2Density(double newMaterial2Density);

    void setInterfacePressure(double newInterfacePressure);

    double getMaterial1VolumeFraction();
    double getInterfaceXVelocity();
    double getInterfaceYVelocity();
    double getInterfaceZVelocity();

    double getMaterial1Density();
    double getMaterial2Density();

    double getInterfacePressure();

private:
    double material1VolumeFraction;
    double interfaceXVelocity;
    double interfaceYVelocity;
    double interfaceZVelocity;

    double material1Density;
    double material2Density;

    double interfacePressure;
};

#endif // EULERREDUCEDSTATEVECTOR_H

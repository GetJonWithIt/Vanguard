#ifndef MHDREDUCEDSTATEVECTOR_H
#define MHDREDUCEDSTATEVECTOR_H

#include "MHD/mhdwavespeeds.h"
using namespace std;

class MHDReducedStateVector
{
public:
    MHDReducedStateVector();
    MHDReducedStateVector(double newMaterial1VolumeFraction, double newInterfaceXVelocity, double newInterfaceYVelocity, double newInterfaceZVelocity, double newMaterial1Density,
                          double newMaterial2Density, double newInterfacePressure, double newInterfaceXMagneticField, double newInterfaceYMagneticField, double newInterfaceZMagneticField,
                          double newInterfaceAuxiliaryField);

    void setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector);
    void setConservedVariableVector(vector<double> newConservedVariableVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    vector<double> computePrimitiveVariableVector();
    vector<double> computeConservedVariableVector(MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static vector<double> computeXFluxVector(vector<double> conservedVariableVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    vector<double> computeXFluxVector(MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static vector<double> computeYFluxVector(vector<double> conservedVariableVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    vector<double> computeYFluxVector(MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static vector<double> computeSourceTermVector(vector<double> conservedVariableVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    vector<double> computeSourceTermVector(MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    double computeMaterial1SpecificInternalEnergy(MHDMaterialParameters material1Parameters);
    double computeMaterial1TotalEnergy(MHDMaterialParameters material1Parameters);
    double computeMaterial1SoundSpeed(MHDMaterialParameters material1Parameters);
    double computeMaterial1Entropy(MHDMaterialParameters material1Parameters);

    double computeMaterial1AlfvenWaveSpeed();

    double computeMaterial1XSlowMagnetoAcousticSpeed(MHDMaterialParameters material1Parameters);
    double computeMaterial1YSlowMagnetoAcousticSpeed(MHDMaterialParameters material1Parameters);

    double computeMaterial1XFastMagnetoAcousticSpeed(MHDMaterialParameters material1Parameters);
    double computeMaterial1YFastMagnetoAcousticSpeed(MHDMaterialParameters material1Parameters);

    double computeMaterial2SpecificInternalEnergy(MHDMaterialParameters material2Parameters);
    double computeMaterial2TotalEnergy(MHDMaterialParameters material2Parameters);
    double computeMaterial2SoundSpeed(MHDMaterialParameters material2Parameters);
    double computeMaterial2Entropy(MHDMaterialParameters material2Parameters);

    double computeMaterial2AlfvenWaveSpeed();

    double computeMaterial2XSlowMagnetoAcousticSpeed(MHDMaterialParameters material2Parameters);
    double computeMaterial2YSlowMagnetoAcousticSpeed(MHDMaterialParameters material2Parameters);

    double computeMaterial2XFastMagnetoAcousticSpeed(MHDMaterialParameters material2Parameters);
    double computeMaterial2YFastMagnetoAcousticSpeed(MHDMaterialParameters material2Parameters);

    double computeTotalDensity();

    void relaxTotalDensity();

    void setMaterial1VolumeFraction(double newMaterial1VolumeFraction);
    void setInterfaceXVelocity(double newInterfaceXVelocity);
    void setInterfaceYVelocity(double newInterfaceYVelocity);
    void setInterfaceZVelocity(double newInterfaceZVelocity);

    void setMaterial1Density(double newMaterial1Density);
    void setMaterial2Density(double newMaterial2Density);

    void setInterfacePressure(double newInterfacePressure);
    void setInterfaceXMagneticField(double newInterfaceXMagneticField);
    void setInterfaceYMagneticField(double newInterfaceYMagneticField);
    void setInterfaceZMagneticField(double newInterfaceZMagneticField);

    void setInterfaceAuxiliaryField(double newInterfaceAuxiliaryField);

    double getMaterial1VolumeFraction();
    double getInterfaceXVelocity();
    double getInterfaceYVelocity();
    double getInterfaceZVelocity();

    double getMaterial1Density();
    double getMaterial2Density();

    double getInterfacePressure();
    double getInterfaceXMagneticField();
    double getInterfaceYMagneticField();
    double getInterfaceZMagneticField();

    double getInterfaceAuxiliaryField();

private:
    double material1VolumeFraction;
    double interfaceXVelocity;
    double interfaceYVelocity;
    double interfaceZVelocity;

    double material1Density;
    double material2Density;

    double interfacePressure;
    double interfaceXMagneticField;
    double interfaceYMagneticField;
    double interfaceZMagneticField;

    double interfaceAuxiliaryField;
};

#endif // MHDREDUCEDSTATEVECTOR_H

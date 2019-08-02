#ifndef MHDMULTIPHYSICSSTATEVECTOR_H
#define MHDMULTIPHYSICSSTATEVECTOR_H

#include "MHD/mhdwavespeeds.h"
using namespace std;

class MHDMultiphysicsStateVector
{
public:
    MHDMultiphysicsStateVector();
    MHDMultiphysicsStateVector(double newMaterial1VolumeFraction, double newInterfaceXVelocity, double newInterfaceYVelocity, double newInterfaceZVelocity, double newMaterial1Density,
                               double newMaterial1Pressure, double newMaterial1XMagneticField, double newMaterial1YMagneticField, double newMaterial1ZMagneticField,
                               double newMaterial1AuxiliaryField, double newMaterial2Density, double newMaterial2Pressure, double newMaterial2XMagneticField, double newMaterial2YMagneticField,
                               double newMaterial2ZMagneticField, double newMaterial2AuxiliaryField);

    void setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector);
    void setConservedVariableVector(vector<double> newConservedVariableVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    vector<double> computePrimitiveVariableVector();
    vector<double> computeConservedVariableVector(MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static vector<double> computeXFluxVector(vector<double> conservedVariableVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    vector<double> computeXFluxVector(MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

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
    double computeTotalPressure();

    double computeTotalXMagneticField();
    double computeTotalYMagneticField();
    double computeTotalZMagneticField();

    double computeTotalAuxiliaryField();

    void relaxTotalDensity();
    void relaxTotalPressure();

    void relaxTotalXMagneticField();
    void relaxTotalYMagneticField();
    void relaxTotalZMagneticField();

    void relaxTotalAuxiliaryField();

    void setMaterial1VolumeFraction(double newMaterial1VolumeFraction);
    void setInterfaceXVelocity(double newInterfaceXVelocity);
    void setInterfaceYVelocity(double newInterfaceYVelocity);
    void setInterfaceZVelocity(double newInterfaceZVelocity);

    void setMaterial1Density(double newMaterial1Density);
    void setMaterial1Pressure(double newMaterial1Pressure);

    void setMaterial1XMagneticField(double newMaterial1XMagneticField);
    void setMaterial1YMagneticField(double newMaterial1YMagneticField);
    void setMaterial1ZMagneticField(double newMaterial1ZMagneticField);

    void setMaterial1AuxiliaryField(double newMaterial1AuxiliaryField);

    void setMaterial2Density(double newMaterial2Density);
    void setMaterial2Pressure(double newMaterial2Pressure);

    void setMaterial2XMagneticField(double newMaterial2XMagneticField);
    void setMaterial2YMagneticField(double newMaterial2YMagneticField);
    void setMaterial2ZMagneticField(double newMaterial2ZMagneticField);

    void setMaterial2AuxiliaryField(double newMaterial2AuxiliaryField);

    double getMaterial1VolumeFraction();
    double getInterfaceXVelocity();
    double getInterfaceYVelocity();
    double getInterfaceZVelocity();

    double getMaterial1Density();
    double getMaterial1Pressure();

    double getMaterial1XMagneticField();
    double getMaterial1YMagneticField();
    double getMaterial1ZMagneticField();

    double getMaterial1AuxiliaryField();

    double getMaterial2Density();
    double getMaterial2Pressure();

    double getMaterial2XMagneticField();
    double getMaterial2YMagneticField();
    double getMaterial2ZMagneticField();

    double getMaterial2AuxiliaryField();

private:
    double material1VolumeFraction;
    double interfaceXVelocity;
    double interfaceYVelocity;
    double interfaceZVelocity;

    double material1Density;
    double material1Pressure;

    double material1XMagneticField;
    double material1YMagneticField;
    double material1ZMagneticField;

    double material1AuxiliaryField;

    double material2Density;
    double material2Pressure;

    double material2XMagneticField;
    double material2YMagneticField;
    double material2ZMagneticField;

    double material2AuxiliaryField;
};

#endif // MHDMULTIPHYSICSSTATEVECTOR_H

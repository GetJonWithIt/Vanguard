#ifndef MHDSTATEVECTOR_H
#define MHDSTATEVECTOR_H

#include "mhdwavespeeds.h"
using namespace std;

class MHDStateVector
{
public:
    MHDStateVector();
    MHDStateVector(double newDensity, double newXVelocity, double newYVelocity, double newZVelocity, double newPressure, double newXMagneticField, double newYMagneticField,
                   double newZMagneticField, double newAuxiliaryField);

    void setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector);
    void setConservedVariableVector(vector<double> newConservedVariableVector, MHDMaterialParameters materialParameters);

    vector<double> computePrimitiveVariableVector();
    vector<double> computeConservedVariableVector(MHDMaterialParameters materialParameters);

    static vector<double> computeXFluxVector(vector<double> conservedVariableVector, MHDMaterialParameters materialParameters);
    vector<double> computeXFluxVector(MHDMaterialParameters materialParameters);

    static vector<double> computeSourceTermVector(vector<double> conservedVariableVector, MHDMaterialParameters materialParameters);
    vector<double> computeSourceTermVector(MHDMaterialParameters materialParameters);

    double computeSpecificInternalEnergy(MHDMaterialParameters materialParameters);
    double computeTotalEnergy(MHDMaterialParameters materialParameters);
    double computeSoundSpeed(MHDMaterialParameters materialParameters);
    double computeEntropy(MHDMaterialParameters materialParameters);

    double computeAlfvenWaveSpeed();
    double computeSlowMagnetoAcousticSpeed(MHDMaterialParameters materialParameters);
    double computeFastMagnetoAcousticSpeed(MHDMaterialParameters materialParameters);

    void setDensity(double newDensity);
    void setXVelocity(double newXVelocity);
    void setYVelocity(double newYVelocity);
    void setZVelocity(double newZVelocity);
    void setPressure(double newPressure);

    void setXMagneticField(double newXMagneticField);
    void setYMagneticField(double newYMagneticField);
    void setZMagneticField(double newZMagneticField);

    void setAuxiliaryField(double newAuxiliaryField);

    double getDensity();
    double getXVelocity();
    double getYVelocity();
    double getZVelocity();
    double getPressure();

    double getXMagneticField();
    double getYMagneticField();
    double getZMagneticField();

    double getAuxiliaryField();

private:
    double density;
    double xVelocity;
    double yVelocity;
    double zVelocity;
    double pressure;

    double xMagneticField;
    double yMagneticField;
    double zMagneticField;

    double auxiliaryField;
};

#endif // MHDSTATEVECTOR_H

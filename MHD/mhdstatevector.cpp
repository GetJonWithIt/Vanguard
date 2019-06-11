#include "mhdstatevector.h"

MHDStateVector::MHDStateVector()
{
    density = 1.0;
    xVelocity = 0.0;
    yVelocity = 0.0;
    zVelocity = 0.0;
    pressure = 1.0;

    xMagneticField = 0.0;
    yMagneticField = 0.0;
    zMagneticField = 0.0;

    auxiliaryField = 0.0;
}

MHDStateVector::MHDStateVector(double newDensity, double newXVelocity, double newYVelocity, double newZVelocity, double newPressure, double newXMagneticField, double newYMagneticField,
                               double newZMagneticField, double newAuxiliaryField)
{
    density = newDensity;
    xVelocity = newXVelocity;
    yVelocity = newYVelocity;
    zVelocity = newZVelocity;
    pressure = newPressure;

    xMagneticField = newXMagneticField;
    yMagneticField = newYMagneticField;
    zMagneticField = newZMagneticField;

    auxiliaryField = newAuxiliaryField;
}

void MHDStateVector::setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector)
{
    density = newPrimitiveVariableVector[0];
    xVelocity = newPrimitiveVariableVector[1];
    yVelocity = newPrimitiveVariableVector[2];
    zVelocity = newPrimitiveVariableVector[3];
    pressure = newPrimitiveVariableVector[4];

    xMagneticField = newPrimitiveVariableVector[5];
    yMagneticField = newPrimitiveVariableVector[6];
    zMagneticField = newPrimitiveVariableVector[7];

    auxiliaryField = newPrimitiveVariableVector[8];
}

void MHDStateVector::setConservedVariableVector(vector<double> newConservedVariableVector, MHDMaterialParameters materialParameters)
{
    density = newConservedVariableVector[0];
    xVelocity = newConservedVariableVector[1] / density;
    yVelocity = newConservedVariableVector[2] / density;
    zVelocity = newConservedVariableVector[3] / density;

    xMagneticField = newConservedVariableVector[5];
    yMagneticField = newConservedVariableVector[6];
    zMagneticField = newConservedVariableVector[7];

    double magneticFieldSquared = (xMagneticField * xMagneticField) + (yMagneticField * yMagneticField) + (zMagneticField * zMagneticField);
    double totalEnergy = (newConservedVariableVector[4] - (0.5 * magneticFieldSquared)) / density;

    double velocitySquared = (xVelocity * xVelocity) + (yVelocity * yVelocity) + (zVelocity * zVelocity);
    double specificInternalEnergy = totalEnergy - (0.5 * velocitySquared);
    pressure = MHDEquationOfState::computePressure(density, specificInternalEnergy, materialParameters);

    auxiliaryField = newConservedVariableVector[8];
}

vector<double> MHDStateVector::computePrimitiveVariableVector()
{
    vector<double> primitiveVariableVector(9);

    primitiveVariableVector[0] = density;
    primitiveVariableVector[1] = xVelocity;
    primitiveVariableVector[2] = yVelocity;
    primitiveVariableVector[3] = zVelocity;
    primitiveVariableVector[4] = pressure;

    primitiveVariableVector[5] = xMagneticField;
    primitiveVariableVector[6] = yMagneticField;
    primitiveVariableVector[7] = zMagneticField;

    primitiveVariableVector[8] = auxiliaryField;

    return primitiveVariableVector;
}

vector<double> MHDStateVector::computeConservedVariableVector(MHDMaterialParameters materialParameters)
{
    vector<double> conservedVariableVector(9);

    conservedVariableVector[0] = density;
    conservedVariableVector[1] = density * xVelocity;
    conservedVariableVector[2] = density * yVelocity;
    conservedVariableVector[3] = density * zVelocity;
    conservedVariableVector[4] = computeTotalEnergy(materialParameters);

    conservedVariableVector[5] = xMagneticField;
    conservedVariableVector[6] = yMagneticField;
    conservedVariableVector[7] = zMagneticField;

    conservedVariableVector[8] = auxiliaryField;

    return conservedVariableVector;
}

vector<double> MHDStateVector::computeXFluxVector(vector<double> conservedVariableVector, MHDMaterialParameters materialParameters)
{
    vector<double> fluxVector(9);

    double computedDensity = conservedVariableVector[0];
    double computedXVelocity = conservedVariableVector[1] / computedDensity;
    double computedYVelocity = conservedVariableVector[2] / computedDensity;
    double computedZVelocity = conservedVariableVector[3] / computedDensity;

    double computedXMagneticField = conservedVariableVector[5];
    double computedYMagneticField = conservedVariableVector[6];
    double computedZMagneticField = conservedVariableVector[7];

    double computedAuxiliaryField = conservedVariableVector[8];

    double magneticFieldSquared = (computedXMagneticField * computedXMagneticField) + (computedYMagneticField * computedYMagneticField) + (computedZMagneticField * computedZMagneticField);
    double totalEnergy = conservedVariableVector[4];
    double hydrodynamicTotalEnergy = (totalEnergy - (0.5 * magneticFieldSquared)) / computedDensity;

    double velocitySquared = (computedXVelocity * computedXVelocity) + (computedYVelocity * computedYVelocity) + (computedZVelocity * computedZVelocity);
    double specificInternalEnergy = hydrodynamicTotalEnergy - (0.5 * velocitySquared);
    double computedPressure = MHDEquationOfState::computePressure(computedDensity, specificInternalEnergy, materialParameters);

    double hyperbolicWaveSpeedSquared = materialParameters.computeHyperbolicWaveSpeedSquared();

    fluxVector[0] = computedDensity * computedXVelocity;
    fluxVector[1] = (computedDensity * (computedXVelocity * computedXVelocity)) + computedPressure + (0.5 * magneticFieldSquared) - (computedXMagneticField * computedXMagneticField);
    fluxVector[2] = (computedDensity * (computedXVelocity * computedYVelocity)) - (computedXMagneticField * computedYMagneticField);
    fluxVector[3] = (computedDensity * (computedXVelocity * computedZVelocity)) - (computedXMagneticField * computedZMagneticField);

    vector<double> velocityVector(3);
    velocityVector[0] = computedXVelocity;
    velocityVector[1] = computedYVelocity;
    velocityVector[2] = computedZVelocity;

    vector<double> magneticFieldVector(3);
    magneticFieldVector[0] = computedXMagneticField;
    magneticFieldVector[1] = computedYMagneticField;
    magneticFieldVector[2] = computedZMagneticField;

    double magneticFieldVelocityVectorProduct = VectorAlgebra::computeDotProduct(velocityVector, magneticFieldVector);
    fluxVector[4] = ((totalEnergy + computedPressure + (0.5 * magneticFieldSquared)) * computedXVelocity) - (magneticFieldVelocityVectorProduct * computedXMagneticField);

    fluxVector[5] = computedAuxiliaryField;
    fluxVector[6] = (computedYMagneticField * computedXVelocity) - (computedXMagneticField * computedYVelocity);
    fluxVector[7] = (computedZMagneticField * computedXVelocity) - (computedXMagneticField * computedZVelocity);

    fluxVector[8] = hyperbolicWaveSpeedSquared * computedXMagneticField;

    return fluxVector;
}

vector<double> MHDStateVector::computeXFluxVector(MHDMaterialParameters materialParameters)
{
    return computeXFluxVector(computeConservedVariableVector(materialParameters), materialParameters);
}

vector<double> MHDStateVector::computeYFluxVector(vector<double> conservedVariableVector, MHDMaterialParameters materialParameters)
{
    vector<double> fluxVector(9);

    double computedDensity = conservedVariableVector[0];
    double computedXVelocity = conservedVariableVector[1] / computedDensity;
    double computedYVelocity = conservedVariableVector[2] / computedDensity;
    double computedZVelocity = conservedVariableVector[3] / computedDensity;

    double computedXMagneticField = conservedVariableVector[5];
    double computedYMagneticField = conservedVariableVector[6];
    double computedZMagneticField = conservedVariableVector[7];

    double computedAuxiliaryField = conservedVariableVector[8];

    double magneticFieldSquared = (computedXMagneticField * computedXMagneticField) + (computedYMagneticField * computedYMagneticField) + (computedZMagneticField * computedZMagneticField);
    double totalEnergy = conservedVariableVector[4];
    double hydrodynamicTotalEnergy = (totalEnergy - (0.5 * magneticFieldSquared)) / computedDensity;

    double velocitySquared = (computedXVelocity * computedXVelocity) + (computedYVelocity * computedYVelocity) + (computedZVelocity * computedZVelocity);
    double specificInternalEnergy = hydrodynamicTotalEnergy - (0.5 * velocitySquared);
    double computedPressure = MHDEquationOfState::computePressure(computedDensity, specificInternalEnergy, materialParameters);

    double hyperbolicWaveSpeedSquared = materialParameters.computeHyperbolicWaveSpeedSquared();

    fluxVector[0] = computedDensity * computedYVelocity;
    fluxVector[1] = (computedDensity * (computedYVelocity * computedXVelocity)) - (computedYMagneticField * computedXMagneticField);
    fluxVector[2] = (computedDensity * (computedYVelocity * computedYVelocity)) + computedPressure + (0.5 * magneticFieldSquared) - (computedYMagneticField * computedYMagneticField);
    fluxVector[3] = (computedDensity * (computedYVelocity * computedZVelocity)) - (computedYMagneticField * computedZMagneticField);

    vector<double> velocityVector(3);
    velocityVector[0] = computedXVelocity;
    velocityVector[1] = computedYVelocity;
    velocityVector[2] = computedZVelocity;

    vector<double> magneticFieldVector(3);
    magneticFieldVector[0] = computedXMagneticField;
    magneticFieldVector[1] = computedYMagneticField;
    magneticFieldVector[2] = computedZMagneticField;

    double magneticFieldVelocityVectorProduct = VectorAlgebra::computeDotProduct(velocityVector, magneticFieldVector);
    fluxVector[4] = ((totalEnergy + computedPressure + (0.5 * magneticFieldSquared)) * computedYVelocity) - (magneticFieldVelocityVectorProduct * computedYMagneticField);

    fluxVector[5] = (computedXMagneticField * computedYVelocity) - (computedYMagneticField * computedXVelocity);
    fluxVector[6] = computedAuxiliaryField;
    fluxVector[7] = (computedZMagneticField * computedYVelocity) - (computedYMagneticField * computedZVelocity);

    fluxVector[8] = hyperbolicWaveSpeedSquared * computedYMagneticField;

    return fluxVector;
}

vector<double> MHDStateVector::computeYFluxVector(MHDMaterialParameters materialParameters)
{
    return computeYFluxVector(computeConservedVariableVector(materialParameters), materialParameters);
}

vector<double> MHDStateVector::computeSourceTermVector(vector<double> conservedVariableVector, MHDMaterialParameters materialParameters)
{
    vector<double> sourceTermVector(9);

    double hyperbolicWaveSpeedSquared = materialParameters.computeHyperbolicWaveSpeedSquared();
    double parabolicDampingSquared = materialParameters.computeParabolicDampingSquared();

    double computedAuxiliaryField = conservedVariableVector[8];

    sourceTermVector[8] = - (hyperbolicWaveSpeedSquared / parabolicDampingSquared) * computedAuxiliaryField;

    return sourceTermVector;
}

vector<double> MHDStateVector::computeSourceTermVector(MHDMaterialParameters materialParameters)
{
    return computeSourceTermVector(computeConservedVariableVector(materialParameters), materialParameters);
}

double MHDStateVector::computeSpecificInternalEnergy(MHDMaterialParameters materialParameters)
{
    return MHDEquationOfState::computeSpecificInternalEnergy(density, pressure, materialParameters);
}

double MHDStateVector::computeTotalEnergy(MHDMaterialParameters materialParameters)
{
    double velocitySquared = (xVelocity * xVelocity) + (yVelocity * yVelocity) + (zVelocity * zVelocity);
    double totalHydrodynamicEnergy = density * ((0.5 * velocitySquared) + computeSpecificInternalEnergy(materialParameters));

    double magneticFieldSquared = (xMagneticField * xMagneticField) + (yMagneticField * yMagneticField) + (zMagneticField * zMagneticField);

    return totalHydrodynamicEnergy + (0.5 * magneticFieldSquared);
}

double MHDStateVector::computeSoundSpeed(MHDMaterialParameters materialParameters)
{
    return MHDEquationOfState::computeSoundSpeed(density, pressure, materialParameters);
}

double MHDStateVector::computeEntropy(MHDMaterialParameters materialParameters)
{
    return MHDEquationOfState::computeEntropy(density, pressure, materialParameters);
}

double MHDStateVector::computeAlfvenWaveSpeed()
{
    return MHDWaveSpeeds::computeAlfvenWaveSpeed(density, xMagneticField, yMagneticField, zMagneticField);
}

double MHDStateVector::computeXSlowMagnetoAcousticSpeed(MHDMaterialParameters materialParameters)
{
    return MHDWaveSpeeds::computeXSlowMagnetoAcousticSpeed(density, pressure, xMagneticField, yMagneticField, zMagneticField, materialParameters);
}

double MHDStateVector::computeYSlowMagnetoAcousticSpeed(MHDMaterialParameters materialParameters)
{
    return MHDWaveSpeeds::computeYSlowMagnetoAcousticSpeed(density, pressure, xMagneticField, yMagneticField, zMagneticField, materialParameters);
}

double MHDStateVector::computeXFastMagnetoAcousticSpeed(MHDMaterialParameters materialParameters)
{
    return MHDWaveSpeeds::computeXFastMagnetoAcousticSpeed(density, pressure, xMagneticField, yMagneticField, zMagneticField, materialParameters);
}

double MHDStateVector::computeYFastMagnetoAcousticSpeed(MHDMaterialParameters materialParameters)
{
    return MHDWaveSpeeds::computeYFastMagnetoAcousticSpeed(density, pressure, xMagneticField, yMagneticField, zMagneticField, materialParameters);
}

void MHDStateVector::setDensity(double newDensity)
{
    density = newDensity;
}

void MHDStateVector::setXVelocity(double newXVelocity)
{
    xVelocity = newXVelocity;
}

void MHDStateVector::setYVelocity(double newYVelocity)
{
    yVelocity = newYVelocity;
}

void MHDStateVector::setZVelocity(double newZVelocity)
{
    zVelocity = newZVelocity;
}

void MHDStateVector::setPressure(double newPressure)
{
    pressure = newPressure;
}

void MHDStateVector::setXMagneticField(double newXMagneticField)
{
    xMagneticField = newXMagneticField;
}

void MHDStateVector::setYMagneticField(double newYMagneticField)
{
    yMagneticField = newYMagneticField;
}

void MHDStateVector::setZMagneticField(double newZMagneticField)
{
    zMagneticField = newZMagneticField;
}

void MHDStateVector::setAuxiliaryField(double newAuxiliaryField)
{
    auxiliaryField = newAuxiliaryField;
}

double MHDStateVector::getDensity()
{
    return density;
}

double MHDStateVector::getXVelocity()
{
    return xVelocity;
}

double MHDStateVector::getYVelocity()
{
    return yVelocity;
}

double MHDStateVector::getZVelocity()
{
    return zVelocity;
}

double MHDStateVector::getPressure()
{
    return pressure;
}

double MHDStateVector::getXMagneticField()
{
    return xMagneticField;
}

double MHDStateVector::getYMagneticField()
{
    return yMagneticField;
}

double MHDStateVector::getZMagneticField()
{
    return zMagneticField;
}

double MHDStateVector::getAuxiliaryField()
{
    return auxiliaryField;
}

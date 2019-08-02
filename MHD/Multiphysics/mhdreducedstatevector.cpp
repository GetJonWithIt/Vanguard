#include "mhdreducedstatevector.h"

MHDReducedStateVector::MHDReducedStateVector()
{
    material1VolumeFraction = 0.999;
    interfaceXVelocity = 0.0;
    interfaceYVelocity = 0.0;
    interfaceZVelocity = 0.0;

    material1Density = 1.0;
    material2Density = 1.0;

    interfacePressure = 1.0;
    interfaceXMagneticField = 0.0;
    interfaceYMagneticField = 0.0;
    interfaceZMagneticField = 0.0;

    interfaceAuxiliaryField = 0.0;
}

MHDReducedStateVector::MHDReducedStateVector(double newMaterial1VolumeFraction, double newInterfaceXVelocity, double newInterfaceYVelocity, double newInterfaceZVelocity,
                                             double newMaterial1Density, double newMaterial2Density, double newInterfacePressure, double newInterfaceXMagneticField,
                                             double newInterfaceYMagneticField, double newInterfaceZMagneticField, double newInterfaceAuxiliaryField)
{
    material1VolumeFraction = newMaterial1VolumeFraction;
    interfaceXVelocity = newInterfaceXVelocity;
    interfaceYVelocity = newInterfaceYVelocity;
    interfaceZVelocity = newInterfaceZVelocity;

    material1Density = newMaterial1Density;
    material2Density = newMaterial2Density;

    interfacePressure = newInterfacePressure;
    interfaceXMagneticField = newInterfaceXMagneticField;
    interfaceYMagneticField = newInterfaceYMagneticField;
    interfaceZMagneticField = newInterfaceZMagneticField;

    interfaceAuxiliaryField = newInterfaceAuxiliaryField;
}

void MHDReducedStateVector::setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector)
{
    material1VolumeFraction = newPrimitiveVariableVector[0];
    interfaceXVelocity = newPrimitiveVariableVector[1];
    interfaceYVelocity = newPrimitiveVariableVector[2];
    interfaceZVelocity = newPrimitiveVariableVector[3];

    material1Density = newPrimitiveVariableVector[4];
    material2Density = newPrimitiveVariableVector[5];

    interfacePressure = newPrimitiveVariableVector[6];
    interfaceXMagneticField = newPrimitiveVariableVector[7];
    interfaceYMagneticField = newPrimitiveVariableVector[8];
    interfaceZMagneticField = newPrimitiveVariableVector[9];

    interfaceAuxiliaryField = newPrimitiveVariableVector[10];
}

void MHDReducedStateVector::setConservedVariableVector(vector<double> newConservedVariableVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double totalDensity = newConservedVariableVector[0];

    if ((newConservedVariableVector[1] / totalDensity) < 0.001)
    {
        material1VolumeFraction = 0.001;
    }
    else if ((newConservedVariableVector[1] / totalDensity) > 0.999)
    {
        material1VolumeFraction = 0.999;
    }
    else
    {
        material1VolumeFraction = newConservedVariableVector[1] / totalDensity;
    }
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    interfaceXVelocity = newConservedVariableVector[2] / totalDensity;
    interfaceYVelocity = newConservedVariableVector[3] / totalDensity;
    interfaceZVelocity = newConservedVariableVector[4] / totalDensity;

    material1Density = newConservedVariableVector[5] / material1VolumeFraction;
    material2Density = newConservedVariableVector[6] / material2VolumeFraction;

    interfaceXMagneticField = newConservedVariableVector[8];
    interfaceYMagneticField = newConservedVariableVector[9];
    interfaceZMagneticField = newConservedVariableVector[10];

    double magneticFieldSquared = (interfaceXMagneticField * interfaceXMagneticField) + (interfaceYMagneticField * interfaceYMagneticField) + (interfaceZMagneticField * interfaceZMagneticField);
    double totalEnergy = (newConservedVariableVector[7] - (0.5 * magneticFieldSquared)) / totalDensity;

    double velocitySquared = (interfaceXVelocity * interfaceXVelocity) + (interfaceYVelocity * interfaceYVelocity) + (interfaceZVelocity * interfaceZVelocity);
    double totalSpecificInternalEnergy = totalEnergy - (0.5 * velocitySquared);

    double material1Pressure = MHDEquationOfState::computePressure(totalDensity, totalSpecificInternalEnergy, material1Parameters);
    double material2Pressure = MHDEquationOfState::computePressure(totalDensity, totalSpecificInternalEnergy, material2Parameters);
    interfacePressure = (material1VolumeFraction * material1Pressure) + (material2VolumeFraction * material2Pressure);

    interfaceAuxiliaryField = newConservedVariableVector[11];
}

vector<double> MHDReducedStateVector::computePrimitiveVariableVector()
{
    vector<double> primitiveVariableVector(11);

    primitiveVariableVector[0] = material1VolumeFraction;
    primitiveVariableVector[1] = interfaceXVelocity;
    primitiveVariableVector[2] = interfaceYVelocity;
    primitiveVariableVector[3] = interfaceZVelocity;

    primitiveVariableVector[4] = material1Density;
    primitiveVariableVector[5] = material2Density;

    primitiveVariableVector[6] = interfacePressure;
    primitiveVariableVector[7] = interfaceXMagneticField;
    primitiveVariableVector[8] = interfaceYMagneticField;
    primitiveVariableVector[9] = interfaceZMagneticField;

    primitiveVariableVector[10] = interfaceAuxiliaryField;

    return primitiveVariableVector;
}

vector<double> MHDReducedStateVector::computeConservedVariableVector(MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> conservedVariableVector(12);

    if (material1VolumeFraction < 0.001)
    {
        material1VolumeFraction = 0.001;
    }
    else if (material1VolumeFraction > 0.999)
    {
        material1VolumeFraction = 0.999;
    }
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    double totalDensity = (material1VolumeFraction * material1Density) + (material2VolumeFraction * material2Density);

    conservedVariableVector[0] = totalDensity;
    conservedVariableVector[1] = totalDensity * material1VolumeFraction;
    conservedVariableVector[2] = totalDensity * interfaceXVelocity;
    conservedVariableVector[3] = totalDensity * interfaceYVelocity;
    conservedVariableVector[4] = totalDensity * interfaceZVelocity;

    conservedVariableVector[5] = material1VolumeFraction * material1Density;
    conservedVariableVector[6] = material2VolumeFraction * material2Density;

    double velocitySquared = (interfaceXVelocity * interfaceXVelocity) + (interfaceYVelocity * interfaceYVelocity) + (interfaceZVelocity * interfaceZVelocity);
    double magneticFieldSquared = (interfaceXMagneticField * interfaceXMagneticField) + (interfaceYMagneticField * interfaceYMagneticField) + (interfaceZMagneticField * interfaceZMagneticField);

    double material1TotalHydrodynamicEnergy = totalDensity * ((0.5 * velocitySquared) + MHDEquationOfState::computeSpecificInternalEnergy(totalDensity, interfacePressure, material1Parameters));
    double material2TotalHydrodynamicEnergy = totalDensity * ((0.5 * velocitySquared) + MHDEquationOfState::computeSpecificInternalEnergy(totalDensity, interfacePressure, material2Parameters));

    double material1TotalEnergy = material1TotalHydrodynamicEnergy + (0.5 * magneticFieldSquared);
    double material2TotalEnergy = material2TotalHydrodynamicEnergy + (0.5 * magneticFieldSquared);

    double computedTotalEnergy = (material1VolumeFraction * material1TotalEnergy) + (material2VolumeFraction * material2TotalEnergy);
    conservedVariableVector[7] = computedTotalEnergy;

    conservedVariableVector[8] = interfaceXMagneticField;
    conservedVariableVector[9] = interfaceYMagneticField;
    conservedVariableVector[10] = interfaceZMagneticField;

    conservedVariableVector[11] = interfaceAuxiliaryField;

    return conservedVariableVector;
}

vector<double> MHDReducedStateVector::computeXFluxVector(vector<double> conservedVariableVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> fluxVector(12);

    double computedTotalDensity = conservedVariableVector[0];

    double computedMaterial1VolumeFraction;
    if ((conservedVariableVector[1] / computedTotalDensity) < 0.001)
    {
        computedMaterial1VolumeFraction = 0.001;
    }
    else if ((conservedVariableVector[1] / computedTotalDensity) > 0.999)
    {
        computedMaterial1VolumeFraction = 0.999;
    }
    else
    {
        computedMaterial1VolumeFraction = conservedVariableVector[1] / computedTotalDensity;
    }
    double computedMaterial2VolumeFraction = 1.0 - computedMaterial1VolumeFraction;

    double computedInterfaceXVelocity = conservedVariableVector[2] / computedTotalDensity;
    double computedInterfaceYVelocity = conservedVariableVector[3] / computedTotalDensity;
    double computedInterfaceZVelocity = conservedVariableVector[4] / computedTotalDensity;

    double computedMaterial1Density = conservedVariableVector[5] / computedMaterial1VolumeFraction;
    double computedMaterial2Density = conservedVariableVector[6] / computedMaterial2VolumeFraction;

    double computedInterfaceXMagneticField = conservedVariableVector[8];
    double computedInterfaceYMagneticField = conservedVariableVector[9];
    double computedInterfaceZMagneticField = conservedVariableVector[10];

    double computedInterfaceAuxiliaryField = conservedVariableVector[11];

    double velocitySquared = (computedInterfaceXVelocity * computedInterfaceXVelocity) + (computedInterfaceYVelocity * computedInterfaceYVelocity) +
            (computedInterfaceZVelocity * computedInterfaceZVelocity);
    double magneticFieldSquared = (computedInterfaceXMagneticField * computedInterfaceXMagneticField) + (computedInterfaceYMagneticField * computedInterfaceYMagneticField) +
            (computedInterfaceZMagneticField * computedInterfaceZMagneticField);

    double totalEnergy = conservedVariableVector[7];
    double hydrodynamicTotalEnergy = (totalEnergy - (0.5 * magneticFieldSquared)) / computedTotalDensity;
    double totalSpecificInternalEnergy = hydrodynamicTotalEnergy - (0.5 * velocitySquared);

    double computedMaterial1Pressure = MHDEquationOfState::computePressure(computedTotalDensity, totalSpecificInternalEnergy, material1Parameters);
    double computedMaterial2Pressure = MHDEquationOfState::computePressure(computedTotalDensity, totalSpecificInternalEnergy, material2Parameters);
    double computedInterfacePressure = (computedMaterial1VolumeFraction * computedMaterial1Pressure) + (computedMaterial2VolumeFraction * computedMaterial2Pressure);

    double material1HyperbolicWaveSpeedSquared = material1Parameters.computeHyperbolicWaveSpeedSquared();
    double material2HyperbolicWaveSpeedSquared = material2Parameters.computeHyperbolicWaveSpeedSquared();
    double hyperbolicWaveSpeedSquared = (computedMaterial1VolumeFraction * material1HyperbolicWaveSpeedSquared) + (computedMaterial2VolumeFraction * material2HyperbolicWaveSpeedSquared);

    fluxVector[0] = computedTotalDensity * computedInterfaceXVelocity;
    fluxVector[1] = computedTotalDensity * (computedInterfaceXVelocity * computedMaterial1VolumeFraction);

    fluxVector[2] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceXVelocity)) + computedInterfacePressure + (0.5 * magneticFieldSquared) -
            (computedInterfaceXMagneticField * computedInterfaceXMagneticField);
    fluxVector[3] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceYVelocity)) - (computedInterfaceXMagneticField * computedInterfaceYMagneticField);
    fluxVector[4] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceZVelocity)) - (computedInterfaceXMagneticField * computedInterfaceZMagneticField);

    fluxVector[5] = computedMaterial1VolumeFraction * (computedMaterial1Density * computedInterfaceXVelocity);
    fluxVector[6] = computedMaterial2VolumeFraction * (computedMaterial2Density * computedInterfaceXVelocity);

    vector<double> velocityVector(3);
    velocityVector[0] = computedInterfaceXVelocity;
    velocityVector[1] = computedInterfaceYVelocity;
    velocityVector[2] = computedInterfaceZVelocity;

    vector<double> magneticFieldVector(3);
    magneticFieldVector[0] = computedInterfaceXMagneticField;
    magneticFieldVector[1] = computedInterfaceYMagneticField;
    magneticFieldVector[2] = computedInterfaceZMagneticField;

    double magneticFieldVelocityVectorProduct = VectorAlgebra::computeDotProduct(velocityVector, magneticFieldVector);
    fluxVector[7] = ((totalEnergy + computedInterfacePressure + (0.5 * magneticFieldSquared)) * computedInterfaceXVelocity) -
            (magneticFieldVelocityVectorProduct * computedInterfaceXMagneticField);

    fluxVector[8] = computedInterfaceAuxiliaryField;
    fluxVector[9] = (computedInterfaceYMagneticField * computedInterfaceXVelocity) - (computedInterfaceXMagneticField * computedInterfaceYVelocity);
    fluxVector[10] = (computedInterfaceZMagneticField * computedInterfaceXVelocity) - (computedInterfaceXMagneticField * computedInterfaceZVelocity);

    fluxVector[11] = hyperbolicWaveSpeedSquared * computedInterfaceXMagneticField;

    return fluxVector;
}

vector<double> MHDReducedStateVector::computeXFluxVector(MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    return computeXFluxVector(computeConservedVariableVector(material1Parameters, material2Parameters), material1Parameters, material2Parameters);
}

vector<double> MHDReducedStateVector::computeYFluxVector(vector<double> conservedVariableVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> fluxVector(12);

    double computedTotalDensity = conservedVariableVector[0];

    double computedMaterial1VolumeFraction;
    if ((conservedVariableVector[1] / computedTotalDensity) < 0.001)
    {
        computedMaterial1VolumeFraction = 0.001;
    }
    else if ((conservedVariableVector[1] / computedTotalDensity) > 0.999)
    {
        computedMaterial1VolumeFraction = 0.999;
    }
    else
    {
        computedMaterial1VolumeFraction = conservedVariableVector[1] / computedTotalDensity;
    }
    double computedMaterial2VolumeFraction = 1.0 - computedMaterial1VolumeFraction;

    double computedInterfaceXVelocity = conservedVariableVector[2] / computedTotalDensity;
    double computedInterfaceYVelocity = conservedVariableVector[3] / computedTotalDensity;
    double computedInterfaceZVelocity = conservedVariableVector[4] / computedTotalDensity;

    double computedMaterial1Density = conservedVariableVector[5] / computedMaterial1VolumeFraction;
    double computedMaterial2Density = conservedVariableVector[6] / computedMaterial2VolumeFraction;

    double computedInterfaceXMagneticField = conservedVariableVector[8];
    double computedInterfaceYMagneticField = conservedVariableVector[9];
    double computedInterfaceZMagneticField = conservedVariableVector[10];

    double computedInterfaceAuxiliaryField = conservedVariableVector[11];

    double velocitySquared = (computedInterfaceXVelocity * computedInterfaceXVelocity) + (computedInterfaceYVelocity * computedInterfaceYVelocity) +
            (computedInterfaceZVelocity * computedInterfaceZVelocity);
    double magneticFieldSquared = (computedInterfaceXMagneticField * computedInterfaceXMagneticField) + (computedInterfaceYMagneticField * computedInterfaceYMagneticField) +
            (computedInterfaceZMagneticField * computedInterfaceZMagneticField);

    double totalEnergy = conservedVariableVector[7];
    double hydrodynamicTotalEnergy = (totalEnergy - (0.5 * magneticFieldSquared)) / computedTotalDensity;
    double totalSpecificInternalEnergy = hydrodynamicTotalEnergy - (0.5 * velocitySquared);

    double computedMaterial1Pressure = MHDEquationOfState::computePressure(computedTotalDensity, totalSpecificInternalEnergy, material1Parameters);
    double computedMaterial2Pressure = MHDEquationOfState::computePressure(computedTotalDensity, totalSpecificInternalEnergy, material2Parameters);
    double computedInterfacePressure = (computedMaterial1VolumeFraction * computedMaterial1Pressure) + (computedMaterial2VolumeFraction * computedMaterial2Pressure);

    double material1HyperbolicWaveSpeedSquared = material1Parameters.computeHyperbolicWaveSpeedSquared();
    double material2HyperbolicWaveSpeedSquared = material2Parameters.computeHyperbolicWaveSpeedSquared();
    double hyperbolicWaveSpeedSquared = (computedMaterial1VolumeFraction * material1HyperbolicWaveSpeedSquared) + (computedMaterial2VolumeFraction * material2HyperbolicWaveSpeedSquared);

    fluxVector[0] = computedTotalDensity * computedInterfaceYVelocity;
    fluxVector[1] = computedTotalDensity * (computedInterfaceYVelocity * computedMaterial1VolumeFraction);

    fluxVector[2] = (computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceXVelocity)) - (computedInterfaceYMagneticField * computedInterfaceXMagneticField);
    fluxVector[3] = (computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceYVelocity)) + computedInterfacePressure + (0.5 * magneticFieldSquared) -
            (computedInterfaceYMagneticField * computedInterfaceYMagneticField);
    fluxVector[4] = (computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceZVelocity)) - (computedInterfaceYMagneticField * computedInterfaceZMagneticField);

    fluxVector[5] = computedMaterial1VolumeFraction * (computedMaterial1Density * computedInterfaceYVelocity);
    fluxVector[6] = computedMaterial2VolumeFraction * (computedMaterial2Density * computedInterfaceYVelocity);

    vector<double> velocityVector(3);
    velocityVector[0] = computedInterfaceXVelocity;
    velocityVector[1] = computedInterfaceYVelocity;
    velocityVector[2] = computedInterfaceZVelocity;

    vector<double> magneticFieldVector(3);
    magneticFieldVector[0] = computedInterfaceXMagneticField;
    magneticFieldVector[1] = computedInterfaceYMagneticField;
    magneticFieldVector[2] = computedInterfaceZMagneticField;

    double magneticFieldVelocityVectorProduct = VectorAlgebra::computeDotProduct(velocityVector, magneticFieldVector);
    fluxVector[7] = ((totalEnergy + computedInterfacePressure + (0.5 * magneticFieldSquared)) * computedInterfaceYVelocity) -
            (magneticFieldVelocityVectorProduct * computedInterfaceYMagneticField);

    fluxVector[8] = (computedInterfaceXMagneticField * computedInterfaceYVelocity) - (computedInterfaceYMagneticField * computedInterfaceXVelocity);
    fluxVector[9] = computedInterfaceAuxiliaryField;
    fluxVector[10] = (computedInterfaceZMagneticField * computedInterfaceYVelocity) - (computedInterfaceYMagneticField * computedInterfaceZVelocity);

    fluxVector[11] = hyperbolicWaveSpeedSquared * computedInterfaceYMagneticField;

    return fluxVector;
}

vector<double> MHDReducedStateVector::computeYFluxVector(MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    return computeYFluxVector(computeConservedVariableVector(material1Parameters, material2Parameters), material1Parameters, material2Parameters);
}

vector<double> MHDReducedStateVector::computeSourceTermVector(vector<double> conservedVariableVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> sourceTermVector(12);

    double computedTotalDensity = conservedVariableVector[0];

    double computedMaterial1VolumeFraction;
    if ((conservedVariableVector[1] / computedTotalDensity) < 0.001)
    {
        computedMaterial1VolumeFraction = 0.001;
    }
    else if ((conservedVariableVector[1] / computedTotalDensity) > 0.999)
    {
        computedMaterial1VolumeFraction = 0.999;
    }
    else
    {
        computedMaterial1VolumeFraction = conservedVariableVector[1] / computedTotalDensity;
    }
    double computedMaterial2VolumeFraction = 1.0 - computedMaterial1VolumeFraction;

    double material1HyperbolicWaveSpeedSquared = material1Parameters.computeHyperbolicWaveSpeedSquared();
    double material2HyperbolicWaveSpeedSquared = material2Parameters.computeHyperbolicWaveSpeedSquared();
    double hyperbolicWaveSpeedSquared = (computedMaterial1VolumeFraction * material1HyperbolicWaveSpeedSquared) + (computedMaterial2VolumeFraction * material2HyperbolicWaveSpeedSquared);

    double material1ParabolicDampingSquared = material1Parameters.computeParabolicDampingSquared();
    double material2ParabolicDampingSquared = material2Parameters.computeParabolicDampingSquared();
    double parabolicDampingSquared = (computedMaterial1VolumeFraction * material1ParabolicDampingSquared) + (computedMaterial2VolumeFraction * material2ParabolicDampingSquared);

    double computedInterfaceAuxiliaryField = conservedVariableVector[11];

    sourceTermVector[11] = - (hyperbolicWaveSpeedSquared / parabolicDampingSquared) * computedInterfaceAuxiliaryField;

    return sourceTermVector;
}

vector<double> MHDReducedStateVector::computeSourceTermVector(MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    return computeSourceTermVector(computeConservedVariableVector(material1Parameters, material2Parameters), material1Parameters, material2Parameters);
}

double MHDReducedStateVector::computeMaterial1SpecificInternalEnergy(MHDMaterialParameters material1Parameters)
{
    return MHDEquationOfState::computeSpecificInternalEnergy(material1Density, interfacePressure, material1Parameters);
}

double MHDReducedStateVector::computeMaterial1TotalEnergy(MHDMaterialParameters material1Parameters)
{
    double velocitySquared = (interfaceXVelocity * interfaceXVelocity) + (interfaceYVelocity * interfaceYVelocity) + (interfaceZVelocity * interfaceZVelocity);
    double totalHydrodynamicEnergy = material1Density * ((0.5 * velocitySquared) + computeMaterial1SpecificInternalEnergy(material1Parameters));

    double magneticFieldSquared = (interfaceXMagneticField * interfaceXMagneticField) + (interfaceYMagneticField * interfaceYMagneticField) + (interfaceZMagneticField * interfaceZMagneticField);

    return totalHydrodynamicEnergy + (0.5 * magneticFieldSquared);
}

double MHDReducedStateVector::computeMaterial1SoundSpeed(MHDMaterialParameters material1Parameters)
{
    return MHDEquationOfState::computeSoundSpeed(material1Density, interfacePressure, material1Parameters);
}

double MHDReducedStateVector::computeMaterial1Entropy(MHDMaterialParameters material1Parameters)
{
    return MHDEquationOfState::computeEntropy(material1Density, interfacePressure, material1Parameters);
}

double MHDReducedStateVector::computeMaterial1AlfvenWaveSpeed()
{
    return MHDWaveSpeeds::computeAlfvenWaveSpeed(material1Density, interfaceXMagneticField, interfaceYMagneticField, interfaceZMagneticField);
}

double MHDReducedStateVector::computeMaterial1XSlowMagnetoAcousticSpeed(MHDMaterialParameters material1Parameters)
{
    return MHDWaveSpeeds::computeXSlowMagnetoAcousticSpeed(material1Density, interfacePressure, interfaceXMagneticField, interfaceYMagneticField, interfaceZMagneticField, material1Parameters);
}

double MHDReducedStateVector::computeMaterial1YSlowMagnetoAcousticSpeed(MHDMaterialParameters material1Parameters)
{
    return MHDWaveSpeeds::computeYSlowMagnetoAcousticSpeed(material1Density, interfacePressure, interfaceXMagneticField, interfaceYMagneticField, interfaceZMagneticField, material1Parameters);
}

double MHDReducedStateVector::computeMaterial1XFastMagnetoAcousticSpeed(MHDMaterialParameters material1Parameters)
{
    return MHDWaveSpeeds::computeXFastMagnetoAcousticSpeed(material2Density, interfacePressure, interfaceXMagneticField, interfaceYMagneticField, interfaceZMagneticField, material1Parameters);
}

double MHDReducedStateVector::computeMaterial1YFastMagnetoAcousticSpeed(MHDMaterialParameters material1Parameters)
{
    return MHDWaveSpeeds::computeYFastMagnetoAcousticSpeed(material2Density, interfacePressure, interfaceXMagneticField, interfaceYMagneticField, interfaceZMagneticField, material1Parameters);
}

double MHDReducedStateVector::computeMaterial2SpecificInternalEnergy(MHDMaterialParameters material2Parameters)
{
    return MHDEquationOfState::computeSpecificInternalEnergy(material2Density, interfacePressure, material2Parameters);
}

double MHDReducedStateVector::computeMaterial2TotalEnergy(MHDMaterialParameters material2Parameters)
{
    double velocitySquared = (interfaceXVelocity * interfaceXVelocity) + (interfaceYVelocity * interfaceYVelocity) + (interfaceZVelocity * interfaceZVelocity);
    double totalHydrodynamicEnergy = material2Density * ((0.5 * velocitySquared) + computeMaterial2SpecificInternalEnergy(material2Parameters));

    double magneticFieldSquared = (interfaceXMagneticField * interfaceXMagneticField) + (interfaceYMagneticField * interfaceYMagneticField) + (interfaceZMagneticField * interfaceZMagneticField);

    return totalHydrodynamicEnergy + (0.5 * magneticFieldSquared);
}

double MHDReducedStateVector::computeMaterial2SoundSpeed(MHDMaterialParameters material2Parameters)
{
    return MHDEquationOfState::computeSoundSpeed(material2Density, interfacePressure, material2Parameters);
}

double MHDReducedStateVector::computeMaterial2Entropy(MHDMaterialParameters material2Parameters)
{
    return MHDEquationOfState::computeEntropy(material2Density, interfacePressure, material2Parameters);
}

double MHDReducedStateVector::computeMaterial2AlfvenWaveSpeed()
{
    return MHDWaveSpeeds::computeAlfvenWaveSpeed(material2Density, interfaceXMagneticField, interfaceYMagneticField, interfaceZMagneticField);
}

double MHDReducedStateVector::computeMaterial2XSlowMagnetoAcousticSpeed(MHDMaterialParameters material2Parameters)
{
    return MHDWaveSpeeds::computeXSlowMagnetoAcousticSpeed(material2Density, interfacePressure, interfaceXMagneticField, interfaceYMagneticField, interfaceZMagneticField, material2Parameters);
}

double MHDReducedStateVector::computeMaterial2YSlowMagnetoAcousticSpeed(MHDMaterialParameters material2Parameters)
{
    return MHDWaveSpeeds::computeYSlowMagnetoAcousticSpeed(material2Density, interfacePressure, interfaceXMagneticField, interfaceYMagneticField, interfaceZMagneticField, material2Parameters);
}

double MHDReducedStateVector::computeMaterial2XFastMagnetoAcousticSpeed(MHDMaterialParameters material2Parameters)
{
    return MHDWaveSpeeds::computeXFastMagnetoAcousticSpeed(material2Density, interfacePressure, interfaceXMagneticField, interfaceYMagneticField, interfaceZMagneticField, material2Parameters);
}

double MHDReducedStateVector::computeMaterial2YFastMagnetoAcousticSpeed(MHDMaterialParameters material2Parameters)
{
    return MHDWaveSpeeds::computeYFastMagnetoAcousticSpeed(material2Density, interfacePressure, interfaceXMagneticField, interfaceYMagneticField, interfaceZMagneticField, material2Parameters);
}

double MHDReducedStateVector::computeTotalDensity()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1Density) + (material2VolumeFraction * material2Density);
}

void MHDReducedStateVector::relaxTotalDensity()
{
    double totalDensity = computeTotalDensity();

    material1Density = totalDensity;
    material2Density = totalDensity;
}

void MHDReducedStateVector::setMaterial1VolumeFraction(double newMaterial1VolumeFraction)
{
    material1VolumeFraction = newMaterial1VolumeFraction;
}

void MHDReducedStateVector::setInterfaceXVelocity(double newInterfaceXVelocity)
{
    interfaceXVelocity = newInterfaceXVelocity;
}

void MHDReducedStateVector::setInterfaceYVelocity(double newInterfaceYVelocity)
{
    interfaceYVelocity = newInterfaceYVelocity;
}

void MHDReducedStateVector::setInterfaceZVelocity(double newInterfaceZVelocity)
{
    interfaceZVelocity = newInterfaceZVelocity;
}

void MHDReducedStateVector::setMaterial1Density(double newMaterial1Density)
{
    material1Density = newMaterial1Density;
}

void MHDReducedStateVector::setMaterial2Density(double newMaterial2Density)
{
    material2Density = newMaterial2Density;
}

void MHDReducedStateVector::setInterfaceXMagneticField(double newInterfaceXMagneticField)
{
    interfaceXMagneticField = newInterfaceXMagneticField;
}

void MHDReducedStateVector::setInterfaceYMagneticField(double newInterfaceYMagneticField)
{
    interfaceYMagneticField = newInterfaceYMagneticField;
}

void MHDReducedStateVector::setInterfaceZMagneticField(double newInterfaceZMagneticField)
{
    interfaceZMagneticField = newInterfaceZMagneticField;
}

void MHDReducedStateVector::setInterfaceAuxiliaryField(double newInterfaceAuxiliaryField)
{
    interfaceAuxiliaryField = newInterfaceAuxiliaryField;
}

double MHDReducedStateVector::getMaterial1VolumeFraction()
{
    return material1VolumeFraction;
}

double MHDReducedStateVector::getInterfaceXVelocity()
{
    return interfaceXVelocity;
}

double MHDReducedStateVector::getInterfaceYVelocity()
{
    return interfaceYVelocity;
}

double MHDReducedStateVector::getInterfaceZVelocity()
{
    return interfaceZVelocity;
}

double MHDReducedStateVector::getMaterial1Density()
{
    return material1Density;
}

double MHDReducedStateVector::getMaterial2Density()
{
    return material2Density;
}

double MHDReducedStateVector::getInterfaceXMagneticField()
{
    return interfaceXMagneticField;
}

double MHDReducedStateVector::getInterfaceYMagneticField()
{
    return interfaceYMagneticField;
}

double MHDReducedStateVector::getInterfaceZMagneticField()
{
    return interfaceZMagneticField;
}

double MHDReducedStateVector::getInterfaceAuxiliaryField()
{
    return interfaceAuxiliaryField;
}

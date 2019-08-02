#include "mhdmultiphysicsstatevector.h"

MHDMultiphysicsStateVector::MHDMultiphysicsStateVector()
{
    material1VolumeFraction = 0.999;
    interfaceXVelocity = 0.0;
    interfaceYVelocity = 0.0;
    interfaceZVelocity = 0.0;

    material1Density = 1.0;
    material1Pressure = 1.0;

    material1XMagneticField = 0.0;
    material1YMagneticField = 0.0;
    material1ZMagneticField = 0.0;

    material1AuxiliaryField = 0.0;

    material2Density = 1.0;
    material2Pressure = 1.0;

    material2XMagneticField = 0.0;
    material2YMagneticField = 0.0;
    material2ZMagneticField = 0.0;

    material2AuxiliaryField = 0.0;
}

MHDMultiphysicsStateVector::MHDMultiphysicsStateVector(double newMaterial1VolumeFraction, double newInterfaceXVelocity, double newInterfaceYVelocity, double newInterfaceZVelocity,
                                                       double newMaterial1Density, double newMaterial1Pressure, double newMaterial1XMagneticField, double newMaterial1YMagneticField,
                                                       double newMaterial1ZMagneticField, double newMaterial1AuxiliaryField, double newMaterial2Density, double newMaterial2Pressure,
                                                       double newMaterial2XMagneticField, double newMaterial2YMagneticField, double newMaterial2ZMagneticField, double newMaterial2AuxiliaryField)
{
    material1VolumeFraction = newMaterial1VolumeFraction;
    interfaceXVelocity = newInterfaceXVelocity;
    interfaceYVelocity = newInterfaceYVelocity;
    interfaceZVelocity = newInterfaceZVelocity;

    material1Density = newMaterial1Density;
    material1Pressure = newMaterial1Pressure;

    material1XMagneticField = newMaterial1XMagneticField;
    material1YMagneticField = newMaterial1YMagneticField;
    material1ZMagneticField = newMaterial1ZMagneticField;

    material1AuxiliaryField = newMaterial1AuxiliaryField;

    material2Density = newMaterial2Density;
    material2Pressure = newMaterial2Pressure;

    material2XMagneticField = newMaterial2XMagneticField;
    material2YMagneticField = newMaterial2YMagneticField;
    material2ZMagneticField = newMaterial2ZMagneticField;

    material2AuxiliaryField = newMaterial2AuxiliaryField;
}

void MHDMultiphysicsStateVector::setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector)
{
    material1VolumeFraction = newPrimitiveVariableVector[0];
    interfaceXVelocity = newPrimitiveVariableVector[1];
    interfaceYVelocity = newPrimitiveVariableVector[2];
    interfaceZVelocity = newPrimitiveVariableVector[3];

    material1Density = newPrimitiveVariableVector[4];
    material1Pressure = newPrimitiveVariableVector[5];

    material1XMagneticField = newPrimitiveVariableVector[6];
    material1YMagneticField = newPrimitiveVariableVector[7];
    material1ZMagneticField = newPrimitiveVariableVector[8];

    material1AuxiliaryField = newPrimitiveVariableVector[9];

    material2Density = newPrimitiveVariableVector[10];
    material2Pressure = newPrimitiveVariableVector[11];

    material2XMagneticField = newPrimitiveVariableVector[12];
    material2YMagneticField = newPrimitiveVariableVector[13];
    material2ZMagneticField = newPrimitiveVariableVector[14];

    material2AuxiliaryField = newPrimitiveVariableVector[15];
}

void MHDMultiphysicsStateVector::setConservedVariableVector(vector<double> newConservedVariableVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
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
    material2Density = newConservedVariableVector[11] / material2VolumeFraction;

    material1XMagneticField = newConservedVariableVector[7] / pow(material1VolumeFraction, 1.0 / 2.0);
    material1YMagneticField = newConservedVariableVector[8] / pow(material1VolumeFraction, 1.0 / 2.0);
    material1ZMagneticField = newConservedVariableVector[9] / pow(material1VolumeFraction, 1.0 / 2.0);

    material2XMagneticField = newConservedVariableVector[13] / pow(material2VolumeFraction, 1.0 / 2.0);
    material2YMagneticField = newConservedVariableVector[14] / pow(material2VolumeFraction, 1.0 / 2.0);
    material2ZMagneticField = newConservedVariableVector[15] / pow(material2VolumeFraction, 1.0 / 2.0);

    double material1MagneticFieldSquared = (material1XMagneticField * material1XMagneticField) + (material1YMagneticField * material1YMagneticField) +
            (material1ZMagneticField * material1ZMagneticField);
    double material2MagneticFieldSquared = (material2XMagneticField * material2XMagneticField) + (material2YMagneticField * material2YMagneticField) +
            (material2ZMagneticField * material2ZMagneticField);
    double velocitySquared = (interfaceXVelocity * interfaceXVelocity) + (interfaceYVelocity * interfaceYVelocity) + (interfaceZVelocity * interfaceZVelocity);

    double material1TotalEnergy = ((newConservedVariableVector[6] / material1VolumeFraction) - (0.5 * material1MagneticFieldSquared)) / material1Density;
    double material1SpecificInternalEnergy = material1TotalEnergy - (0.5 * velocitySquared);
    material1Pressure = MHDEquationOfState::computePressure(material1Density, material1SpecificInternalEnergy, material1Parameters);

    double material2TotalEnergy = ((newConservedVariableVector[12] / material2VolumeFraction) - (0.5 * material2MagneticFieldSquared)) / material2Density;
    double material2SpecificInternalEnergy = material2TotalEnergy - (0.5 * velocitySquared);
    material2Pressure = MHDEquationOfState::computePressure(material2Density, material2SpecificInternalEnergy, material2Parameters);

    material1AuxiliaryField = newConservedVariableVector[10] / material1VolumeFraction;
    material2AuxiliaryField = newConservedVariableVector[16] / material2VolumeFraction;
}

vector<double> MHDMultiphysicsStateVector::computePrimitiveVariableVector()
{
    vector<double> primitiveVariableVector(16);

    primitiveVariableVector[0] = material1VolumeFraction;
    primitiveVariableVector[1] = interfaceXVelocity;
    primitiveVariableVector[2] = interfaceYVelocity;
    primitiveVariableVector[3] = interfaceZVelocity;

    primitiveVariableVector[4] = material1Density;
    primitiveVariableVector[5] = material1Pressure;

    primitiveVariableVector[6] = material1XMagneticField;
    primitiveVariableVector[7] = material1YMagneticField;
    primitiveVariableVector[8] = material1ZMagneticField;

    primitiveVariableVector[9] = material1AuxiliaryField;

    primitiveVariableVector[10] = material2Density;
    primitiveVariableVector[11] = material2Pressure;

    primitiveVariableVector[12] = material2XMagneticField;
    primitiveVariableVector[13] = material2YMagneticField;
    primitiveVariableVector[14] = material2ZMagneticField;

    primitiveVariableVector[15] = material2AuxiliaryField;

    return primitiveVariableVector;
}

vector<double> MHDMultiphysicsStateVector::computeConservedVariableVector(MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> conservedVariableVector(17);

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
    conservedVariableVector[11] = material2VolumeFraction * material2Density;

    double material1MagneticFieldSquared = (material1XMagneticField * material1XMagneticField) + (material1YMagneticField * material1YMagneticField) +
            (material1ZMagneticField * material1ZMagneticField);
    double material2MagneticFieldSquared = (material2XMagneticField * material2XMagneticField) + (material2YMagneticField * material2YMagneticField) +
            (material2ZMagneticField * material2ZMagneticField);
    double velocitySquared = (interfaceXVelocity * interfaceXVelocity) + (interfaceYVelocity * interfaceYVelocity) + (interfaceZVelocity * interfaceZVelocity);

    double material1TotalHydrodynamicEnergy = material1Density * ((0.5 * velocitySquared) + MHDEquationOfState::computeSpecificInternalEnergy(material1Density, material1Pressure,
                                                                                                                                              material1Parameters));
    double material2TotalHydrodynamicEnergy = material2Density * ((0.5 * velocitySquared) + MHDEquationOfState::computeSpecificInternalEnergy(material2Density, material2Pressure,
                                                                                                                                              material2Parameters));

    double material1TotalEnergy = material1TotalHydrodynamicEnergy + (0.5 * material1MagneticFieldSquared);
    double material2TotalEnergy = material2TotalHydrodynamicEnergy + (0.5 * material2MagneticFieldSquared);

    conservedVariableVector[6] = material1VolumeFraction * material1TotalEnergy;
    conservedVariableVector[12] = material2VolumeFraction * material2TotalEnergy;

    conservedVariableVector[7] = pow(material1VolumeFraction, 1.0 / 2.0) * material1XMagneticField;
    conservedVariableVector[8] = pow(material1VolumeFraction, 1.0 / 2.0) * material1YMagneticField;
    conservedVariableVector[9] = pow(material1VolumeFraction, 1.0 / 2.0) * material1ZMagneticField;

    conservedVariableVector[10] = material1VolumeFraction * material1AuxiliaryField;

    conservedVariableVector[13] = pow(material2VolumeFraction, 1.0 / 2.0) * material2XMagneticField;
    conservedVariableVector[14] = pow(material2VolumeFraction, 1.0 / 2.0) * material2YMagneticField;
    conservedVariableVector[15] = pow(material2VolumeFraction, 1.0 / 2.0) * material2ZMagneticField;

    conservedVariableVector[16] = material2VolumeFraction * material2AuxiliaryField;

    return conservedVariableVector;
}

vector<double> MHDMultiphysicsStateVector::computeXFluxVector(vector<double> conservedVariableVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> fluxVector(17);

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
    double computedMaterial2Density = conservedVariableVector[11] / computedMaterial2VolumeFraction;

    double computedMaterial1XMagneticField = conservedVariableVector[7] / pow(computedMaterial1VolumeFraction, 1.0 / 2.0);
    double computedMaterial1YMagneticField = conservedVariableVector[8] / pow(computedMaterial1VolumeFraction, 1.0 / 2.0);
    double computedMaterial1ZMagneticField = conservedVariableVector[9] / pow(computedMaterial1VolumeFraction, 1.0 / 2.0);

    double computedMaterial2XMagneticField = conservedVariableVector[13] / pow(computedMaterial2VolumeFraction, 1.0 / 2.0);
    double computedMaterial2YMagneticField = conservedVariableVector[14] / pow(computedMaterial2VolumeFraction, 1.0 / 2.0);
    double computedMaterial2ZMagneticField = conservedVariableVector[15] / pow(computedMaterial2VolumeFraction, 1.0 / 2.0);

    /*
    double computedInterfaceXMagneticField = (computedMaterial1VolumeFraction * computedMaterial1XMagneticField) + (computedMaterial2VolumeFraction * computedMaterial2XMagneticField);
    double computedInterfaceYMagneticField = (computedMaterial1VolumeFraction * computedMaterial1YMagneticField) + (computedMaterial2VolumeFraction * computedMaterial2YMagneticField);
    double computedInterfaceZMagneticField = (computedMaterial1VolumeFraction * computedMaterial1ZMagneticField) + (computedMaterial2VolumeFraction * computedMaterial2ZMagneticField);
    */
    double computedInterfaceXMagneticFieldSquared = (computedMaterial1VolumeFraction * (computedMaterial1XMagneticField * computedMaterial1XMagneticField)) +
            (computedMaterial2VolumeFraction * (computedMaterial2XMagneticField * computedMaterial2XMagneticField));
    double computedInterfaceYMagneticFieldSquared = (computedMaterial1VolumeFraction * (computedMaterial1YMagneticField * computedMaterial1YMagneticField)) +
            (computedMaterial2VolumeFraction * (computedMaterial2YMagneticField * computedMaterial2YMagneticField));
    double computedInterfaceZMagneticFieldSquared = (computedMaterial1VolumeFraction * (computedMaterial1ZMagneticField * computedMaterial1ZMagneticField)) +
            (computedMaterial2VolumeFraction * (computedMaterial2ZMagneticField * computedMaterial2ZMagneticField));

    double computedInterfaceXMagneticField = pow(computedInterfaceXMagneticFieldSquared, 1.0 / 2.0);
    double computedInterfaceYMagneticField = pow(computedInterfaceYMagneticFieldSquared, 1.0 / 2.0);
    double computedInterfaceZMagneticField = pow(computedInterfaceZMagneticFieldSquared, 1.0 / 2.0);

    double computedMaterial1AuxiliaryField = conservedVariableVector[10] / computedMaterial1VolumeFraction;
    double computedMaterial2AuxiliaryField = conservedVariableVector[16] / computedMaterial2VolumeFraction;

    double velocitySquared = (computedInterfaceXVelocity * computedInterfaceXVelocity) + (computedInterfaceYVelocity *  computedInterfaceYVelocity) +
            (computedInterfaceZVelocity * computedInterfaceZVelocity);

    double material1MagneticFieldSquared = (computedMaterial1XMagneticField * computedMaterial1XMagneticField) + (computedMaterial1YMagneticField * computedMaterial1YMagneticField) +
            (computedMaterial1ZMagneticField * computedMaterial1ZMagneticField);
    double material2MagneticFieldSquared = (computedMaterial2XMagneticField * computedMaterial2XMagneticField) + (computedMaterial2YMagneticField * computedMaterial2YMagneticField) +
            (computedMaterial2ZMagneticField * computedMaterial2ZMagneticField);

    /*
    double interfaceMagneticFieldSquared = (computedInterfaceXMagneticField * computedInterfaceXMagneticField) + (computedInterfaceYMagneticField * computedInterfaceYMagneticField) +
            (computedInterfaceZMagneticField * computedInterfaceZMagneticField);
    */
    double interfaceMagneticFieldSquared = computedInterfaceXMagneticFieldSquared + computedInterfaceYMagneticFieldSquared + computedInterfaceZMagneticFieldSquared;

    double material1TotalEnergy = conservedVariableVector[6] / computedMaterial1VolumeFraction;
    double material1TotalHydrodynamicEnergy = (material1TotalEnergy - (0.5 * material1MagneticFieldSquared)) / computedMaterial1Density;
    double material1SpecificInternalEnergy = material1TotalHydrodynamicEnergy - (0.5 * velocitySquared);

    double material2TotalEnergy = conservedVariableVector[12] / computedMaterial2VolumeFraction;
    double material2TotalHydrodynamicEnergy = (material2TotalEnergy - (0.5 * material2MagneticFieldSquared)) / computedMaterial2Density;
    double material2SpecificInternalEnergy = material2TotalHydrodynamicEnergy - (0.5 * velocitySquared);

    double computedMaterial1Pressure = MHDEquationOfState::computePressure(computedMaterial1Density, material1SpecificInternalEnergy, material1Parameters);
    double computedMaterial2Pressure = MHDEquationOfState::computePressure(computedMaterial2Density, material2SpecificInternalEnergy, material2Parameters);

    double material1HyperbolicWaveSpeedSquared = material1Parameters.computeHyperbolicWaveSpeedSquared();
    double material2HyperbolicWaveSpeedSquared = material2Parameters.computeHyperbolicWaveSpeedSquared();

    fluxVector[0] = computedTotalDensity * computedInterfaceXVelocity;
    fluxVector[1] = computedTotalDensity * (computedInterfaceXVelocity * computedMaterial1VolumeFraction);

    /*
    fluxVector[2] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceXVelocity)) + ((computedMaterial1VolumeFraction * computedMaterial1Pressure) +
                                                                                                          (computedMaterial2VolumeFraction * computedMaterial2Pressure)) +
            (0.5 * interfaceMagneticFieldSquared) - (computedInterfaceXMagneticField * computedInterfaceXMagneticField);
    */
    fluxVector[2] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceXVelocity)) + ((computedMaterial1VolumeFraction * computedMaterial1Pressure) +
                                                                                                          (computedMaterial2VolumeFraction * computedMaterial2Pressure)) +
            (0.5 * interfaceMagneticFieldSquared) - computedInterfaceXMagneticFieldSquared;
    fluxVector[3] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceYVelocity)) - (computedInterfaceXMagneticField * computedInterfaceYMagneticField);
    fluxVector[4] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceZVelocity)) - (computedInterfaceXMagneticField * computedInterfaceZMagneticField);

    fluxVector[5] = computedMaterial1VolumeFraction * (computedMaterial1Density * computedInterfaceXVelocity);
    fluxVector[11] = computedMaterial2VolumeFraction * (computedMaterial2Density * computedInterfaceXVelocity);

    vector<double> velocityVector(3);
    velocityVector[0] = computedInterfaceXVelocity;
    velocityVector[1] = computedInterfaceYVelocity;
    velocityVector[2] = computedInterfaceZVelocity;

    vector<double> material1MagneticFieldVector(3);
    material1MagneticFieldVector[0] = computedMaterial1XMagneticField;
    material1MagneticFieldVector[1] = computedMaterial1YMagneticField;
    material1MagneticFieldVector[2] = computedMaterial1ZMagneticField;

    double material1MagneticFieldVelocityVectorProduct = VectorAlgebra::computeDotProduct(velocityVector, material1MagneticFieldVector);
    fluxVector[6] = (computedMaterial1VolumeFraction * ((material1TotalEnergy + computedMaterial1Pressure + (0.5 * material1MagneticFieldSquared)) * computedInterfaceXVelocity)) -
            (computedMaterial1VolumeFraction * (material1MagneticFieldVelocityVectorProduct * computedMaterial1XMagneticField));

    vector<double> material2MagneticFieldVector(3);
    material2MagneticFieldVector[0] = computedMaterial2XMagneticField;
    material2MagneticFieldVector[1] = computedMaterial2YMagneticField;
    material2MagneticFieldVector[2] = computedMaterial2ZMagneticField;

    double material2MagneticFieldVelocityVectorProduct = VectorAlgebra::computeDotProduct(velocityVector, material2MagneticFieldVector);
    fluxVector[12] = (computedMaterial2VolumeFraction * ((material2TotalEnergy + computedMaterial2Pressure + (0.5 * material2MagneticFieldSquared)) * computedInterfaceXVelocity)) -
            (computedMaterial2VolumeFraction * (material2MagneticFieldVelocityVectorProduct * computedMaterial2XMagneticField));

    fluxVector[7] = pow(computedMaterial1VolumeFraction, 1.0 / 2.0) * computedMaterial1AuxiliaryField;
    fluxVector[8] = pow(computedMaterial1VolumeFraction, 1.0 / 2.0) * ((computedMaterial1YMagneticField * computedInterfaceXVelocity) - (computedMaterial1XMagneticField * computedInterfaceYVelocity));
    fluxVector[9] = pow(computedMaterial1VolumeFraction, 1.0 / 2.0) * ((computedMaterial1ZMagneticField * computedInterfaceXVelocity) - (computedMaterial1XMagneticField * computedInterfaceZVelocity));

    fluxVector[10] = computedMaterial1VolumeFraction * (material1HyperbolicWaveSpeedSquared * computedMaterial1XMagneticField);

    fluxVector[13] = pow(computedMaterial2VolumeFraction, 1.0 / 2.0) * computedMaterial2AuxiliaryField;
    fluxVector[14] = pow(computedMaterial2VolumeFraction, 1.0 / 2.0) * ((computedMaterial2YMagneticField * computedInterfaceXVelocity) - (computedMaterial2XMagneticField * computedInterfaceYVelocity));
    fluxVector[15] = pow(computedMaterial2VolumeFraction, 1.0 / 2.0) * ((computedMaterial2ZMagneticField * computedInterfaceXVelocity) - (computedMaterial2XMagneticField * computedInterfaceZVelocity));

    fluxVector[16] = computedMaterial2VolumeFraction * (material2HyperbolicWaveSpeedSquared * computedMaterial2XMagneticField);

    return fluxVector;
}

vector<double> MHDMultiphysicsStateVector::computeXFluxVector(MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    return computeXFluxVector(computeConservedVariableVector(material1Parameters, material2Parameters), material1Parameters, material2Parameters);
}

double MHDMultiphysicsStateVector::computeMaterial1SpecificInternalEnergy(MHDMaterialParameters material1Parameters)
{
    return MHDEquationOfState::computeSpecificInternalEnergy(material1Density, material1Pressure, material1Parameters);
}

double MHDMultiphysicsStateVector::computeMaterial1TotalEnergy(MHDMaterialParameters material1Parameters)
{
    double velocitySquared = (interfaceXVelocity * interfaceXVelocity) + (interfaceYVelocity * interfaceYVelocity) + (interfaceZVelocity * interfaceZVelocity);
    double totalHydrodynamicEnergy = material1Density * ((0.5 * velocitySquared) + computeMaterial1SpecificInternalEnergy(material1Parameters));

    double magneticFieldSquared = (material1XMagneticField * material1XMagneticField) + (material1YMagneticField * material1YMagneticField) + (material1ZMagneticField * material1ZMagneticField);

    return totalHydrodynamicEnergy + (0.5 * magneticFieldSquared);
}

double MHDMultiphysicsStateVector::computeMaterial1SoundSpeed(MHDMaterialParameters material1Parameters)
{
    return MHDEquationOfState::computeSoundSpeed(material1Density, material1Pressure, material1Parameters);
}

double MHDMultiphysicsStateVector::computeMaterial1Entropy(MHDMaterialParameters material1Parameters)
{
    return MHDEquationOfState::computeEntropy(material1Density, material1Pressure, material1Parameters);
}

double MHDMultiphysicsStateVector::computeMaterial1AlfvenWaveSpeed()
{
    return MHDWaveSpeeds::computeAlfvenWaveSpeed(material1Density, material1XMagneticField, material1YMagneticField, material1ZMagneticField);
}

double MHDMultiphysicsStateVector::computeMaterial1XSlowMagnetoAcousticSpeed(MHDMaterialParameters material1Parameters)
{
    return MHDWaveSpeeds::computeXSlowMagnetoAcousticSpeed(material1Density, material1Pressure, material1XMagneticField, material1YMagneticField, material1ZMagneticField, material1Parameters);
}

double MHDMultiphysicsStateVector::computeMaterial1XFastMagnetoAcousticSpeed(MHDMaterialParameters material1Parameters)
{
    return MHDWaveSpeeds::computeXFastMagnetoAcousticSpeed(material1Density, material1Pressure, material1XMagneticField, material1YMagneticField, material1ZMagneticField, material1Parameters);
}

double MHDMultiphysicsStateVector::computeMaterial2SpecificInternalEnergy(MHDMaterialParameters material2Parameters)
{
    return MHDEquationOfState::computeSpecificInternalEnergy(material2Density, material2Pressure, material2Parameters);
}

double MHDMultiphysicsStateVector::computeMaterial2TotalEnergy(MHDMaterialParameters material2Parameters)
{
    double velocitySquared = (interfaceXVelocity * interfaceXVelocity) + (interfaceYVelocity * interfaceYVelocity) + (interfaceZVelocity * interfaceZVelocity);
    double totalHydrodynamicEnergy = material2Density * ((0.5 * velocitySquared) + computeMaterial2SpecificInternalEnergy(material2Parameters));

    double magneticFieldSquared = (material2XMagneticField * material2XMagneticField) + (material2YMagneticField * material2YMagneticField) + (material2ZMagneticField * material2ZMagneticField);

    return totalHydrodynamicEnergy - (0.5 * magneticFieldSquared);
}

double MHDMultiphysicsStateVector::computeMaterial2SoundSpeed(MHDMaterialParameters material2Parameters)
{
    return MHDEquationOfState::computeSoundSpeed(material2Density, material2Pressure, material2Parameters);
}

double MHDMultiphysicsStateVector::computeMaterial2Entropy(MHDMaterialParameters material2Parameters)
{
    return MHDEquationOfState::computeEntropy(material2Density, material2Pressure, material2Parameters);
}

double MHDMultiphysicsStateVector::computeMaterial2AlfvenWaveSpeed()
{
    return MHDWaveSpeeds::computeAlfvenWaveSpeed(material2Density, material2XMagneticField, material2YMagneticField, material2ZMagneticField);
}

double MHDMultiphysicsStateVector::computeMaterial2XSlowMagnetoAcousticSpeed(MHDMaterialParameters material2Parameters)
{
    return MHDWaveSpeeds::computeXSlowMagnetoAcousticSpeed(material2Density, material2Pressure, material2XMagneticField, material2YMagneticField, material2ZMagneticField, material2Parameters);
}

double MHDMultiphysicsStateVector::computeMaterial2XFastMagnetoAcousticSpeed(MHDMaterialParameters material2Parameters)
{
    return MHDWaveSpeeds::computeXFastMagnetoAcousticSpeed(material2Density, material2Pressure, material2XMagneticField, material2YMagneticField, material2ZMagneticField, material2Parameters);
}

double MHDMultiphysicsStateVector::computeTotalDensity()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1Density) + (material2VolumeFraction * material2Density);
}

double MHDMultiphysicsStateVector::computeTotalPressure()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1Pressure) + (material2VolumeFraction * material2Pressure);
}

double MHDMultiphysicsStateVector::computeTotalXMagneticField()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1XMagneticField) + (material2VolumeFraction * material2XMagneticField);
}

double MHDMultiphysicsStateVector::computeTotalYMagneticField()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1YMagneticField) + (material2VolumeFraction * material2YMagneticField);
}

double MHDMultiphysicsStateVector::computeTotalZMagneticField()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1ZMagneticField) + (material2VolumeFraction * material2ZMagneticField);
}

double MHDMultiphysicsStateVector::computeTotalAuxiliaryField()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1AuxiliaryField) + (material2VolumeFraction * material2AuxiliaryField);
}

void MHDMultiphysicsStateVector::relaxTotalDensity()
{
    double totalDensity = computeTotalDensity();

    material1Density = totalDensity;
    material2Density = totalDensity;
}

void MHDMultiphysicsStateVector::relaxTotalPressure()
{
    double totalPressure = computeTotalPressure();

    material1Pressure = totalPressure;
    material2Pressure = totalPressure;
}

void MHDMultiphysicsStateVector::relaxTotalXMagneticField()
{
    double totalXMagneticField = computeTotalXMagneticField();

    material1XMagneticField = totalXMagneticField;
    material2XMagneticField = totalXMagneticField;
}

void MHDMultiphysicsStateVector::relaxTotalYMagneticField()
{
    double totalYMagneticField = computeTotalYMagneticField();

    material1YMagneticField = totalYMagneticField;
    material2YMagneticField = totalYMagneticField;
}

void MHDMultiphysicsStateVector::relaxTotalZMagneticField()
{
    double totalZMagneticField = computeTotalZMagneticField();

    material1ZMagneticField = totalZMagneticField;
    material2ZMagneticField = totalZMagneticField;
}

void MHDMultiphysicsStateVector::relaxTotalAuxiliaryField()
{
    double totalAuxiliaryField = computeTotalAuxiliaryField();

    material1AuxiliaryField = totalAuxiliaryField;
    material2AuxiliaryField = totalAuxiliaryField;
}

void MHDMultiphysicsStateVector::setMaterial1VolumeFraction(double newMaterial1VolumeFraction)
{
    material1VolumeFraction = newMaterial1VolumeFraction;
}

void MHDMultiphysicsStateVector::setInterfaceXVelocity(double newInterfaceXVelocity)
{
    interfaceXVelocity = newInterfaceXVelocity;
}

void MHDMultiphysicsStateVector::setInterfaceYVelocity(double newInterfaceYVelocity)
{
    interfaceYVelocity = newInterfaceYVelocity;
}

void MHDMultiphysicsStateVector::setInterfaceZVelocity(double newInterfaceZVelocity)
{
    interfaceZVelocity = newInterfaceZVelocity;
}

void MHDMultiphysicsStateVector::setMaterial1Density(double newMaterial1Density)
{
    material1Density = newMaterial1Density;
}

void MHDMultiphysicsStateVector::setMaterial1Pressure(double newMaterial1Pressure)
{
    material1Pressure = newMaterial1Pressure;
}

void MHDMultiphysicsStateVector::setMaterial1XMagneticField(double newMaterial1XMagneticField)
{
    material1XMagneticField = newMaterial1XMagneticField;
}

void MHDMultiphysicsStateVector::setMaterial1YMagneticField(double newMaterial1YMagneticField)
{
    material1YMagneticField = newMaterial1YMagneticField;
}

void MHDMultiphysicsStateVector::setMaterial1ZMagneticField(double newMaterial1ZMagneticField)
{
    material1ZMagneticField = newMaterial1ZMagneticField;
}

void MHDMultiphysicsStateVector::setMaterial1AuxiliaryField(double newMaterial1AuxiliaryField)
{
    material1AuxiliaryField = newMaterial1AuxiliaryField;
}

void MHDMultiphysicsStateVector::setMaterial2Density(double newMaterial2Density)
{
    material2Density = newMaterial2Density;
}

void MHDMultiphysicsStateVector::setMaterial2Pressure(double newMaterial2Pressure)
{
    material2Pressure = newMaterial2Pressure;
}

void MHDMultiphysicsStateVector::setMaterial2XMagneticField(double newMaterial2XMagneticField)
{
    material2XMagneticField = newMaterial2XMagneticField;
}

void MHDMultiphysicsStateVector::setMaterial2YMagneticField(double newMaterial2YMagneticField)
{
    material2YMagneticField = newMaterial2YMagneticField;
}

void MHDMultiphysicsStateVector::setMaterial2ZMagneticField(double newMaterial2ZMagneticField)
{
    material2ZMagneticField = newMaterial2ZMagneticField;
}

void MHDMultiphysicsStateVector::setMaterial2AuxiliaryField(double newMaterial2AuxiliaryField)
{
    material2AuxiliaryField = newMaterial2AuxiliaryField;
}

double MHDMultiphysicsStateVector::getMaterial1VolumeFraction()
{
    return material1VolumeFraction;
}

double MHDMultiphysicsStateVector::getInterfaceXVelocity()
{
    return interfaceXVelocity;
}

double MHDMultiphysicsStateVector::getInterfaceYVelocity()
{
    return interfaceYVelocity;
}

double MHDMultiphysicsStateVector::getInterfaceZVelocity()
{
    return interfaceZVelocity;
}

double MHDMultiphysicsStateVector::getMaterial1Density()
{
    return material1Density;
}

double MHDMultiphysicsStateVector::getMaterial1Pressure()
{
    return material1Pressure;
}

double MHDMultiphysicsStateVector::getMaterial1XMagneticField()
{
    return material1XMagneticField;
}

double MHDMultiphysicsStateVector::getMaterial1YMagneticField()
{
    return material1YMagneticField;
}

double MHDMultiphysicsStateVector::getMaterial1ZMagneticField()
{
    return material1ZMagneticField;
}

double MHDMultiphysicsStateVector::getMaterial1AuxiliaryField()
{
    return material1AuxiliaryField;
}

double MHDMultiphysicsStateVector::getMaterial2Density()
{
    return material2Density;
}

double MHDMultiphysicsStateVector::getMaterial2Pressure()
{
    return material2Pressure;
}

double MHDMultiphysicsStateVector::getMaterial2XMagneticField()
{
    return material2XMagneticField;
}

double MHDMultiphysicsStateVector::getMaterial2YMagneticField()
{
    return material2YMagneticField;
}

double MHDMultiphysicsStateVector::getMaterial2ZMagneticField()
{
    return material2ZMagneticField;
}

double MHDMultiphysicsStateVector::getMaterial2AuxiliaryField()
{
    return material2AuxiliaryField;
}

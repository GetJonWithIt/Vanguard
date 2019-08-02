#include "mhdintermediatestatevector.h"

MHDIntermediateStateVector::MHDIntermediateStateVector()
{
    material1VolumeFraction = 0.999;
    interfaceXVelocity = 0.0;
    interfaceYVelocity = 0.0;
    interfaceZVelocity = 0.0;

    material1Density = 1.0;
    material1Pressure = 1.0;

    material2Density = 1.0;
    material2Pressure = 1.0;

    interfaceXMagneticField = 0.0;
    interfaceYMagneticField = 0.0;
    interfaceZMagneticField = 0.0;

    interfaceAuxiliaryField = 0.0;
}

MHDIntermediateStateVector::MHDIntermediateStateVector(double newMaterial1VolumeFraction, double newInterfaceXVelocity, double newInterfaceYVelocity, double newInterfaceZVelocity,
                                                       double newMaterial1Density, double newMaterial1Pressure, double newMaterial2Density, double newMaterial2Pressure,
                                                       double newInterfaceXMagneticField, double newInterfaceYMagneticField, double newInterfaceZMagneticField, double newInterfaceAuxiliaryField)
{
    material1VolumeFraction = newMaterial1VolumeFraction;
    interfaceXVelocity = newInterfaceXVelocity;
    interfaceYVelocity = newInterfaceYVelocity;
    interfaceZVelocity = newInterfaceZVelocity;

    material1Density = newMaterial1Density;
    material1Pressure = newMaterial1Pressure;

    material2Density = newMaterial2Density;
    material2Pressure = newMaterial2Pressure;

    interfaceXMagneticField = newInterfaceXMagneticField;
    interfaceYMagneticField = newInterfaceYMagneticField;
    interfaceZMagneticField = newInterfaceZMagneticField;

    interfaceAuxiliaryField = newInterfaceAuxiliaryField;
}

void MHDIntermediateStateVector::setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector)
{
    material1VolumeFraction = newPrimitiveVariableVector[0];
    interfaceXVelocity = newPrimitiveVariableVector[1];
    interfaceYVelocity = newPrimitiveVariableVector[2];
    interfaceZVelocity = newPrimitiveVariableVector[3];

    material1Density = newPrimitiveVariableVector[4];
    material1Pressure = newPrimitiveVariableVector[5];

    material2Density = newPrimitiveVariableVector[6];
    material2Pressure = newPrimitiveVariableVector[7];

    interfaceXMagneticField = newPrimitiveVariableVector[8];
    interfaceYMagneticField = newPrimitiveVariableVector[9];
    interfaceZMagneticField = newPrimitiveVariableVector[10];

    interfaceAuxiliaryField = newPrimitiveVariableVector[11];
}

void MHDIntermediateStateVector::setConservedVariableVector(vector<double> newConservedVariableVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
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
    material2Density = newConservedVariableVector[7] / material2VolumeFraction;

    interfaceXMagneticField = newConservedVariableVector[9];
    interfaceYMagneticField = newConservedVariableVector[10];
    interfaceZMagneticField = newConservedVariableVector[11];

    double magneticFieldSquared = (interfaceXMagneticField * interfaceXMagneticField) + (interfaceYMagneticField * interfaceYMagneticField) + (interfaceZMagneticField * interfaceZMagneticField);
    double velocitySquared = (interfaceXVelocity * interfaceXVelocity) + (interfaceYVelocity * interfaceYVelocity) + (interfaceZVelocity * interfaceZVelocity);

    double material1TotalEnergy = ((newConservedVariableVector[6] / material1VolumeFraction) - (0.5 * magneticFieldSquared)) / material1Density;
    double material1SpecificInternalEnergy = material1TotalEnergy - (0.5 * velocitySquared);
    material1Pressure = MHDEquationOfState::computePressure(material1Density, material1SpecificInternalEnergy, material1Parameters);

    double material2TotalEnergy = ((newConservedVariableVector[8] / material2VolumeFraction) - (0.5 * magneticFieldSquared)) / material2Density;
    double material2SpecificInternalEnergy = material2TotalEnergy - (0.5 * velocitySquared);
    material2Pressure = MHDEquationOfState::computePressure(material2Density, material2SpecificInternalEnergy, material2Parameters);

    interfaceAuxiliaryField = newConservedVariableVector[12];
}

vector<double> MHDIntermediateStateVector::computePrimitiveVariableVector()
{
    vector<double> primitiveVariableVector(12);

    primitiveVariableVector[0] = material1VolumeFraction;
    primitiveVariableVector[1] = interfaceXVelocity;
    primitiveVariableVector[2] = interfaceYVelocity;
    primitiveVariableVector[3] = interfaceZVelocity;

    primitiveVariableVector[4] = material1Density;
    primitiveVariableVector[5] = material1Pressure;

    primitiveVariableVector[6] = material2Density;
    primitiveVariableVector[7] = material2Pressure;

    primitiveVariableVector[8] = interfaceXMagneticField;
    primitiveVariableVector[9] = interfaceYMagneticField;
    primitiveVariableVector[10] = interfaceZMagneticField;

    primitiveVariableVector[11] = interfaceAuxiliaryField;

    return primitiveVariableVector;
}

vector<double> MHDIntermediateStateVector::computeConservedVariableVector(MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> conservedVariableVector(13);

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
    conservedVariableVector[7] = material2VolumeFraction * material2Density;

    double magneticFieldSquared = (interfaceXMagneticField * interfaceXMagneticField) + (interfaceYMagneticField * interfaceYMagneticField) + (interfaceZMagneticField * interfaceZMagneticField);
    double velocitySquared = (interfaceXVelocity * interfaceXVelocity) + (interfaceYVelocity * interfaceYVelocity) + (interfaceZVelocity * interfaceZVelocity);

    double material1TotalHydrodynamicEnergy = material1Density * ((0.5 * velocitySquared) + MHDEquationOfState::computeSpecificInternalEnergy(material1Density, material1Pressure, material1Parameters));
    double material2TotalHydrodynamicEnergy = material2Density * ((0.5 * velocitySquared) + MHDEquationOfState::computeSpecificInternalEnergy(material2Density, material2Pressure, material2Parameters));

    double material1TotalEnergy = material1TotalHydrodynamicEnergy + (0.5 * magneticFieldSquared);
    double material2TotalEnergy = material2TotalHydrodynamicEnergy + (0.5 * magneticFieldSquared);

    conservedVariableVector[6] = material1VolumeFraction * material1TotalEnergy;
    conservedVariableVector[8] = material2VolumeFraction * material2TotalEnergy;

    conservedVariableVector[9] = interfaceXMagneticField;
    conservedVariableVector[10] = interfaceYMagneticField;
    conservedVariableVector[11] = interfaceZMagneticField;

    conservedVariableVector[12] = interfaceAuxiliaryField;

    return conservedVariableVector;
}

vector<double> MHDIntermediateStateVector::computeXFluxVector(vector<double> conservedVariableVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> fluxVector(13);

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
    double computedMaterial2Density = conservedVariableVector[7] / computedMaterial2VolumeFraction;

    double computedInterfaceXMagneticField = conservedVariableVector[9];
    double computedInterfaceYMagneticField = conservedVariableVector[10];
    double computedInterfaceZMagneticField = conservedVariableVector[11];

    double computedInterfaceAuxiliaryField = conservedVariableVector[12];

    double velocitySquared = (computedInterfaceXVelocity * computedInterfaceXVelocity) + (computedInterfaceYVelocity * computedInterfaceYVelocity) +
            (computedInterfaceZVelocity * computedInterfaceZVelocity);
    double magneticFieldSquared = (computedInterfaceXMagneticField * computedInterfaceXMagneticField) + (computedInterfaceYMagneticField * computedInterfaceYMagneticField) +
            (computedInterfaceZMagneticField * computedInterfaceZMagneticField);

    double material1TotalEnergy = conservedVariableVector[6] / computedMaterial1VolumeFraction;
    double material1TotalHydrodynamicEnergy = (material1TotalEnergy - (0.5 * magneticFieldSquared)) / computedMaterial1Density;
    double material1SpecificInternalEnergy = material1TotalHydrodynamicEnergy - (0.5 * velocitySquared);

    double material2TotalEnergy = conservedVariableVector[8] / computedMaterial2VolumeFraction;
    double material2TotalHydrodynamicEnergy = (material2TotalEnergy - (0.5 * magneticFieldSquared)) / computedMaterial2Density;
    double material2SpecificInternalEnergy = material2TotalHydrodynamicEnergy - (0.5 * velocitySquared);

    double computedMaterial1Pressure = MHDEquationOfState::computePressure(computedMaterial1Density, material1SpecificInternalEnergy, material1Parameters);
    double computedMaterial2Pressure = MHDEquationOfState::computePressure(computedMaterial2Density, material2SpecificInternalEnergy, material2Parameters);

    double material1HyperbolicWaveSpeedSquared = material1Parameters.computeHyperbolicWaveSpeedSquared();
    double material2HyperbolicWaveSpeedSquared = material2Parameters.computeHyperbolicWaveSpeedSquared();
    double hyperbolicWaveSpeedSquared = (computedMaterial1VolumeFraction * material1HyperbolicWaveSpeedSquared) + (computedMaterial2VolumeFraction * material2HyperbolicWaveSpeedSquared);

    fluxVector[0] = computedTotalDensity * computedInterfaceXVelocity;
    fluxVector[1] = computedTotalDensity * (computedInterfaceXVelocity * computedMaterial1VolumeFraction);

    fluxVector[2] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceXVelocity)) + ((computedMaterial1VolumeFraction * computedMaterial1Pressure) +
                                                                                                          (computedMaterial2VolumeFraction * computedMaterial2Pressure)) +
            (0.5 * magneticFieldSquared) - (computedInterfaceXMagneticField * computedInterfaceXMagneticField);
    fluxVector[3] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceYVelocity)) - (computedInterfaceXMagneticField * computedInterfaceYMagneticField);
    fluxVector[4] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceZVelocity)) - (computedInterfaceXMagneticField * computedInterfaceZMagneticField);

    fluxVector[5] = computedMaterial1VolumeFraction * (computedMaterial1Density * computedInterfaceXVelocity);
    fluxVector[7] = computedMaterial2VolumeFraction * (computedMaterial2Density * computedInterfaceXVelocity);

    vector<double> velocityVector(3);
    velocityVector[0] = computedInterfaceXVelocity;
    velocityVector[1] = computedInterfaceYVelocity;
    velocityVector[2] = computedInterfaceZVelocity;

    vector<double> magneticFieldVector(3);
    magneticFieldVector[0] = computedInterfaceXMagneticField;
    magneticFieldVector[1] = computedInterfaceYMagneticField;
    magneticFieldVector[2] = computedInterfaceZMagneticField;

    double magneticFieldVelocityVectorProduct = VectorAlgebra::computeDotProduct(velocityVector, magneticFieldVector);
    fluxVector[6] = (computedMaterial1VolumeFraction * ((material1TotalEnergy + computedMaterial1Pressure + (0.5 * magneticFieldSquared)) * computedInterfaceXVelocity)) -
            (computedMaterial1VolumeFraction * (magneticFieldVelocityVectorProduct * computedInterfaceXMagneticField));

    fluxVector[8] = (computedMaterial2VolumeFraction * ((material2TotalEnergy + computedMaterial2Pressure + (0.5 * magneticFieldSquared)) * computedInterfaceXVelocity)) -
            (computedMaterial2VolumeFraction * (magneticFieldVelocityVectorProduct * computedInterfaceXMagneticField));

    fluxVector[9] = computedInterfaceAuxiliaryField;
    fluxVector[10] = (computedInterfaceYMagneticField * computedInterfaceXVelocity) - (computedInterfaceXMagneticField * computedInterfaceYVelocity);
    fluxVector[11] = (computedInterfaceZMagneticField * computedInterfaceXVelocity) - (computedInterfaceXMagneticField * computedInterfaceZVelocity);

    fluxVector[12] = hyperbolicWaveSpeedSquared * computedInterfaceXMagneticField;

    return fluxVector;
}

vector<double> MHDIntermediateStateVector::computeXFluxVector(MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    return computeXFluxVector(computeConservedVariableVector(material1Parameters, material2Parameters), material1Parameters, material2Parameters);
}

vector<double> MHDIntermediateStateVector::computeYFluxVector(vector<double> conservedVariableVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> fluxVector(13);

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
    double computedMaterial2Density = conservedVariableVector[7] / computedMaterial2VolumeFraction;

    double computedInterfaceXMagneticField = conservedVariableVector[9];
    double computedInterfaceYMagneticField = conservedVariableVector[10];
    double computedInterfaceZMagneticField = conservedVariableVector[11];

    double computedInterfaceAuxiliaryField = conservedVariableVector[12];

    double velocitySquared = (computedInterfaceXVelocity * computedInterfaceXVelocity) + (computedInterfaceYVelocity * computedInterfaceYVelocity) +
            (computedInterfaceZVelocity * computedInterfaceZVelocity);
    double magneticFieldSquared = (computedInterfaceXMagneticField * computedInterfaceXMagneticField) + (computedInterfaceYMagneticField * computedInterfaceYMagneticField) +
            (computedInterfaceZMagneticField * computedInterfaceZMagneticField);

    double material1TotalEnergy = conservedVariableVector[6] / computedMaterial1VolumeFraction;
    double material1TotalHydrodynamicEnergy = (material1TotalEnergy - (0.5 * magneticFieldSquared)) / computedMaterial1Density;
    double material1SpecificInternalEnergy = material1TotalHydrodynamicEnergy - (0.5 * velocitySquared);

    double material2TotalEnergy = conservedVariableVector[8] / computedMaterial2VolumeFraction;
    double material2TotalHydrodynamicEnergy = (material2TotalEnergy - (0.5 * magneticFieldSquared)) / computedMaterial2Density;
    double material2SpecificInternalEnergy = material2TotalHydrodynamicEnergy - (0.5 * velocitySquared);

    double computedMaterial1Pressure = MHDEquationOfState::computePressure(computedMaterial1Density, material1SpecificInternalEnergy, material1Parameters);
    double computedMaterial2Pressure = MHDEquationOfState::computePressure(computedMaterial2Density, material2SpecificInternalEnergy, material2Parameters);

    double material1HyperbolicWaveSpeedSquared = material1Parameters.computeHyperbolicWaveSpeedSquared();
    double material2HyperbolicWaveSpeedSquared = material2Parameters.computeHyperbolicWaveSpeedSquared();
    double hyperbolicWaveSpeedSquared = (computedMaterial1VolumeFraction * material1HyperbolicWaveSpeedSquared) + (computedMaterial2VolumeFraction * material2HyperbolicWaveSpeedSquared);

    fluxVector[0] = computedTotalDensity * computedInterfaceYVelocity;
    fluxVector[1] = computedTotalDensity * (computedInterfaceYVelocity * computedMaterial1VolumeFraction);

    fluxVector[2] = (computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceXVelocity)) - (computedInterfaceYMagneticField * computedInterfaceXMagneticField);
    fluxVector[3] = (computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceYVelocity)) + ((computedMaterial1VolumeFraction * computedMaterial1Pressure) +
                                                                                                          (computedMaterial2VolumeFraction * computedMaterial2Pressure)) +
            (0.5 * magneticFieldSquared) - (computedInterfaceYMagneticField * computedInterfaceYMagneticField);
    fluxVector[4] = (computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceZVelocity)) - (computedInterfaceYMagneticField * computedInterfaceZMagneticField);

    fluxVector[5] = computedMaterial1VolumeFraction * (computedMaterial1Density * computedInterfaceYVelocity);
    fluxVector[7] = computedMaterial2VolumeFraction * (computedMaterial2Density * computedInterfaceYVelocity);

    vector<double> velocityVector(3);
    velocityVector[0] = computedInterfaceXVelocity;
    velocityVector[1] = computedInterfaceYVelocity;
    velocityVector[2] = computedInterfaceZVelocity;

    vector<double> magneticFieldVector(3);
    magneticFieldVector[0] = computedInterfaceXMagneticField;
    magneticFieldVector[1] = computedInterfaceYMagneticField;
    magneticFieldVector[2] = computedInterfaceZMagneticField;

    double magneticFieldVelocityVectorProduct = VectorAlgebra::computeDotProduct(velocityVector, magneticFieldVector);
    fluxVector[6] = (computedMaterial1VolumeFraction * ((material1TotalEnergy + computedMaterial1Pressure + (0.5 * magneticFieldSquared)) * computedInterfaceYVelocity)) -
            (computedMaterial1VolumeFraction * (magneticFieldVelocityVectorProduct * computedInterfaceYMagneticField));

    fluxVector[8] = (computedMaterial2VolumeFraction * ((material2TotalEnergy + computedMaterial2Pressure + (0.5 * magneticFieldSquared)) * computedInterfaceYVelocity)) -
            (computedMaterial2VolumeFraction * (magneticFieldVelocityVectorProduct * computedInterfaceYMagneticField));

    fluxVector[9] = (computedInterfaceXMagneticField * computedInterfaceYVelocity) - (computedInterfaceYMagneticField * computedInterfaceXVelocity);
    fluxVector[10] = computedInterfaceAuxiliaryField;
    fluxVector[11] = (computedInterfaceZMagneticField * computedInterfaceYVelocity) - (computedInterfaceYMagneticField * computedInterfaceZVelocity);

    fluxVector[12] = hyperbolicWaveSpeedSquared * computedInterfaceYMagneticField;

    return fluxVector;
}

vector<double> MHDIntermediateStateVector::computeYFluxVector(MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    return computeYFluxVector(computeConservedVariableVector(material1Parameters, material2Parameters), material1Parameters, material2Parameters);
}

vector<double> MHDIntermediateStateVector::computeSourceTermVector(vector<double> conservedVariableVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> sourceTermVector(13);

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

    double computedInterfaceAuxiliaryField = conservedVariableVector[12];

    sourceTermVector[12] = - (hyperbolicWaveSpeedSquared / parabolicDampingSquared) * computedInterfaceAuxiliaryField;

    return sourceTermVector;
}

vector<double> MHDIntermediateStateVector::computeSourceTermVector(MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    return computeSourceTermVector(computeConservedVariableVector(material1Parameters, material2Parameters), material1Parameters, material2Parameters);
}

double MHDIntermediateStateVector::computeMaterial1SpecificInternalEnergy(MHDMaterialParameters material1Parameters)
{
    return MHDEquationOfState::computeSpecificInternalEnergy(material1Density, material1Pressure, material1Parameters);
}

double MHDIntermediateStateVector::computeMaterial1TotalEnergy(MHDMaterialParameters material1Parameters)
{
    double velocitySquared = (interfaceXVelocity * interfaceXVelocity) + (interfaceYVelocity * interfaceYVelocity) + (interfaceZVelocity * interfaceZVelocity);
    double totalHydrodynamicEnergy = material1Density * ((0.5 * velocitySquared) + computeMaterial1SpecificInternalEnergy(material1Parameters));

    double magneticFieldSquared = (interfaceXMagneticField * interfaceXMagneticField) + (interfaceYMagneticField * interfaceYMagneticField) + (interfaceZMagneticField * interfaceZMagneticField);

    return totalHydrodynamicEnergy + (0.5 * magneticFieldSquared);
}

double MHDIntermediateStateVector::computeMaterial1SoundSpeed(MHDMaterialParameters material1Parameters)
{
    return MHDEquationOfState::computeSoundSpeed(material1Density, material1Pressure, material1Parameters);
}

double MHDIntermediateStateVector::computeMaterial1Entropy(MHDMaterialParameters material1Parameters)
{
    return MHDEquationOfState::computeEntropy(material1Density, material1Pressure, material1Parameters);
}

double MHDIntermediateStateVector::computeMaterial1AlfvenWaveSpeed()
{
    return MHDWaveSpeeds::computeAlfvenWaveSpeed(material1Density, interfaceXMagneticField, interfaceYMagneticField, interfaceZMagneticField);
}

double MHDIntermediateStateVector::computeMaterial1XSlowMagnetoAcousticSpeed(MHDMaterialParameters material1Parameters)
{
    return MHDWaveSpeeds::computeXSlowMagnetoAcousticSpeed(material1Density, material1Pressure, interfaceXMagneticField, interfaceYMagneticField, interfaceZMagneticField, material1Parameters);
}

double MHDIntermediateStateVector::computeMaterial1YSlowMagnetoAcousticSpeed(MHDMaterialParameters material1Parameters)
{
    return MHDWaveSpeeds::computeYSlowMagnetoAcousticSpeed(material1Density, material1Pressure, interfaceXMagneticField, interfaceYMagneticField, interfaceZMagneticField, material1Parameters);
}

double MHDIntermediateStateVector::computeMaterial1XFastMagnetoAcousticSpeed(MHDMaterialParameters material1Parameters)
{
    return MHDWaveSpeeds::computeXFastMagnetoAcousticSpeed(material1Density, material1Pressure, interfaceXMagneticField, interfaceYMagneticField, interfaceZMagneticField, material1Parameters);
}

double MHDIntermediateStateVector::computeMaterial1YFastMagnetoAcousticSpeed(MHDMaterialParameters material1Parameters)
{
    return MHDWaveSpeeds::computeYFastMagnetoAcousticSpeed(material1Density, material1Pressure, interfaceXMagneticField, interfaceYMagneticField, interfaceZMagneticField, material1Parameters);
}

double MHDIntermediateStateVector::computeMaterial2SpecificInternalEnergy(MHDMaterialParameters material2Parameters)
{
    return MHDEquationOfState::computeSpecificInternalEnergy(material2Density, material2Pressure, material2Parameters);
}

double MHDIntermediateStateVector::computeMaterial2TotalEnergy(MHDMaterialParameters material2Parameters)
{
    double velocitySquared = (interfaceXVelocity * interfaceXVelocity) + (interfaceYVelocity * interfaceYVelocity) + (interfaceZVelocity * interfaceZVelocity);
    double totalHydrodynamicEnergy = material2Density * ((0.5 * velocitySquared) + computeMaterial2SpecificInternalEnergy(material2Parameters));

    double magneticFieldSquared = (interfaceXMagneticField * interfaceXMagneticField) + (interfaceYMagneticField * interfaceYMagneticField) + (interfaceZMagneticField * interfaceZMagneticField);

    return totalHydrodynamicEnergy + (0.5 * magneticFieldSquared);
}

double MHDIntermediateStateVector::computeMaterial2SoundSpeed(MHDMaterialParameters material2Parameters)
{
    return MHDEquationOfState::computeSoundSpeed(material2Density, material2Pressure, material2Parameters);
}

double MHDIntermediateStateVector::computeMaterial2Entropy(MHDMaterialParameters material2Parameters)
{
    return MHDEquationOfState::computeEntropy(material2Density, material2Pressure, material2Parameters);
}

double MHDIntermediateStateVector::computeMaterial2AlfvenWaveSpeed()
{
    return MHDWaveSpeeds::computeAlfvenWaveSpeed(material2Density, interfaceXMagneticField, interfaceYMagneticField, interfaceZMagneticField);
}

double MHDIntermediateStateVector::computeMaterial2XSlowMagnetoAcousticSpeed(MHDMaterialParameters material2Parameters)
{
    return MHDWaveSpeeds::computeXSlowMagnetoAcousticSpeed(material2Density, material2Pressure, interfaceXMagneticField, interfaceYMagneticField, interfaceZMagneticField, material2Parameters);
}

double MHDIntermediateStateVector::computeMaterial2YSlowMagnetoAcousticSpeed(MHDMaterialParameters material2Parameters)
{
    return MHDWaveSpeeds::computeYSlowMagnetoAcousticSpeed(material2Density, material2Pressure, interfaceXMagneticField, interfaceYMagneticField, interfaceZMagneticField, material2Parameters);
}

double MHDIntermediateStateVector::computeMaterial2XFastMagnetoAcousticSpeed(MHDMaterialParameters material2Parameters)
{
    return MHDWaveSpeeds::computeXFastMagnetoAcousticSpeed(material2Density, material2Pressure, interfaceXMagneticField, interfaceYMagneticField, interfaceZMagneticField, material2Parameters);
}

double MHDIntermediateStateVector::computeMaterial2YFastMagnetoAcousticSpeed(MHDMaterialParameters material2Parameters)
{
    return MHDWaveSpeeds::computeYFastMagnetoAcousticSpeed(material2Density, material2Pressure, interfaceXMagneticField, interfaceYMagneticField, interfaceZMagneticField, material2Parameters);
}

double MHDIntermediateStateVector::computeTotalDensity()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1Density) + (material2VolumeFraction * material2Density);
}

double MHDIntermediateStateVector::computeTotalPressure()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1Pressure) + (material2VolumeFraction * material2Pressure);
}

void MHDIntermediateStateVector::relaxTotalDensity()
{
    double totalDensity = computeTotalDensity();

    material1Density = totalDensity;
    material2Density = totalDensity;
}

void MHDIntermediateStateVector::relaxTotalPressure()
{
    double totalPressure = computeTotalPressure();

    material1Pressure = totalPressure;
    material2Pressure = totalPressure;
}

void MHDIntermediateStateVector::setMaterial1VolumeFraction(double newMaterial1VolumeFraction)
{
    material1VolumeFraction = newMaterial1VolumeFraction;
}

void MHDIntermediateStateVector::setInterfaceXVelocity(double newInterfaceXVelocity)
{
    interfaceXVelocity = newInterfaceXVelocity;
}

void MHDIntermediateStateVector::setInterfaceYVelocity(double newInterfaceYVelocity)
{
    interfaceYVelocity = newInterfaceYVelocity;
}

void MHDIntermediateStateVector::setInterfaceZVelocity(double newInterfaceZVelocity)
{
    interfaceZVelocity = newInterfaceZVelocity;
}

void MHDIntermediateStateVector::setMaterial1Density(double newMaterial1Density)
{
    material1Density = newMaterial1Density;
}

void MHDIntermediateStateVector::setMaterial1Pressure(double newMaterial1Pressure)
{
    material1Pressure = newMaterial1Pressure;
}

void MHDIntermediateStateVector::setMaterial2Density(double newMaterial2Density)
{
    material2Density = newMaterial2Density;
}

void MHDIntermediateStateVector::setMaterial2Pressure(double newMaterial2Pressure)
{
    material2Pressure = newMaterial2Pressure;
}

void MHDIntermediateStateVector::setInterfaceXMagneticField(double newInterfaceXMagneticField)
{
    interfaceXMagneticField = newInterfaceXMagneticField;
}

void MHDIntermediateStateVector::setInterfaceYMagneticField(double newInterfaceYMagneticField)
{
    interfaceYMagneticField = newInterfaceYMagneticField;
}

void MHDIntermediateStateVector::setInterfaceZMagneticField(double newInterfaceZMagneticField)
{
    interfaceZMagneticField = newInterfaceZMagneticField;
}

void MHDIntermediateStateVector::setInterfaceAuxiliaryField(double newInterfaceAuxiliaryField)
{
    interfaceAuxiliaryField = newInterfaceAuxiliaryField;
}

double MHDIntermediateStateVector::getMaterial1VolumeFraction()
{
    return material1VolumeFraction;
}

double MHDIntermediateStateVector::getInterfaceXVelocity()
{
    return interfaceXVelocity;
}

double MHDIntermediateStateVector::getInterfaceYVelocity()
{
    return interfaceYVelocity;
}

double MHDIntermediateStateVector::getInterfaceZVelocity()
{
    return interfaceZVelocity;
}

double MHDIntermediateStateVector::getMaterial1Density()
{
    return material1Density;
}

double MHDIntermediateStateVector::getMaterial1Pressure()
{
    return material1Pressure;
}

double MHDIntermediateStateVector::getMaterial2Density()
{
    return material2Density;
}

double MHDIntermediateStateVector::getMaterial2Pressure()
{
    return material2Pressure;
}

double MHDIntermediateStateVector::getInterfaceXMagneticField()
{
    return interfaceXMagneticField;
}

double MHDIntermediateStateVector::getInterfaceYMagneticField()
{
    return interfaceYMagneticField;
}

double MHDIntermediateStateVector::getInterfaceZMagneticField()
{
    return interfaceZMagneticField;
}

double MHDIntermediateStateVector::getInterfaceAuxiliaryField()
{
    return interfaceAuxiliaryField;
}

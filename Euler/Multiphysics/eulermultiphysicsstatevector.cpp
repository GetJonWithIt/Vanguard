#include "eulermultiphysicsstatevector.h"

EulerMultiphysicsStateVector::EulerMultiphysicsStateVector()
{
    material1VolumeFraction = 0.999;
    interfaceXVelocity = 0.0;
    interfaceYVelocity = 0.0;
    interfaceZVelocity = 0.0;

    material1Density = 1.0;
    material1Pressure = 1.0;

    material2Density = 1.0;
    material2Pressure = 1.0;
}

EulerMultiphysicsStateVector::EulerMultiphysicsStateVector(double newMaterial1VolumeFraction, double newInterfaceXVelocity, double newInterfaceYVelocity, double newInterfaceZVelocity,
                                                           double newMaterial1Density, double newMaterial1Pressure, double newMaterial2Density, double newMaterial2Pressure)
{
    material1VolumeFraction = newMaterial1VolumeFraction;
    interfaceXVelocity = newInterfaceXVelocity;
    interfaceYVelocity = newInterfaceYVelocity;
    interfaceZVelocity = newInterfaceZVelocity;

    material1Density = newMaterial1Density;
    material1Pressure = newMaterial1Pressure;

    material2Density = newMaterial2Density;
    material2Pressure = newMaterial2Pressure;
}

void EulerMultiphysicsStateVector::setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector)
{
    material1VolumeFraction = newPrimitiveVariableVector[0];
    interfaceXVelocity = newPrimitiveVariableVector[1];
    interfaceYVelocity = newPrimitiveVariableVector[2];
    interfaceZVelocity = newPrimitiveVariableVector[3];

    material1Density = newPrimitiveVariableVector[4];
    material1Pressure = newPrimitiveVariableVector[5];

    material2Density = newPrimitiveVariableVector[6];
    material2Pressure = newPrimitiveVariableVector[7];
}

void EulerMultiphysicsStateVector::setConservedVariableVector(vector<double> newConservedVariableVector, EulerMaterialParameters material1Parameters,
                                                              EulerMaterialParameters material2Parameters)
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

    double velocitySquared = (interfaceXVelocity * interfaceXVelocity) + (interfaceYVelocity * interfaceYVelocity) + (interfaceZVelocity + interfaceZVelocity);
    double material1SpecificInternalEnergy = (newConservedVariableVector[6] / (material1Density * material1VolumeFraction)) - (0.5 * velocitySquared);
    material1Pressure = EulerEquationOfState::computePressure(material1Density, material1SpecificInternalEnergy, material1Parameters);

    material2Density = newConservedVariableVector[7] / material2VolumeFraction;

    double material2SpecificInternalEnergy = (newConservedVariableVector[8] / (material2Density * material2VolumeFraction)) - (0.5 * velocitySquared);
    material2Pressure = EulerEquationOfState::computePressure(material2Density, material2SpecificInternalEnergy, material2Parameters);
}

vector<double> EulerMultiphysicsStateVector::computePrimitiveVariableVector()
{
    vector<double> primitiveVariableVector(8);

    primitiveVariableVector[0] = material1VolumeFraction;
    primitiveVariableVector[1] = interfaceXVelocity;
    primitiveVariableVector[2] = interfaceYVelocity;
    primitiveVariableVector[3] = interfaceZVelocity;

    primitiveVariableVector[4] = material1Density;
    primitiveVariableVector[5] = material1Pressure;

    primitiveVariableVector[6] = material2Density;
    primitiveVariableVector[7] = material2Pressure;

    return primitiveVariableVector;
}

vector<double> EulerMultiphysicsStateVector::computeConservedVariableVector(EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    vector<double> conservedVariableVector(9);

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
    conservedVariableVector[6] = material1VolumeFraction * computeMaterial1TotalEnergy(material1Parameters);

    conservedVariableVector[7] = material2VolumeFraction * material2Density;
    conservedVariableVector[8] = material2VolumeFraction * computeMaterial2TotalEnergy(material2Parameters);

    return conservedVariableVector;
}

vector<double> EulerMultiphysicsStateVector::computeXFluxVector(vector<double> conservedVariableVector, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    vector<double> fluxVector(9);

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

    double velocitySquared = (computedInterfaceXVelocity * computedInterfaceXVelocity) + (computedInterfaceYVelocity * computedInterfaceYVelocity) +
            (computedInterfaceZVelocity * computedInterfaceZVelocity);
    double material1TotalEnergy = conservedVariableVector[6] / (computedMaterial1VolumeFraction * computedMaterial1Density);
    double material1SpecificInternalEnergy = material1TotalEnergy - (0.5 * velocitySquared);
    double computedMaterial1Pressure = EulerEquationOfState::computePressure(computedMaterial1Density, material1SpecificInternalEnergy, material1Parameters);

    double computedMaterial2Density = conservedVariableVector[7] / computedMaterial2VolumeFraction;

    double material2TotalEnergy = conservedVariableVector[8] / (computedMaterial2VolumeFraction * computedMaterial2Density);
    double material2SpecificInternalEnergy = material2TotalEnergy - (0.5 * velocitySquared);
    double computedMaterial2Pressure = EulerEquationOfState::computePressure(computedMaterial2Density, material2SpecificInternalEnergy, material2Parameters);

    fluxVector[0] = computedTotalDensity * computedInterfaceXVelocity;
    fluxVector[1] = computedTotalDensity * (computedInterfaceXVelocity * computedMaterial1VolumeFraction);
    fluxVector[2] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceXVelocity)) + ((computedMaterial1VolumeFraction * computedMaterial1Pressure) +
                                                                                                          (computedMaterial2VolumeFraction * computedMaterial2Pressure));
    fluxVector[3] = computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceYVelocity);
    fluxVector[4] = computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceZVelocity);

    fluxVector[5] = computedMaterial1VolumeFraction * (computedMaterial1Density * computedInterfaceXVelocity);
    fluxVector[6] = computedMaterial1VolumeFraction * ((computedMaterial1Density * (computedInterfaceXVelocity * material1TotalEnergy)) + (computedInterfaceXVelocity * computedMaterial1Pressure));

    fluxVector[7] = computedMaterial2VolumeFraction * (computedMaterial2Density * computedInterfaceXVelocity);
    fluxVector[8] = computedMaterial2VolumeFraction * ((computedMaterial2Density * (computedInterfaceXVelocity * material2TotalEnergy)) + (computedInterfaceXVelocity * computedMaterial2Pressure));

    return fluxVector;
}

vector<double> EulerMultiphysicsStateVector::computeXFluxVector(EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    return computeXFluxVector(computeConservedVariableVector(material1Parameters, material2Parameters), material1Parameters, material2Parameters);
}

vector<double> EulerMultiphysicsStateVector::computeYFluxVector(vector<double> conservedVariableVector, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    vector<double> fluxVector(9);

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

    double velocitySquared = (computedInterfaceXVelocity * computedInterfaceXVelocity) + (computedInterfaceYVelocity * computedInterfaceYVelocity) +
            (computedInterfaceZVelocity * computedInterfaceZVelocity);
    double material1TotalEnergy = conservedVariableVector[6] / (computedMaterial1VolumeFraction * computedMaterial1Density);
    double material1SpecificInternalEnergy = material1TotalEnergy - (0.5 * velocitySquared);
    double computedMaterial1Pressure = EulerEquationOfState::computePressure(computedMaterial1Density, material1SpecificInternalEnergy, material1Parameters);

    double computedMaterial2Density = conservedVariableVector[7] / computedMaterial2VolumeFraction;

    double material2TotalEnergy = conservedVariableVector[8] / (computedMaterial2VolumeFraction * computedMaterial2Density);
    double material2SpecificInternalEnergy = material2TotalEnergy - (0.5 * velocitySquared);
    double computedMaterial2Pressure = EulerEquationOfState::computePressure(computedMaterial2Density, material2SpecificInternalEnergy, material2Parameters);

    fluxVector[0] = computedTotalDensity * computedInterfaceYVelocity;
    fluxVector[1] = computedTotalDensity * (computedInterfaceYVelocity * computedMaterial1VolumeFraction);
    fluxVector[2] = computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceXVelocity);
    fluxVector[3] = (computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceYVelocity)) + ((computedMaterial1VolumeFraction * computedMaterial1Pressure) +
                                                                                                          (computedMaterial2VolumeFraction * computedMaterial2Pressure));
    fluxVector[4] = computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceZVelocity);

    fluxVector[5] = computedMaterial1VolumeFraction * (computedMaterial1Density * computedInterfaceYVelocity);
    fluxVector[6] = computedMaterial1VolumeFraction * ((computedMaterial1Density * (computedInterfaceYVelocity * material1TotalEnergy)) + (computedInterfaceYVelocity * computedMaterial1Pressure));

    fluxVector[7] = computedMaterial2VolumeFraction * (computedMaterial2Density * computedInterfaceYVelocity);
    fluxVector[8] = computedMaterial2VolumeFraction * ((computedMaterial2Density * (computedInterfaceYVelocity * material2TotalEnergy)) + (computedInterfaceYVelocity * computedMaterial2Pressure));

    return fluxVector;
}

vector<double> EulerMultiphysicsStateVector::computeYFluxVector(EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    return computeYFluxVector(computeConservedVariableVector(material1Parameters, material2Parameters), material1Parameters, material2Parameters);
}

double EulerMultiphysicsStateVector::computeMaterial1SpecificInternalEnergy(EulerMaterialParameters material1Parameters)
{
    return EulerEquationOfState::computeSpecificInternalEnergy(material1Density, material1Pressure, material1Parameters);
}

double EulerMultiphysicsStateVector::computeMaterial1TotalEnergy(EulerMaterialParameters material1Parameters)
{
    double velocitySquared = (interfaceXVelocity * interfaceXVelocity) + (interfaceYVelocity * interfaceYVelocity) + (interfaceZVelocity * interfaceZVelocity);

    return material1Density * ((0.5 * velocitySquared) + computeMaterial1SpecificInternalEnergy(material1Parameters));
}

double EulerMultiphysicsStateVector::computeMaterial1SoundSpeed(EulerMaterialParameters material1Parameters)
{
    return EulerEquationOfState::computeSoundSpeed(material1Density, material1Pressure, material1Parameters);
}

double EulerMultiphysicsStateVector::computeMaterial1Entropy(EulerMaterialParameters material1Parameters)
{
    return EulerEquationOfState::computeEntropy(material1Density, material1Pressure, material1Parameters);
}

double EulerMultiphysicsStateVector::computeMaterial2SpecificInternalEnergy(EulerMaterialParameters material2Parameters)
{
    return EulerEquationOfState::computeSpecificInternalEnergy(material2Density, material2Pressure, material2Parameters);
}

double EulerMultiphysicsStateVector::computeMaterial2TotalEnergy(EulerMaterialParameters material2Parameters)
{
    double velocitySquared = (interfaceXVelocity * interfaceXVelocity) + (interfaceYVelocity * interfaceYVelocity) + (interfaceZVelocity * interfaceZVelocity);

    return material2Density * ((0.5 * velocitySquared) + computeMaterial2SpecificInternalEnergy(material2Parameters));
}

double EulerMultiphysicsStateVector::computeMaterial2SoundSpeed(EulerMaterialParameters material2Parameters)
{
    return EulerEquationOfState::computeSoundSpeed(material2Density, material2Pressure, material2Parameters);
}

double EulerMultiphysicsStateVector::computeMaterial2Entropy(EulerMaterialParameters material2Parameters)
{
    return EulerEquationOfState::computeEntropy(material2Density, material2Pressure, material2Parameters);
}

double EulerMultiphysicsStateVector::computeTotalDensity()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1Density) + (material2VolumeFraction * material2Density);
}

double EulerMultiphysicsStateVector::computeTotalPressure()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1Pressure) + (material2VolumeFraction * material2Pressure);
}

void EulerMultiphysicsStateVector::relaxTotalDensity()
{
    double totalDensity = computeTotalDensity();

    material1Density = totalDensity;
    material2Density = totalDensity;
}

void EulerMultiphysicsStateVector::relaxTotalPressure(EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;
    double totalDensity = computeTotalDensity();
    double totalSpecificInternalEnergy = ((material1VolumeFraction * material1Density * computeMaterial1SpecificInternalEnergy(material1Parameters)) +
            (material2VolumeFraction * material2Density * computeMaterial2SpecificInternalEnergy(material2Parameters))) / totalDensity;

    material1Pressure = EulerEquationOfState::computePressure(material1Density, totalSpecificInternalEnergy, material1Parameters);
    material2Pressure = EulerEquationOfState::computePressure(material2Density, totalSpecificInternalEnergy, material2Parameters);
}

void EulerMultiphysicsStateVector::setMaterial1VolumeFraction(double newMaterial1VolumeFraction)
{
    material1VolumeFraction = newMaterial1VolumeFraction;
}

void EulerMultiphysicsStateVector::setInterfaceXVelocity(double newInterfaceXVelocity)
{
    interfaceXVelocity = newInterfaceXVelocity;
}

void EulerMultiphysicsStateVector::setInterfaceYVelocity(double newInterfaceYVelocity)
{
    interfaceYVelocity = newInterfaceYVelocity;
}

void EulerMultiphysicsStateVector::setInterfaceZVelocity(double newInterfaceZVelocity)
{
    interfaceZVelocity = newInterfaceZVelocity;
}

void EulerMultiphysicsStateVector::setMaterial1Density(double newMaterial1Density)
{
    material1Density = newMaterial1Density;
}

void EulerMultiphysicsStateVector::setMaterial1Pressure(double newMaterial1Pressure)
{
    material1Pressure = newMaterial1Pressure;
}

void EulerMultiphysicsStateVector::setMaterial2Density(double newMaterial2Density)
{
    material2Density = newMaterial2Density;
}

void EulerMultiphysicsStateVector::setMaterial2Pressure(double newMaterial2Pressure)
{
    material2Pressure = newMaterial2Pressure;
}

double EulerMultiphysicsStateVector::getMaterial1VolumeFraction()
{
    return material1VolumeFraction;
}

double EulerMultiphysicsStateVector::getInterfaceXVelocity()
{
    return interfaceXVelocity;
}

double EulerMultiphysicsStateVector::getInterfaceYVelocity()
{
    return interfaceYVelocity;
}

double EulerMultiphysicsStateVector::getInterfaceZVelocity()
{
    return interfaceZVelocity;
}

double EulerMultiphysicsStateVector::getMaterial1Density()
{
    return material1Density;
}

double EulerMultiphysicsStateVector::getMaterial1Pressure()
{
    return material1Pressure;
}

double EulerMultiphysicsStateVector::getMaterial2Density()
{
    return material2Density;
}

double EulerMultiphysicsStateVector::getMaterial2Pressure()
{
    return material2Pressure;
}

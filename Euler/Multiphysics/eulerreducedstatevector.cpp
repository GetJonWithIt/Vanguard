#include "eulerreducedstatevector.h"

EulerReducedStateVector::EulerReducedStateVector()
{
    material1VolumeFraction = 0.999;
    interfaceXVelocity = 0.0;
    interfaceYVelocity = 0.0;
    interfaceZVelocity = 0.0;

    material1Density = 1.0;
    material2Density = 1.0;

    interfacePressure = 1.0;
}

EulerReducedStateVector::EulerReducedStateVector(double newMaterial1VolumeFraction, double newInterfaceXVelocity, double newInterfaceYVelocity, double newInterfaceZVelocity,
                                                 double newMaterial1Density, double newMaterial2Density, double newInterfacePressure)
{
    material1VolumeFraction = newMaterial1VolumeFraction;
    interfaceXVelocity = newInterfaceXVelocity;
    interfaceYVelocity = newInterfaceYVelocity;
    interfaceZVelocity = newInterfaceZVelocity;

    material1Density = newMaterial1Density;
    material2Density = newMaterial2Density;

    interfacePressure = newInterfacePressure;
}

void EulerReducedStateVector::setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector)
{
    material1VolumeFraction = newPrimitiveVariableVector[0];
    interfaceXVelocity = newPrimitiveVariableVector[1];
    interfaceYVelocity = newPrimitiveVariableVector[2];
    interfaceZVelocity = newPrimitiveVariableVector[3];

    material1Density = newPrimitiveVariableVector[4];
    material2Density = newPrimitiveVariableVector[5];

    interfacePressure = newPrimitiveVariableVector[6];
}

void EulerReducedStateVector::setConservedVariableVector(vector<double> newConservedVariableVector, EulerMaterialParameters material1Parameters,
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
    material2Density = newConservedVariableVector[6] / material2VolumeFraction;

    double velocitySquared = (interfaceXVelocity * interfaceXVelocity) + (interfaceYVelocity * interfaceYVelocity) + (interfaceZVelocity * interfaceZVelocity);
    double totalSpecificInternalEnergy = (newConservedVariableVector[7] / totalDensity) - (0.5 * velocitySquared);

    double material1Pressure = EulerEquationOfState::computePressure(totalDensity, totalSpecificInternalEnergy, material1Parameters);
    double material2Pressure = EulerEquationOfState::computePressure(totalDensity, totalSpecificInternalEnergy, material2Parameters);
    interfacePressure = (material1VolumeFraction * material1Pressure) + (material2VolumeFraction * material2Pressure);
}

vector<double> EulerReducedStateVector::computePrimitiveVariableVector()
{
    vector<double> primitiveVariableVector(7);

    primitiveVariableVector[0] = material1VolumeFraction;
    primitiveVariableVector[1] = interfaceXVelocity;
    primitiveVariableVector[2] = interfaceYVelocity;
    primitiveVariableVector[3] = interfaceZVelocity;

    primitiveVariableVector[4] = material1Density;
    primitiveVariableVector[5] = material2Density;

    primitiveVariableVector[6] = interfacePressure;

    return primitiveVariableVector;
}

vector<double> EulerReducedStateVector::computeConservedVariableVector(EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    vector<double> conservedVariableVector(8);

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

    double material1TotalEnergy = EulerEquationOfState::computeSpecificInternalEnergy(totalDensity, interfacePressure, material1Parameters) + (0.5 * velocitySquared);
    double material2TotalEnergy = EulerEquationOfState::computeSpecificInternalEnergy(totalDensity, interfacePressure, material2Parameters) + (0.5 * velocitySquared);
    double computedTotalEnergy = (material1VolumeFraction * material1TotalEnergy) + (material2VolumeFraction * material2TotalEnergy);

    conservedVariableVector[7] = totalDensity * computedTotalEnergy;
    return conservedVariableVector;
}

vector<double> EulerReducedStateVector::computeXFluxVector(vector<double> conservedVariableVector, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    vector<double> fluxVector(8);

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

    double computedTotalEnergy = conservedVariableVector[7] / computedTotalDensity;
    double velocitySquared = (computedInterfaceXVelocity * computedInterfaceXVelocity) + (computedInterfaceYVelocity * computedInterfaceYVelocity) +
            (computedInterfaceZVelocity * computedInterfaceZVelocity);
    double totalSpecificInternalEnergy = computedTotalEnergy - (0.5 * velocitySquared);

    double computedMaterial1Pressure = EulerEquationOfState::computePressure(computedTotalDensity, totalSpecificInternalEnergy, material1Parameters);
    double computedMaterial2Pressure = EulerEquationOfState::computePressure(computedTotalDensity, totalSpecificInternalEnergy, material2Parameters);
    double computedInterfacePressure = (computedMaterial1VolumeFraction * computedMaterial1Pressure) + (computedMaterial2VolumeFraction * computedMaterial2Pressure);

    fluxVector[0] = computedTotalDensity * computedInterfaceXVelocity;
    fluxVector[1] = computedTotalDensity * (computedInterfaceXVelocity * computedMaterial1VolumeFraction);
    fluxVector[2] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceXVelocity)) + computedInterfacePressure;
    fluxVector[3] = computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceYVelocity);
    fluxVector[4] = computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceZVelocity);

    fluxVector[5] = computedMaterial1VolumeFraction * (computedMaterial1Density * computedInterfaceXVelocity);
    fluxVector[6] = computedMaterial2VolumeFraction * (computedMaterial2Density * computedInterfaceXVelocity);

    fluxVector[7] = (computedTotalDensity * (computedInterfaceXVelocity * computedTotalEnergy)) + (computedInterfaceXVelocity * computedInterfacePressure);

    return fluxVector;
}

vector<double> EulerReducedStateVector::computeXFluxVector(EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    return computeXFluxVector(computeConservedVariableVector(material1Parameters, material2Parameters), material1Parameters, material2Parameters);
}

vector<double> EulerReducedStateVector::computeYFluxVector(vector<double> conservedVariableVector, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    vector<double> fluxVector(8);

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

    double computedTotalEnergy = conservedVariableVector[7] / computedTotalDensity;
    double velocitySquared = (computedInterfaceXVelocity * computedInterfaceXVelocity) + (computedInterfaceYVelocity * computedInterfaceYVelocity) +
            (computedInterfaceZVelocity * computedInterfaceZVelocity);
    double totalSpecificInternalEnergy = computedTotalEnergy - (0.5 * velocitySquared);

    double computedMaterial1Pressure = EulerEquationOfState::computePressure(computedTotalDensity, totalSpecificInternalEnergy, material1Parameters);
    double computedMaterial2Pressure = EulerEquationOfState::computePressure(computedTotalDensity, totalSpecificInternalEnergy, material2Parameters);
    double computedInterfacePressure = (computedMaterial1VolumeFraction * computedMaterial1Pressure) + (computedMaterial2VolumeFraction * computedMaterial2Pressure);

    fluxVector[0] = computedTotalDensity * computedInterfaceYVelocity;
    fluxVector[1] = computedTotalDensity * (computedInterfaceYVelocity * computedMaterial1VolumeFraction);
    fluxVector[2] = computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceXVelocity);
    fluxVector[3] = (computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceYVelocity)) + computedInterfacePressure;
    fluxVector[4] = computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceZVelocity);

    fluxVector[5] = computedMaterial1VolumeFraction * (computedMaterial1Density * computedInterfaceYVelocity);
    fluxVector[6] = computedMaterial2VolumeFraction * (computedMaterial2Density * computedInterfaceYVelocity);

    fluxVector[7] = (computedTotalDensity * (computedInterfaceYVelocity * computedTotalEnergy)) + (computedInterfaceYVelocity * computedInterfacePressure);

    return fluxVector;
}

vector<double> EulerReducedStateVector::computeYFluxVector(EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    return computeYFluxVector(computeConservedVariableVector(material1Parameters, material2Parameters), material1Parameters, material2Parameters);
}

double EulerReducedStateVector::computeMaterial1SpecificInternalEnergy(EulerMaterialParameters material1Parameters)
{
    return EulerEquationOfState::computeSpecificInternalEnergy(material1Density, interfacePressure, material1Parameters);
}

double EulerReducedStateVector::computeMaterial1TotalEnergy(EulerMaterialParameters material1Parameters)
{
    double velocitySquared = (interfaceXVelocity * interfaceXVelocity) + (interfaceYVelocity * interfaceYVelocity) + (interfaceZVelocity * interfaceZVelocity);

    return (0.5 * velocitySquared) + computeMaterial1SpecificInternalEnergy(material1Parameters);
}

double EulerReducedStateVector::computeMaterial1SoundSpeed(EulerMaterialParameters material1Parameters)
{
    return EulerEquationOfState::computeSoundSpeed(material1Density, interfacePressure, material1Parameters);
}

double EulerReducedStateVector::computeMaterial1Entropy(EulerMaterialParameters material1Parameters)
{
    return EulerEquationOfState::computeEntropy(material1Density, interfacePressure, material1Parameters);
}

double EulerReducedStateVector::computeMaterial2SpecificInternalEnergy(EulerMaterialParameters material2Parameters)
{
    return EulerEquationOfState::computeSpecificInternalEnergy(material2Density, interfacePressure, material2Parameters);
}

double EulerReducedStateVector::computeMaterial2TotalEnergy(EulerMaterialParameters material2Parameters)
{
    double velocitySquared = (interfaceXVelocity * interfaceXVelocity) + (interfaceYVelocity * interfaceYVelocity) + (interfaceZVelocity * interfaceZVelocity);

    return (0.5 * velocitySquared) + computeMaterial2SpecificInternalEnergy(material2Parameters);
}

double EulerReducedStateVector::computeMaterial2SoundSpeed(EulerMaterialParameters material2Parameters)
{
    return EulerEquationOfState::computeSoundSpeed(material2Density, interfacePressure, material2Parameters);
}

double EulerReducedStateVector::computeMaterial2Entropy(EulerMaterialParameters material2Parameters)
{
    return EulerEquationOfState::computeEntropy(material2Density, interfacePressure, material2Parameters);
}

double EulerReducedStateVector::computeTotalDensity()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1Density) + (material2VolumeFraction * material2Density);
}

void EulerReducedStateVector::relaxTotalDensity()
{
    double totalDensity = computeTotalDensity();

    material1Density = totalDensity;
    material2Density = totalDensity;
}

void EulerReducedStateVector::setMaterial1VolumeFraction(double newMaterial1VolumeFraction)
{
    material1VolumeFraction = newMaterial1VolumeFraction;
}

void EulerReducedStateVector::setInterfaceXVelocity(double newInterfaceXVelocity)
{
    interfaceXVelocity = newInterfaceXVelocity;
}

void EulerReducedStateVector::setInterfaceYVelocity(double newInterfaceYVelocity)
{
    interfaceYVelocity = newInterfaceYVelocity;
}

void EulerReducedStateVector::setInterfaceZVelocity(double newInterfaceZVelocity)
{
    interfaceZVelocity = newInterfaceZVelocity;
}

void EulerReducedStateVector::setMaterial1Density(double newMaterial1Density)
{
    material1Density = newMaterial1Density;
}

void EulerReducedStateVector::setMaterial2Density(double newMaterial2Density)
{
    material2Density = newMaterial2Density;
}

void EulerReducedStateVector::setInterfacePressure(double newInterfacePressure)
{
    interfacePressure = newInterfacePressure;
}

double EulerReducedStateVector::getMaterial1VolumeFraction()
{
    return material1VolumeFraction;
}

double EulerReducedStateVector::getInterfaceXVelocity()
{
    return interfaceXVelocity;
}

double EulerReducedStateVector::getInterfaceYVelocity()
{
    return interfaceYVelocity;
}

double EulerReducedStateVector::getInterfaceZVelocity()
{
    return interfaceZVelocity;
}

double EulerReducedStateVector::getMaterial1Density()
{
    return material1Density;
}

double EulerReducedStateVector::getMaterial2Density()
{
    return material2Density;
}

double EulerReducedStateVector::getInterfacePressure()
{
    return interfacePressure;
}

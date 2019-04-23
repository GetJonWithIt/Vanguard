#include "hprmaterialparameters.h"

HPRMaterialParameters::HPRMaterialParameters()
{
    equationOfState = "StiffenedGas";
    isThermal = false;
    isPlastic = false;

    idealGasConstant = 8.31445985;
    referenceDensity = 1.0;
    referencePressure = 1.0 / 1.4;
    referenceTemperature = 0.0;
    initialTemperature = 0.0;
    specificHeatCapacity = 1.0;
    adiabaticIndex = 1.4;
    stiffeningParameter = 0.0;
    referenceSoundSpeed = 0.0;
    referenceGruneisenCoefficient = 0.0;
    hugoniotSlopeCoefficient = 0.0;
    referenceInternalEnergy = 0.0;
    transverseWaveSpeed = 1.0;

    dynamicViscosityCoefficient = pow(10.0, -2.0);
    prandtlNumber = 0.75;
    elasticPlasticTransitionParameter = 0.0;
    powerLawIndex = 1.0;
    heatWaveSpeed = pow(10.0, -16.0);

    alphaParameter = 0.0;
    betaParameter = 0.0;
    gammaParameter = 0.0;
}

HPRMaterialParameters::HPRMaterialParameters(string newEquationOfState, bool newIsThermal, bool newIsPlastic, double newIdealGasConstant, double newReferenceDensity, double newReferencePressure,
                                             double newReferenceTemperature, double newInitialTemperature, double newSpecificHeatCapacity, double newAdiabaticIndex, double newStiffeningParameter,
                                             double newReferenceSoundSpeed, double newReferenceGruneisenCoefficient, double newHugoniotSlopeCoefficient, double newReferenceInternalEnergy,
                                             double newTransverseWaveSpeed, double newStrainDissipationTime, double newThermalImpulseRelaxationTime, double newDynamicViscosityCoefficient,
                                             double newPrandtlNumber, double newElasticPlasticTransitionParameter, double newPowerLawIndex, double newHeatWaveSpeed, double newAlphaParameter,
                                             double newBetaParameter, double newGammaParameter)
{
    equationOfState = newEquationOfState;
    isThermal = newIsThermal;
    isPlastic = newIsPlastic;

    idealGasConstant = newIdealGasConstant;
    referenceDensity = newReferenceDensity;
    referencePressure = newReferencePressure;
    referenceTemperature = newReferenceTemperature;
    initialTemperature = newInitialTemperature;
    specificHeatCapacity = newSpecificHeatCapacity;
    adiabaticIndex = newAdiabaticIndex;
    stiffeningParameter = newStiffeningParameter;
    referenceSoundSpeed = newReferenceSoundSpeed;
    referenceGruneisenCoefficient = newReferenceGruneisenCoefficient;
    hugoniotSlopeCoefficient = newHugoniotSlopeCoefficient;
    transverseWaveSpeed = newTransverseWaveSpeed;

    strainDissipationTime = newStrainDissipationTime;
    thermalImpulseRelaxationTime = newThermalImpulseRelaxationTime;
    dynamicViscosityCoefficient = newDynamicViscosityCoefficient;
    prandtlNumber = newPrandtlNumber;
    elasticPlasticTransitionParameter = newElasticPlasticTransitionParameter;
    powerLawIndex = newPowerLawIndex;
    heatWaveSpeed = newHeatWaveSpeed;

    alphaParameter = newAlphaParameter;
    betaParameter = newBetaParameter;
    gammaParameter = newGammaParameter;
}

void HPRMaterialParameters::configureStrainDissipationTime()
{
    strainDissipationTime = computeStrainDissipationTime();
}

void HPRMaterialParameters::configureThermalImpulseRelaxationTime()
{
    thermalImpulseRelaxationTime = computeThermalImpulseRelaxationTime();
}

void HPRMaterialParameters::configureThermalImpulseRelaxationTime(double heatConductionCoefficient)
{
    thermalImpulseRelaxationTime = computeThermalImpulseRelaxationTime(heatConductionCoefficient);
}

double HPRMaterialParameters::computeReferenceSoundSpeedSquared()
{
    return referenceSoundSpeed * referenceSoundSpeed;
}

double HPRMaterialParameters::computeTransverseWaveSpeedSquared()
{
    return transverseWaveSpeed * transverseWaveSpeed;
}

double HPRMaterialParameters::computeHeatWaveSpeedSquared()
{
    return heatWaveSpeed * heatWaveSpeed;
}

double HPRMaterialParameters::computeStrainDissipationTime()
{
    return 6.0 * (dynamicViscosityCoefficient / (referencePressure * (transverseWaveSpeed * transverseWaveSpeed)));
}

double HPRMaterialParameters::computeHeatConductionCoefficient()
{
    return dynamicViscosityCoefficient * adiabaticIndex * (specificHeatCapacity / prandtlNumber);
}

double HPRMaterialParameters::computeThermalImpulseRelaxationTime()
{
    return computeThermalImpulseRelaxationTime(computeHeatConductionCoefficient());
}

double HPRMaterialParameters::computeThermalImpulseRelaxationTime(double heatConductionCoefficient)
{
    return heatConductionCoefficient * (referenceDensity / (initialTemperature * (heatWaveSpeed * heatWaveSpeed)));
}

void HPRMaterialParameters::setEquationOfState(string newEquationOfState)
{
    equationOfState = newEquationOfState;
}

void HPRMaterialParameters::setIsThermal(bool newIsThermal)
{
    isThermal = newIsThermal;
}

void HPRMaterialParameters::setIsPlastic(bool newIsPlastic)
{
    isPlastic = newIsPlastic;
}

void HPRMaterialParameters::setIdealGasConstant(double newIdealGasConstant)
{
    idealGasConstant = newIdealGasConstant;
}

void HPRMaterialParameters::setReferenceDensity(double newReferenceDensity)
{
    referenceDensity = newReferenceDensity;
}

void HPRMaterialParameters::setReferencePressure(double newReferencePressure)
{
    referencePressure = newReferencePressure;
}

void HPRMaterialParameters::setReferenceTemperature(double newReferenceTemperature)
{
    referenceTemperature = newReferenceTemperature;
}

void HPRMaterialParameters::setInitialTemperature(double newInitialTemperature)
{
    initialTemperature = newInitialTemperature;
}

void HPRMaterialParameters::setSpecificHeatCapacity(double newSpecificHeatCapacity)
{
    specificHeatCapacity = newSpecificHeatCapacity;
}

void HPRMaterialParameters::setAdiabaticIndex(double newAdiabaticIndex)
{
    adiabaticIndex = newAdiabaticIndex;
}

void HPRMaterialParameters::setStiffeningParameter(double newStiffeningParameter)
{
    stiffeningParameter = newStiffeningParameter;
}

void HPRMaterialParameters::setReferenceSoundSpeed(double newReferenceSoundSpeed)
{
    referenceSoundSpeed = newReferenceSoundSpeed;
}

void HPRMaterialParameters::setReferenceGruneisenCoefficient(double newReferenceGruneisenCoefficient)
{
    referenceGruneisenCoefficient = newReferenceGruneisenCoefficient;
}

void HPRMaterialParameters::setHugoniotSlopeCoefficient(double newHugoniotSlopeCoefficient)
{
    hugoniotSlopeCoefficient = newHugoniotSlopeCoefficient;
}

void HPRMaterialParameters::setTransverseWaveSpeed(double newTransverseWaveSpeed)
{
    transverseWaveSpeed = newTransverseWaveSpeed;
}

void HPRMaterialParameters::setStrainDissipationTime(double newStrainDissipationTime)
{
    strainDissipationTime = newStrainDissipationTime;
}

void HPRMaterialParameters::setThermalImpulseRelaxationTime(double newThermalImpulseRelaxationTime)
{
    thermalImpulseRelaxationTime = newThermalImpulseRelaxationTime;
}

void HPRMaterialParameters::setDynamicViscosityCoefficient(double newDynamicViscosityCoefficient)
{
    dynamicViscosityCoefficient = newDynamicViscosityCoefficient;
}

void HPRMaterialParameters::setPrandtlNumber(double newPrandtlNumber)
{
    prandtlNumber = newPrandtlNumber;
}

void HPRMaterialParameters::setElasticPlasticTransitionParameter(double newElasticPlasticTransitionParameter)
{
    elasticPlasticTransitionParameter = newElasticPlasticTransitionParameter;
}

void HPRMaterialParameters::setPowerLawIndex(double newPowerLawIndex)
{
    powerLawIndex = newPowerLawIndex;
}

void HPRMaterialParameters::setHeatWaveSpeed(double newHeatWaveSpeed)
{
    heatWaveSpeed = newHeatWaveSpeed;
}

void HPRMaterialParameters::setAlphaParameter(double newAlphaParameter)
{
    alphaParameter = newAlphaParameter;
}

void HPRMaterialParameters::setBetaParameter(double newBetaParameter)
{
    betaParameter = newBetaParameter;
}

void HPRMaterialParameters::setGammaParameter(double newGammaParameter)
{
    gammaParameter = newGammaParameter;
}

string HPRMaterialParameters::getEquationOfState()
{
    return equationOfState;
}

bool HPRMaterialParameters::getIsThermal()
{
    return isThermal;
}

bool HPRMaterialParameters::getIsPlastic()
{
    return isPlastic;
}

double HPRMaterialParameters::getIdealGasConstant()
{
    return idealGasConstant;
}

double HPRMaterialParameters::getReferenceDensity()
{
    return referenceDensity;
}

double HPRMaterialParameters::getReferencePressure()
{
    return referencePressure;
}

double HPRMaterialParameters::getReferenceTemperature()
{
    return referenceTemperature;
}

double HPRMaterialParameters::getInitialTemperature()
{
    return initialTemperature;
}

double HPRMaterialParameters::getSpecificHeatCapacity()
{
    return specificHeatCapacity;
}

double HPRMaterialParameters::getAdiabaticIndex()
{
    return adiabaticIndex;
}

double HPRMaterialParameters::getStiffeningParameter()
{
    return stiffeningParameter;
}

double HPRMaterialParameters::getReferenceSoundSpeed()
{
    return referenceSoundSpeed;
}

double HPRMaterialParameters::getReferenceGruneisenCoefficient()
{
    return referenceGruneisenCoefficient;
}

double HPRMaterialParameters::getHugoniotSlopeCoefficient()
{
    return hugoniotSlopeCoefficient;
}

double HPRMaterialParameters::getTransverseWaveSpeed()
{
    return transverseWaveSpeed;
}

double HPRMaterialParameters::getStrainDissipationTime()
{
    return strainDissipationTime;
}

double HPRMaterialParameters::getThermalImpulseRelaxationTime()
{
    return thermalImpulseRelaxationTime;
}

double HPRMaterialParameters::getDynamicViscosityCoefficient()
{
    return dynamicViscosityCoefficient;
}

double HPRMaterialParameters::getPrandtlNumber()
{
    return prandtlNumber;
}

double HPRMaterialParameters::getElasticPlasticTransitionParameter()
{
    return elasticPlasticTransitionParameter;
}

double HPRMaterialParameters::getPowerLawIndex()
{
    return powerLawIndex;
}

double HPRMaterialParameters::getHeatWaveSpeed()
{
    return heatWaveSpeed;
}

double HPRMaterialParameters::getAlphaParameter()
{
    return alphaParameter;
}

double HPRMaterialParameters::getBetaParameter()
{
    return betaParameter;
}

double HPRMaterialParameters::getGammaParameter()
{
    return gammaParameter;
}

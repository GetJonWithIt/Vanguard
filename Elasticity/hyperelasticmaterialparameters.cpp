#include "hyperelasticmaterialparameters.h"

HyperelasticMaterialParameters::HyperelasticMaterialParameters()
{
    referenceMassDensity = 8.9;
    soundSpeed = 4.6;
    shearWaveSpeed = 2.1;
    specificHeatCapacity = 0.0004;
    initialTemperature = 300.0;

    alphaParameter = 1.0;
    betaParameter = 3.0;
    gammaParameter = 2.0;
}

HyperelasticMaterialParameters::HyperelasticMaterialParameters(double newReferenceMassDensity, double newSoundSpeed, double newShearWaveSpeed, double newSpecificHeatCapacity,
                                                               double newInitialTemperature, double newAlphaParameter, double newBetaParameter, double newGammaParameter)
{
    referenceMassDensity = newReferenceMassDensity;
    soundSpeed = newSoundSpeed;
    shearWaveSpeed = newShearWaveSpeed;
    specificHeatCapacity = newSpecificHeatCapacity;
    initialTemperature = newInitialTemperature;

    alphaParameter = newAlphaParameter;
    betaParameter = newBetaParameter;
    gammaParameter = newGammaParameter;
}

double HyperelasticMaterialParameters::computeBulkSoundSpeedSquared()
{
    return (soundSpeed * soundSpeed) - ((4.0 / 3.0) * (shearWaveSpeed * shearWaveSpeed));
}

double HyperelasticMaterialParameters::computeShearWaveSpeedSquared()
{
    return shearWaveSpeed * shearWaveSpeed;
}

void HyperelasticMaterialParameters::setReferenceMassDensity(double newReferenceMassDensity)
{
    referenceMassDensity = newReferenceMassDensity;
}

void HyperelasticMaterialParameters::setSoundSpeed(double newSoundSpeed)
{
    soundSpeed = newSoundSpeed;
}

void HyperelasticMaterialParameters::setShearWaveSpeed(double newShearWaveSpeed)
{
    shearWaveSpeed = newShearWaveSpeed;
}

void HyperelasticMaterialParameters::setSpecificHeatCapacity(double newSpecificHeatCapacity)
{
    specificHeatCapacity = newSpecificHeatCapacity;
}

void HyperelasticMaterialParameters::setInitialTemperature(double newInitialTemperature)
{
    initialTemperature = newInitialTemperature;
}

void HyperelasticMaterialParameters::setAlphaParameter(double newAlphaParameter)
{
    alphaParameter = newAlphaParameter;
}

void HyperelasticMaterialParameters::setBetaParameter(double newBetaParameter)
{
    betaParameter = newBetaParameter;
}

void HyperelasticMaterialParameters::setGammaParameter(double newGammaParameter)
{
    gammaParameter = newGammaParameter;
}

double HyperelasticMaterialParameters::getReferenceMassDensity()
{
    return referenceMassDensity;
}

double HyperelasticMaterialParameters::getSoundSpeed()
{
    return soundSpeed;
}

double HyperelasticMaterialParameters::getShearWaveSpeed()
{
    return shearWaveSpeed;
}

double HyperelasticMaterialParameters::getSpecificHeatCapacity()
{
    return specificHeatCapacity;
}

double HyperelasticMaterialParameters::getInitialTemperature()
{
    return initialTemperature;
}

double HyperelasticMaterialParameters::getAlphaParameter()
{
    return alphaParameter;
}

double HyperelasticMaterialParameters::getBetaParameter()
{
    return betaParameter;
}

double HyperelasticMaterialParameters::getGammaParameter()
{
    return gammaParameter;
}

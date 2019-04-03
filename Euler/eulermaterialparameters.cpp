#include "eulermaterialparameters.h"

EulerMaterialParameters::EulerMaterialParameters()
{
    adiabaticIndex = 1.0;
    stiffeningParameter = 0.0;

    specificHeatCapacity = 25;
    energyOfFormation = 0.0;
    ignitionTemperature = 0.25;
    reactionRate = 250.0;
}

EulerMaterialParameters::EulerMaterialParameters(double newAdiabaticIndex)
{
    adiabaticIndex = newAdiabaticIndex;
    stiffeningParameter = 0.0;

    specificHeatCapacity = 25;
    energyOfFormation = 0.0;
    ignitionTemperature = 0.25;
    reactionRate = 250.0;
}

EulerMaterialParameters::EulerMaterialParameters(double newAdiabaticIndex, double newStiffeningParameter)
{
    adiabaticIndex = newAdiabaticIndex;
    stiffeningParameter = newStiffeningParameter;

    specificHeatCapacity = 25;
    energyOfFormation = 0.0;
    ignitionTemperature = 0.25;
    reactionRate = 250.0;
}

EulerMaterialParameters::EulerMaterialParameters(double newAdiabaticIndex, double newStiffeningParameter, double newSpecificHeatCapacity, double newEnergyOfFormation, double newIgnitionTemperature,
                                                 double newReactionRate)
{
    adiabaticIndex = newAdiabaticIndex;
    stiffeningParameter = newStiffeningParameter;

    specificHeatCapacity = newSpecificHeatCapacity;
    energyOfFormation = newEnergyOfFormation;
    ignitionTemperature = newIgnitionTemperature;
    reactionRate = newReactionRate;
}

void EulerMaterialParameters::setAdiabaticIndex(double newAdiabaticIndex)
{
    adiabaticIndex = newAdiabaticIndex;
}

void EulerMaterialParameters::setStiffeningParameter(double newStiffeningParameter)
{
    stiffeningParameter = newStiffeningParameter;
}

void EulerMaterialParameters::setSpecificHeatCapacity(double newSpecificHeatCapacity)
{
    specificHeatCapacity = newSpecificHeatCapacity;
}

void EulerMaterialParameters::setEnergyOfFormation(double newEnergyOfFormation)
{
    energyOfFormation = newEnergyOfFormation;
}

void EulerMaterialParameters::setIgnitionTemperature(double newIgnitionTemperature)
{
    ignitionTemperature = newIgnitionTemperature;
}

void EulerMaterialParameters::setReactionRate(double newReactionRate)
{
    reactionRate = newReactionRate;
}

double EulerMaterialParameters::getAdiabaticIndex()
{
    return adiabaticIndex;
}

double EulerMaterialParameters::getStiffeningParameter()
{
    return stiffeningParameter;
}

double EulerMaterialParameters::getSpecificHeatCapacity()
{
    return specificHeatCapacity;
}

double EulerMaterialParameters::getEnergyOfFormation()
{
    return energyOfFormation;
}

double EulerMaterialParameters::getIgnitionTemperature()
{
    return ignitionTemperature;
}

double EulerMaterialParameters::getReactionRate()
{
    return reactionRate;
}

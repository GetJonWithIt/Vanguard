#include "eulerequationofstate.h"

EulerEquationOfState::EulerEquationOfState()
{
}

double EulerEquationOfState::computeSpecificInternalEnergy(double density, double pressure, EulerMaterialParameters materialParameters)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    double stiffeningParameter = materialParameters.getStiffeningParameter();

    return (pressure + (adiabaticIndex * stiffeningParameter)) / ((adiabaticIndex - 1) * density);
}

double EulerEquationOfState::computePressure(double density, double specificInternalEnergy, EulerMaterialParameters materialParameters)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    double stiffeningParameter = materialParameters.getStiffeningParameter();

    return (specificInternalEnergy * (adiabaticIndex - 1) * density) - (adiabaticIndex * stiffeningParameter);
}

double EulerEquationOfState::computeSoundSpeed(double density, double pressure, EulerMaterialParameters materialParameters)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    double stiffeningParameter = materialParameters.getStiffeningParameter();

    return sqrt((adiabaticIndex * (pressure + stiffeningParameter)) / density);
}

double EulerEquationOfState::computeEntropy(double density, double pressure, EulerMaterialParameters materialParameters)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();

    return pressure / pow(density, adiabaticIndex);
}

double EulerEquationOfState::computeTemperature(double specificInternalEnergy, EulerMaterialParameters materialParameters)
{
    double specificHeatCapacity = materialParameters.getSpecificHeatCapacity();

    return specificInternalEnergy / specificHeatCapacity;
}

double EulerEquationOfState::computeReactionRate(double specificInternalEnergy, EulerMaterialParameters materialParameters)
{
    double temperature = computeTemperature(specificInternalEnergy, materialParameters);
    double ignitionTemperature = materialParameters.getIgnitionTemperature();
    double reactionRate = materialParameters.getReactionRate();

    if (temperature > ignitionTemperature)
    {
        return reactionRate;
    }
    else
    {
        return 0.0;
    }
}

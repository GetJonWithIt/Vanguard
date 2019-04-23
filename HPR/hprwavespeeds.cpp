#include "hprwavespeeds.h"

HPRWaveSpeeds::HPRWaveSpeeds()
{
}

double HPRWaveSpeeds::computeAdiabaticSoundSpeed(double density, double pressure, HPRMaterialParameters materialParameters)
{
    double totalEnergyDerivativeDensity = HPRDerivatives::computeTotalEnergyDerivativeDensity(density, pressure, materialParameters);
    double totalEnergyDerivativePressure = HPRDerivatives::computeTotalEnergyDerivativePressure(density, materialParameters);

    return sqrt(((pressure / (density * density)) - totalEnergyDerivativeDensity) / totalEnergyDerivativePressure);
}

double HPRWaveSpeeds::computeHeatCharacteristicSpeed(double density, double temperature, HPRMaterialParameters materialParameters)
{
    double heatWaveSpeedSquared = materialParameters.computeHeatWaveSpeedSquared();
    double specificHeatCapacity = materialParameters.getSpecificHeatCapacity();

    return (sqrt(heatWaveSpeedSquared * (temperature / specificHeatCapacity))) / density;
}

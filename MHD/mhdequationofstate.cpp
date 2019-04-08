#include "mhdequationofstate.h"

MHDEquationOfState::MHDEquationOfState()
{
}

double MHDEquationOfState::computeSpecificInternalEnergy(double density, double pressure, MHDMaterialParameters materialParameters)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    double stiffeningParameter = materialParameters.getStiffeningParameter();

    return (pressure + (adiabaticIndex * stiffeningParameter)) / ((adiabaticIndex - 1) * density);
}

double MHDEquationOfState::computePressure(double density, double specificInternalEnergy, MHDMaterialParameters materialParameters)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    double stiffeningParameter = materialParameters.getStiffeningParameter();

    return (specificInternalEnergy * (adiabaticIndex - 1) * density) - (adiabaticIndex * stiffeningParameter);
}

double MHDEquationOfState::computeSoundSpeed(double density, double pressure, MHDMaterialParameters materialParameters)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();
    double stiffeningParameter = materialParameters.getStiffeningParameter();

    return sqrt((adiabaticIndex * (pressure + stiffeningParameter)) / density);
}

double MHDEquationOfState::computeEntropy(double density, double pressure, MHDMaterialParameters materialParameters)
{
    double adiabaticIndex = materialParameters.getAdiabaticIndex();

    return pressure / pow(density, adiabaticIndex);
}

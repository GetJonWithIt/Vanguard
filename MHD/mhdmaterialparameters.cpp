#include "mhdmaterialparameters.h"

MHDMaterialParameters::MHDMaterialParameters()
{
    adiabaticIndex = 1.0;
    stiffeningParameter = 0.0;

    hyperbolicWaveSpeed = 2.0;
    parabolicDamping = sqrt(0.18);
}

MHDMaterialParameters::MHDMaterialParameters(double newAdiabaticIndex)
{
    adiabaticIndex = newAdiabaticIndex;
    stiffeningParameter = 0.0;

    hyperbolicWaveSpeed = 0.0;
    parabolicDamping = 0.0;
}

void MHDMaterialParameters::setAdiabaticIndex(double newAdiabaticIndex)
{
    adiabaticIndex = newAdiabaticIndex;
}

void MHDMaterialParameters::setStiffeningParameter(double newStiffeningParameter)
{
    stiffeningParameter = newStiffeningParameter;
}

void MHDMaterialParameters::setHyperbolicWaveSpeed(double newHyperbolicWaveSpeed)
{
    hyperbolicWaveSpeed = newHyperbolicWaveSpeed;
}

void MHDMaterialParameters::setParabolicDamping(double newParabolicDamping)
{
    parabolicDamping = newParabolicDamping;
}

double MHDMaterialParameters::getAdiabaticIndex()
{
    return adiabaticIndex;
}

double MHDMaterialParameters::getStiffeningParameter()
{
    return stiffeningParameter;
}

double MHDMaterialParameters::getHyperbolicWaveSpeed()
{
    return hyperbolicWaveSpeed;
}

double MHDMaterialParameters::getParabolicDamping()
{
    return parabolicDamping;
}

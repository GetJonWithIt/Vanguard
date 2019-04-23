#include "hprmiegruneisen.h"

HPRMieGruneisen::HPRMieGruneisen()
{
}

double HPRMieGruneisen::computeGruneisenCoefficient(double density, HPRMaterialParameters materialParameters)
{
    string equationOfState = materialParameters.getEquationOfState();

    if (equationOfState == "StiffenedGas")
    {
        double adiabaticIndex = materialParameters.getAdiabaticIndex();

        return adiabaticIndex - 1.0;
    }
    else if (equationOfState == "ShockMieGruneisen")
    {
        double referenceGruneisenCoefficient = materialParameters.getReferenceGruneisenCoefficient();
        double referenceDensity = materialParameters.getReferenceDensity();

        return referenceGruneisenCoefficient * (referenceDensity / density);
    }
    else if (equationOfState == "GodunovRomenski")
    {
        return materialParameters.getAdiabaticIndex();
    }
}

double HPRMieGruneisen::computeReferencePressure(double density, HPRMaterialParameters materialParameters)
{
    string equationOfState = materialParameters.getEquationOfState();

    if (equationOfState == "StiffenedGas")
    {
        double stiffeningParameter = materialParameters.getStiffeningParameter();

        return -stiffeningParameter;
    }
    else if (equationOfState == "ShockMieGruneisen")
    {
        double referenceSoundSpeedSquared = materialParameters.computeReferenceSoundSpeedSquared();
        double referenceDensity = materialParameters.getReferenceDensity();
        double hugoniotSlopeCoefficient = materialParameters.getHugoniotSlopeCoefficient();

        if (density > referenceDensity)
        {
            return referenceSoundSpeedSquared * (((1.0 / referenceDensity) - (1.0 / density)) / pow(((1.0 / referenceDensity) - (hugoniotSlopeCoefficient * ((1.0 / referenceDensity) -
                                                                                                                                                             (1.0 / density)))), 2.0));
        }
        else
        {
            return referenceSoundSpeedSquared * (density - referenceDensity);
        }
    }
    else if (equationOfState == "GodunovRomenski")
    {
        double referenceSoundSpeedSquared = materialParameters.computeReferenceSoundSpeedSquared();
        double alphaParameter = materialParameters.getAlphaParameter();
        double referenceDensity = materialParameters.getReferenceDensity();

        double coefficient = pow((density / referenceDensity), alphaParameter);

        return referenceSoundSpeedSquared * (density / alphaParameter) * (coefficient - 1.0) * coefficient;
    }
}

double HPRMieGruneisen::computeReferenceInternalEnergy(double density, HPRMaterialParameters materialParameters)
{
    string equationOfState = materialParameters.getEquationOfState();

    if (equationOfState == "StiffenedGas")
    {
        double stiffeningParameter = materialParameters.getStiffeningParameter();

        return stiffeningParameter / density;
    }
    else if (equationOfState == "ShockMieGruneisen")
    {
        double referenceDensity = materialParameters.getReferenceDensity();
        double referencePressure = materialParameters.getReferencePressure();

        if (density > referenceDensity)
        {
            return 0.5 * referencePressure * ((1.0 / referenceDensity) - (1.0 / density));
        }
        else
        {
            return 0.0;
        }
    }
    else if (equationOfState == "GodunovRomenski")
    {
        double referenceSoundSpeedSquared = materialParameters.computeReferenceSoundSpeedSquared();
        double alphaParameter = materialParameters.getAlphaParameter();
        double referenceDensity = materialParameters.getReferenceDensity();

        double coefficient = pow((density / referenceDensity), alphaParameter);

        return (referenceSoundSpeedSquared / (2.0 * (alphaParameter * alphaParameter))) * ((coefficient - 1.0) * (coefficient - 1.0));
    }
}

double HPRMieGruneisen::computeGruneisenCoefficientDerivative(double density, HPRMaterialParameters materialParameters)
{
    string equationOfState = materialParameters.getEquationOfState();

    if (equationOfState == "StiffenedGas")
    {
        return 0.0;
    }
    else if (equationOfState == "ShockMieGruneisen")
    {
        double referenceGruneisenCoefficient = materialParameters.getReferenceGruneisenCoefficient();
        double referenceDensity = materialParameters.getReferenceDensity();

        return -referenceGruneisenCoefficient * (referenceDensity / (density * density));
    }
}

double HPRMieGruneisen::computeReferencePressureDerivative(double density, HPRMaterialParameters materialParameters)
{
    string equationOfState = materialParameters.getEquationOfState();

    if (equationOfState == "StiffenedGas")
    {
        return 0.0;
    }
    else if (equationOfState == "ShockMieGruneisen")
    {
        double referenceSoundSpeedSquared = materialParameters.computeReferenceSoundSpeedSquared();
        double referenceDensity = materialParameters.getReferenceDensity();
        double hugoniotSlopeCoefficient = materialParameters.getHugoniotSlopeCoefficient();

        if (density > referenceDensity)
        {
            return referenceSoundSpeedSquared * (referenceDensity * referenceDensity) * (((hugoniotSlopeCoefficient * (referenceDensity - density)) - density) /
                                                                                         pow((hugoniotSlopeCoefficient * (density - referenceDensity)) - density, 3.0));
        }
        else
        {
            return referenceSoundSpeedSquared;
        }
    }
    else if (equationOfState == "GodunovRomenski")
    {
        double referenceDensity = materialParameters.getReferenceDensity();
        double referenceSoundSpeedSquared = materialParameters.computeReferenceSoundSpeedSquared();
        double alphaParameter = materialParameters.getAlphaParameter();

        double coefficient = pow((density / referenceDensity), alphaParameter);

        return (referenceSoundSpeedSquared / alphaParameter) * coefficient * (((1.0 + alphaParameter) * (coefficient - 1.0)) + (alphaParameter * coefficient));
    }
}

double HPRMieGruneisen::computeReferenceInternalEnergyDerivative(double density, HPRMaterialParameters materialParameters)
{
    string equationOfState = materialParameters.getEquationOfState();

    if (equationOfState == "StiffenedGas")
    {
        double stiffeningParameter = materialParameters.getStiffeningParameter();

        return -(stiffeningParameter / (density * density));
    }
    else if (equationOfState == "ShockMieGruneisen")
    {
        double referenceSoundSpeedSquared = materialParameters.computeReferenceSoundSpeedSquared();
        double referenceDensity = materialParameters.getReferenceDensity();
        double hugoniotSlopeCoefficient = materialParameters.getHugoniotSlopeCoefficient();

        if (density > referenceDensity)
        {
            return -((density - referenceDensity) * referenceDensity * referenceSoundSpeedSquared) / pow(((hugoniotSlopeCoefficient * (density - referenceDensity)) - density), 3.0);
        }
        else
        {
            return 0.0;
        }
    }
    else if (equationOfState == "GodunovRomenski")
    {
        double referenceSoundSpeedSquared = materialParameters.computeReferenceSoundSpeedSquared();
        double alphaParameter = materialParameters.getAlphaParameter();
        double referenceDensity = materialParameters.getReferenceDensity();

        double coefficient = pow((density / referenceDensity), alphaParameter);

        return (referenceSoundSpeedSquared / (density * alphaParameter)) * (coefficient - 1.0) * coefficient;
    }
}

double HPRMieGruneisen::computeInternalEnergy(double density, double pressure, HPRMaterialParameters materialParameters)
{
    double gruneisenCoefficient = computeGruneisenCoefficient(density, materialParameters);
    double referencePressure = computeReferencePressure(density, materialParameters);
    double referenceInternalEnergy = computeReferenceInternalEnergy(density, materialParameters);

    return referenceInternalEnergy + ((pressure - referencePressure) / (density * gruneisenCoefficient));
}

double HPRMieGruneisen::computePressure(double density, double internalEnergy, HPRMaterialParameters materialParameters)
{
    double gruneisenCoefficient = computeGruneisenCoefficient(density, materialParameters);
    double referencePressure = computeReferencePressure(density, materialParameters);
    double referenceInternalEnergy = computeReferenceInternalEnergy(density, materialParameters);

    return ((internalEnergy - referenceInternalEnergy) * density * gruneisenCoefficient) + referencePressure;
}

double HPRMieGruneisen::computeTemperature(double density, double pressure, HPRMaterialParameters materialParameters)
{
    double specificHeatCapacity = materialParameters.getSpecificHeatCapacity();
    double referenceTemperature = materialParameters.getReferenceTemperature();

    double gruneisenCoefficient = computeGruneisenCoefficient(density, materialParameters);
    double referencePressure = computeReferencePressure(density, materialParameters);

    return referenceTemperature + ((pressure - referencePressure) / (density * gruneisenCoefficient * specificHeatCapacity));
}

double HPRMieGruneisen::computeInternalEnergyDerivativeDensity(double density, double pressure, HPRMaterialParameters materialParameters)
{
    double gruneisenCoefficient = computeGruneisenCoefficient(density, materialParameters);
    double gruneisenCoefficientDerivative = computeGruneisenCoefficientDerivative(density, materialParameters);

    double referencePressure = computeReferencePressure(density, materialParameters);
    double referencePressureDerivative = computeReferencePressureDerivative(density, materialParameters);
    double referenceInternalEnergyDerivative = computeReferenceInternalEnergyDerivative(density, materialParameters);

    return referenceInternalEnergyDerivative - (((referencePressureDerivative * density * gruneisenCoefficient) + ((gruneisenCoefficient + (density * gruneisenCoefficientDerivative)) *
                                                                                                                   (pressure - referencePressure))) / ((density * gruneisenCoefficient) *
                                                                                                                                                       (density * gruneisenCoefficient)));
}

double HPRMieGruneisen::computeInternalEnergyDerivativePressure(double density, HPRMaterialParameters materialParameters)
{
    double gruneisenCoefficient = computeGruneisenCoefficient(density, materialParameters);

    return 1.0 / (density * gruneisenCoefficient);
}

double HPRMieGruneisen::computeTemperatureDerivativeDensity(double density, double pressure, HPRMaterialParameters materialParameters)
{
    double specificHeatCapacity = materialParameters.getSpecificHeatCapacity();

    double gruneisenCoefficient = computeGruneisenCoefficient(density, materialParameters);
    double gruneisenCoefficientDerivative = computeGruneisenCoefficientDerivative(density, materialParameters);

    double referencePressure = computeReferencePressure(density, materialParameters);
    double referencePressureDerivative = computeReferencePressureDerivative(density, materialParameters);

    return -((referencePressureDerivative * density * gruneisenCoefficient) + ((gruneisenCoefficient + (density * gruneisenCoefficientDerivative)) * (pressure - referencePressure))) /
            (((density * gruneisenCoefficient) * (density * gruneisenCoefficient)) / specificHeatCapacity);
}

double HPRMieGruneisen::computeTemperatureDerivativePressure(double density, HPRMaterialParameters materialParameters)
{
    double specificHeatCapacity = materialParameters.getSpecificHeatCapacity();

    return computeInternalEnergyDerivativePressure(density, materialParameters) / specificHeatCapacity;
}

#ifndef HPRMIEGRUNEISEN_H
#define HPRMIEGRUNEISEN_H

#include "hprmaterialparameters.h"
using namespace std;

class HPRMieGruneisen
{
public:
    HPRMieGruneisen();

    static double computeGruneisenCoefficient(double density, HPRMaterialParameters materialParameters);
    static double computeReferencePressure(double density, HPRMaterialParameters materialParameters);
    static double computeReferenceInternalEnergy(double density, HPRMaterialParameters materialParameters);

    static double computeGruneisenCoefficientDerivative(double density, HPRMaterialParameters materialParameters);
    static double computeReferencePressureDerivative(double density, HPRMaterialParameters materialParameters);
    static double computeReferenceInternalEnergyDerivative(double density, HPRMaterialParameters materialParameters);

    static double computeInternalEnergy(double density, double pressure, HPRMaterialParameters materialParameters);
    static double computePressure(double density, double internalEnergy, HPRMaterialParameters materialParameters);
    static double computeTemperature(double density, double pressure, HPRMaterialParameters materialParameters);

    static double computeInternalEnergyDerivativeDensity(double density, double pressure, HPRMaterialParameters materialParameters);
    static double computeInternalEnergyDerivativePressure(double density, HPRMaterialParameters materialParameters);

    static double computeTemperatureDerivativeDensity(double density, double pressure, HPRMaterialParameters materialParameters);
    static double computeTemperatureDerivativePressure(double density, HPRMaterialParameters materialParameters);
};

#endif // HPRMIEGRUNEISEN_H

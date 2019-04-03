#ifndef EULEREQUATIONOFSTATE_H
#define EULEREQUATIONOFSTATE_H

#include "eulermaterialparameters.h"
#include <cmath>
using namespace std;

class EulerEquationOfState
{
public:
    EulerEquationOfState();

    static double computeSpecificInternalEnergy(double density, double pressure, EulerMaterialParameters materialParameters);
    static double computePressure(double density, double specificInternalEnergy, EulerMaterialParameters materialParameters);
    static double computeSoundSpeed(double density, double pressure, EulerMaterialParameters materialParameters);
    static double computeEntropy(double density, double pressure, EulerMaterialParameters materialParameters);

    static double computeTemperature(double specificInternalEnergy, EulerMaterialParameters materialParameters);
    static double computeReactionRate(double specificInternalEnergy, EulerMaterialParameters materialParameters);
};

#endif // EULEREQUATIONOFSTATE_H

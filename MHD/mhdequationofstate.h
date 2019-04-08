#ifndef MHDEQUATIONOFSTATE_H
#define MHDEQUATIONOFSTATE_H

#include "mhdmaterialparameters.h"
#include <cmath>
using namespace std;

class MHDEquationOfState
{
public:
    MHDEquationOfState();

    static double computeSpecificInternalEnergy(double density, double pressure, MHDMaterialParameters materialParameters);
    static double computePressure(double density, double specificInternalEnergy, MHDMaterialParameters materialParameters);
    static double computeSoundSpeed(double density, double pressure, MHDMaterialParameters materialParameters);
    static double computeEntropy(double density, double pressure, MHDMaterialParameters materialParameters);
};

#endif // MHDEQUATIONOFSTATE_H

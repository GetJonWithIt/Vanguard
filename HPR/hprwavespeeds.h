#ifndef HPRWAVESPEEDS_H
#define HPRWAVESPEEDS_H

#include "hprstatevector.h"
using namespace std;

class HPRWaveSpeeds
{
public:
    HPRWaveSpeeds();

    static double computeAdiabaticSoundSpeed(double density, double pressure, HPRMaterialParameters materialParameters);
    static double computeHeatCharacteristicSpeed(double density, double temperature, HPRMaterialParameters materialParameters);
};

#endif // HPRWAVESPEEDS_H

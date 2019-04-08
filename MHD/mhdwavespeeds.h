#ifndef MHDWAVESPEEDS_H
#define MHDWAVESPEEDS_H

#include "mhdequationofstate.h"
#include "Mathematics/vectoralgebra.h"
using namespace std;

class MHDWaveSpeeds
{
public:
    MHDWaveSpeeds();

    static double computeAlfvenWaveSpeed(double density, double xMagneticField, double yMagneticField, double zMagneticField);

    static double computeSlowMagnetoAcousticSpeed(double density, double pressure, double xMagneticField, double yMagneticField, double zMagneticField, MHDMaterialParameters materialParameters);
    static double computeFastMagnetoAcousticSpeed(double density, double pressure, double xMagneticField, double yMagneticField, double zMagneticField, MHDMaterialParameters materialParameters);
};

#endif // MHDWAVESPEEDS_H

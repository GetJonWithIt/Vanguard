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

    static double computeXSlowMagnetoAcousticSpeed(double density, double pressure, double xMagneticField, double yMagneticField, double zMagneticField, MHDMaterialParameters materialParameters);
    static double computeYSlowMagnetoAcousticSpeed(double density, double pressure, double xMagneticField, double yMagneticField, double zMagneticField, MHDMaterialParameters materialParameters);

    static double computeXFastMagnetoAcousticSpeed(double density, double pressure, double xMagneticField, double yMagneticField, double zMagneticField, MHDMaterialParameters materialParameters);
    static double computeYFastMagnetoAcousticSpeed(double density, double pressure, double xMagneticField, double yMagneticField, double zMagneticField, MHDMaterialParameters materialParameters);
};

#endif // MHDWAVESPEEDS_H

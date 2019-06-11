#include "mhdwavespeeds.h"

MHDWaveSpeeds::MHDWaveSpeeds()
{
}

double MHDWaveSpeeds::computeAlfvenWaveSpeed(double density, double xMagneticField, double yMagneticField, double zMagneticField)
{
    vector<double> magneticFieldVector(3);
    magneticFieldVector[0] = xMagneticField;
    magneticFieldVector[1] = yMagneticField;
    magneticFieldVector[2] = zMagneticField;

    return (VectorAlgebra::computeNorm(magneticFieldVector) / sqrt(density));
}

double MHDWaveSpeeds::computeXSlowMagnetoAcousticSpeed(double density, double pressure, double xMagneticField, double yMagneticField, double zMagneticField, MHDMaterialParameters materialParameters)
{
    double soundSpeed = MHDEquationOfState::computeSoundSpeed(density, pressure, materialParameters);
    double alfvenWaveSpeed = computeAlfvenWaveSpeed(density, xMagneticField, yMagneticField, zMagneticField);

    double soundSpeedSquared = (soundSpeed * soundSpeed);
    double alfvenWaveSpeedSquared = (alfvenWaveSpeed * alfvenWaveSpeed);

    return sqrt(0.5 * (soundSpeedSquared + alfvenWaveSpeedSquared - sqrt(((soundSpeedSquared + alfvenWaveSpeedSquared) * (soundSpeedSquared + alfvenWaveSpeedSquared)) -
                                                                         4.0 * ((soundSpeedSquared * (xMagneticField * xMagneticField)) / density))));
}

double MHDWaveSpeeds::computeYSlowMagnetoAcousticSpeed(double density, double pressure, double xMagneticField, double yMagneticField, double zMagneticField, MHDMaterialParameters materialParameters)
{
    double soundSpeed = MHDEquationOfState::computeSoundSpeed(density, pressure, materialParameters);
    double alfvenWaveSpeed = computeAlfvenWaveSpeed(density, xMagneticField, yMagneticField, zMagneticField);

    double soundSpeedSquared = (soundSpeed * soundSpeed);
    double alfvenWaveSpeedSquared = (alfvenWaveSpeed * alfvenWaveSpeed);

    return sqrt(0.5 * (soundSpeedSquared + alfvenWaveSpeedSquared - sqrt(((soundSpeedSquared + alfvenWaveSpeedSquared) * (soundSpeedSquared + alfvenWaveSpeedSquared)) -
                                                                         4.0 * ((soundSpeedSquared * (yMagneticField * yMagneticField)) / density))));
}

double MHDWaveSpeeds::computeXFastMagnetoAcousticSpeed(double density, double pressure, double xMagneticField, double yMagneticField, double zMagneticField, MHDMaterialParameters materialParameters)
{
    double soundSpeed = MHDEquationOfState::computeSoundSpeed(density, pressure, materialParameters);
    double alfvenWaveSpeed = computeAlfvenWaveSpeed(density, xMagneticField, yMagneticField, zMagneticField);

    double soundSpeedSquared = (soundSpeed * soundSpeed);
    double alfvenWaveSpeedSquared = (alfvenWaveSpeed * alfvenWaveSpeed);

    return sqrt(0.5 * (soundSpeedSquared + alfvenWaveSpeedSquared + sqrt(((soundSpeedSquared + alfvenWaveSpeedSquared) * (soundSpeedSquared + alfvenWaveSpeedSquared)) -
                                                                         4.0 * ((soundSpeedSquared * (xMagneticField * xMagneticField)) / density))));
}

double MHDWaveSpeeds::computeYFastMagnetoAcousticSpeed(double density, double pressure, double xMagneticField, double yMagneticField, double zMagneticField, MHDMaterialParameters materialParameters)
{
    double soundSpeed = MHDEquationOfState::computeSoundSpeed(density, pressure, materialParameters);
    double alfvenWaveSpeed = computeAlfvenWaveSpeed(density, xMagneticField, yMagneticField, zMagneticField);

    double soundSpeedSquared = (soundSpeed * soundSpeed);
    double alfvenWaveSpeedSquared = (alfvenWaveSpeed * alfvenWaveSpeed);

    return sqrt(0.5 * (soundSpeedSquared + alfvenWaveSpeedSquared + sqrt(((soundSpeedSquared + alfvenWaveSpeedSquared) * (soundSpeedSquared + alfvenWaveSpeedSquared)) -
                                                                         4.0 * ((soundSpeedSquared * (yMagneticField * yMagneticField)) / density))));
}

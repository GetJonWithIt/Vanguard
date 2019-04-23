#include "hpracoustictensor.h"

HPRAcousticTensor::HPRAcousticTensor()
{
}

vector<vector<double> > HPRAcousticTensor::computeAcousticTensorComponent1(HPRStateVector stateVector, HPRMaterialParameters materialParameters, int direction)
{
    vector<vector<double> > acousticTensorComponent1(4, vector<double>(5));

    double density = stateVector.getDensity();
    acousticTensorComponent1[0][1] = 1.0 / density;

    vector<vector<double> > shearStressTensorDerivativeDensity = stateVector.computeShearStressTensorDerivativeDensity(materialParameters);
    vector<vector<vector<vector<double> > > > shearStressTensorDerivativeDistortionTensor = stateVector.computeShearStressTensorDerivativeDistortionTensor(materialParameters);

    for (int i = 0; i < 3; i++)
    {
        acousticTensorComponent1[i][0] = -(1.0 / density) * shearStressTensorDerivativeDensity[direction][i];
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            acousticTensorComponent1[i][2 + j] = -(1.0 / density) * shearStressTensorDerivativeDistortionTensor[direction][i][j][direction];
        }
    }

    if (materialParameters.getIsThermal())
    {
        double temperatureDerivativeDensity = stateVector.computeTemperatureDerivativeDensity(materialParameters);
        double temperatureDerivativePressure = stateVector.computeTemperatureDerivativePressure(materialParameters);

        acousticTensorComponent1[3][0] = temperatureDerivativeDensity / density;
        acousticTensorComponent1[3][1] = temperatureDerivativePressure / density;
    }

    return acousticTensorComponent1;
}

vector<vector<double> > HPRAcousticTensor::computeAcousticTensorComponent2(HPRStateVector stateVector, HPRMaterialParameters materialParameters, int direction)
{
    vector<vector<double> > acousticTensorComponent2(5, vector<double>(4));

    double density = stateVector.getDensity();
    double pressure = stateVector.getPressure();
    vector<vector<double> > distortionTensor = stateVector.getDistortionTensor();

    double adiabaticSoundSpeed = HPRWaveSpeeds::computeAdiabaticSoundSpeed(density, pressure, materialParameters);

    acousticTensorComponent2[0][0] = density;
    acousticTensorComponent2[1][direction] = density * (adiabaticSoundSpeed * adiabaticSoundSpeed);

    vector<vector<double> > shearStressTensor = stateVector.computeShearStressTensor(materialParameters);
    vector<vector<double> > shearStressTensorDerivativeDensity = stateVector.computeShearStressTensorDerivativeDensity(materialParameters);

    for (int i = 0; i < 3; i++)
    {
        acousticTensorComponent2[1][i] += shearStressTensor[direction][i] - (density * shearStressTensorDerivativeDensity[direction][i]);
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            acousticTensorComponent2[2 + i][j] = distortionTensor[i][j];
        }
    }

    if (materialParameters.getIsThermal())
    {
        double temperature = stateVector.computeTemperature(materialParameters);
        double temperatureDerivativePressure = stateVector.computeTemperatureDerivativePressure(materialParameters);

        double heatCharacteristicSpeed = HPRWaveSpeeds::computeHeatCharacteristicSpeed(density, temperature, materialParameters);

        acousticTensorComponent2[1][3] = density * ((heatCharacteristicSpeed * heatCharacteristicSpeed) / temperatureDerivativePressure);
    }

    return acousticTensorComponent2;
}

vector<vector<double> > HPRAcousticTensor::computeAcousticTensor(HPRStateVector stateVector, HPRMaterialParameters materialParameters, int direction)
{
    vector<vector<double> > acousticTensorComponent1 = computeAcousticTensorComponent1(stateVector, materialParameters, direction);
    vector<vector<double> > acousticTensorComponent2 = computeAcousticTensorComponent2(stateVector, materialParameters, direction);

    return MatrixAlgebra::multiplyMatrices(acousticTensorComponent1, acousticTensorComponent2);
}

double HPRAcousticTensor::computeMaximumWaveSpeed(HPRStateVector stateVector, HPRMaterialParameters materialParameters, int direction)
{
    vector<vector<double> > acousticTensor = computeAcousticTensor(stateVector, materialParameters, direction);

    return sqrt(MatrixAlgebra::computeLargestEigenvalue(acousticTensor));
}

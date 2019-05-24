#include "hprintermediateacoustictensor.h"

HPRIntermediateAcousticTensor::HPRIntermediateAcousticTensor()
{
}

vector<vector<double> > HPRIntermediateAcousticTensor::computeMaterial1AcousticTensorComponent1(HPRIntermediateStateVector stateVector, HPRMaterialParameters material1Parameters, int direction)
{
    vector<vector<double> > material1AcousticTensorComponent1(4, vector<double>(5));

    double material1Density = stateVector.getMaterial1Density();
    material1AcousticTensorComponent1[0][1] = 1.0 / material1Density;

    vector<vector<double> > material1ShearStressTensorDerivativeDensity = stateVector.computeShearStressTensorDerivativeDensity(material1Parameters);
    vector<vector<vector<vector<double> > > > material1ShearStressTensorDerivativeDistortionTensor = stateVector.computeMaterial1ShearStressTensorDerivativeDistortionTensor(material1Parameters);

    for (int i = 0; i < 3; i++)
    {
        material1AcousticTensorComponent1[i][0] = -(1.0 / material1Density) * material1ShearStressTensorDerivativeDensity[direction][i];
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            material1AcousticTensorComponent1[i][2 + j] = -(1.0 / material1Density) * material1ShearStressTensorDerivativeDistortionTensor[direction][i][j][direction];
        }
    }

    if (material1Parameters.getIsThermal())
    {
        double material1TemperatureDerivativeDensity = stateVector.computeMaterial1TemperatureDerivativeDensity(material1Parameters);
        double material1TemperatureDerivativePressure = stateVector.computeMaterial1TemperatureDerivativePressure(material1Parameters);

        material1AcousticTensorComponent1[3][0] = material1TemperatureDerivativeDensity / material1Density;
        material1AcousticTensorComponent1[3][1] = material1TemperatureDerivativePressure / material1Density;
    }

    return material1AcousticTensorComponent1;
}

vector<vector<double> > HPRIntermediateAcousticTensor::computeMaterial2AcousticTensorComponent1(HPRIntermediateStateVector stateVector, HPRMaterialParameters material2Parameters, int direction)
{
    vector<vector<double> > material2AcousticTensorComponent1(4, vector<double>(5));

    double material2Density = stateVector.getMaterial2Density();
    material2AcousticTensorComponent1[0][1] = 1.0 / material2Density;

    vector<vector<double> > material2ShearStressTensorDerivativeDensity = stateVector.computeShearStressTensorDerivativeDensity(material2Parameters);
    vector<vector<vector<vector<double> > > > material2ShearStressTensorDerivativeDistortionTensor = stateVector.computeMaterial2ShearStressTensorDerivativeDistortionTensor(material2Parameters);

    for (int i = 0; i < 3; i++)
    {
        material2AcousticTensorComponent1[i][0] = -(1.0 / material2Density) * material2ShearStressTensorDerivativeDensity[direction][i];
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            material2AcousticTensorComponent1[i][2 + j] = -(1.0 / material2Density) * material2ShearStressTensorDerivativeDistortionTensor[direction][i][j][direction];
        }
    }

    if (material2Parameters.getIsThermal())
    {
        double material2TemperatureDerivativeDensity = stateVector.computeMaterial2TemperatureDerivativeDensity(material2Parameters);
        double material2TemperatureDerivativePressure = stateVector.computeMaterial2TemperatureDerivativePressure(material2Parameters);

        material2AcousticTensorComponent1[3][0] = material2TemperatureDerivativeDensity / material2Density;
        material2AcousticTensorComponent1[3][1] = material2TemperatureDerivativePressure / material2Density;
    }

    return material2AcousticTensorComponent1;
}

vector<vector<double> > HPRIntermediateAcousticTensor::computeMaterial1AcousticTensorComponent2(HPRIntermediateStateVector stateVector, HPRMaterialParameters material1Parameters, int direction)
{
    vector<vector<double> > material1AcousticTensorComponent2(5, vector<double>(4));

    double material1Density = stateVector.getMaterial1Density();
    double material1Pressure = stateVector.getMaterial1Pressure();
    vector<vector<double> > interfaceDistortionTensor = stateVector.getInterfaceDistortionTensor();

    double material1AdiabaticSoundSpeed = HPRWaveSpeeds::computeAdiabaticSoundSpeed(material1Density, material1Pressure, material1Parameters);

    material1AcousticTensorComponent2[0][0] = material1Density;
    material1AcousticTensorComponent2[1][direction] = material1Density * (material1AdiabaticSoundSpeed * material1AdiabaticSoundSpeed);

    vector<vector<double> > material1ShearStressTensor = stateVector.computeMaterial1ShearStressTensor(material1Parameters);
    vector<vector<double> > material1ShearStressTensorDerivativeDensity = stateVector.computeShearStressTensorDerivativeDensity(material1Parameters);

    for (int i = 0; i < 3; i++)
    {
        material1AcousticTensorComponent2[1][i] += material1ShearStressTensor[direction][i] - (material1Density * material1ShearStressTensorDerivativeDensity[direction][i]);
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            material1AcousticTensorComponent2[2 + i][j] = interfaceDistortionTensor[i][j];
        }
    }

    if (material1Parameters.getIsThermal())
    {
        double material1Temperature = stateVector.computeMaterial1Temperature(material1Parameters);
        double material1TemperatureDerivativePressure = stateVector.computeMaterial1TemperatureDerivativePressure(material1Parameters);

        double material1HeatCharacteristicSpeed = HPRWaveSpeeds::computeHeatCharacteristicSpeed(material1Density, material1Temperature, material1Parameters);

        material1AcousticTensorComponent2[1][3] = material1Density * ((material1HeatCharacteristicSpeed * material1HeatCharacteristicSpeed) / material1TemperatureDerivativePressure);
    }

    return material1AcousticTensorComponent2;
}

vector<vector<double> > HPRIntermediateAcousticTensor::computeMaterial2AcousticTensorComponent2(HPRIntermediateStateVector stateVector, HPRMaterialParameters material2Parameters, int direction)
{
    vector<vector<double> > material2AcousticTensorComponent2(5, vector<double>(4));

    double material2Density = stateVector.getMaterial2Density();
    double material2Pressure = stateVector.getMaterial2Pressure();
    vector<vector<double> > interfaceDistortionTensor = stateVector.getInterfaceDistortionTensor();

    double material2AdiabaticSoundSpeed = HPRWaveSpeeds::computeAdiabaticSoundSpeed(material2Density, material2Pressure, material2Parameters);

    material2AcousticTensorComponent2[0][0] = material2Density;
    material2AcousticTensorComponent2[1][direction] = material2Density * (material2AdiabaticSoundSpeed * material2AdiabaticSoundSpeed);

    vector<vector<double> > material2ShearStressTensor = stateVector.computeMaterial2ShearStressTensor(material2Parameters);
    vector<vector<double> > material2ShearStressTensorDerivativeDensity = stateVector.computeShearStressTensorDerivativeDensity(material2Parameters);

    for (int i = 0; i < 3; i++)
    {
        material2AcousticTensorComponent2[1][i] += material2ShearStressTensor[direction][i] - (material2Density * material2ShearStressTensorDerivativeDensity[direction][i]);
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            material2AcousticTensorComponent2[2 + i][j] = interfaceDistortionTensor[i][j];
        }
    }

    if (material2Parameters.getIsThermal())
    {
        double material2Temperature = stateVector.computeMaterial2Temperature(material2Parameters);
        double material2TemperatureDerivativePressure = stateVector.computeMaterial2TemperatureDerivativePressure(material2Parameters);

        double material2HeatCharacteristicSpeed = HPRWaveSpeeds::computeHeatCharacteristicSpeed(material2Density, material2Temperature, material2Parameters);

        material2AcousticTensorComponent2[1][3] = material2Density * ((material2HeatCharacteristicSpeed * material2HeatCharacteristicSpeed) / material2TemperatureDerivativePressure);
    }

    return material2AcousticTensorComponent2;
}

vector<vector<double> > HPRIntermediateAcousticTensor::computeMaterial1AcousticTensor(HPRIntermediateStateVector stateVector, HPRMaterialParameters material1Parameters, int direction)
{
    vector<vector<double> > material1AcousticTensorComponent1 = computeMaterial1AcousticTensorComponent1(stateVector, material1Parameters, direction);
    vector<vector<double> > material1AcousticTensorComponent2 = computeMaterial1AcousticTensorComponent2(stateVector, material1Parameters, direction);

    return MatrixAlgebra::multiplyMatrices(material1AcousticTensorComponent1, material1AcousticTensorComponent2);
}

vector<vector<double> > HPRIntermediateAcousticTensor::computeMaterial2AcousticTensor(HPRIntermediateStateVector stateVector, HPRMaterialParameters material2Parameters, int direction)
{
    vector<vector<double> > material2AcousticTensorComponent1 = computeMaterial2AcousticTensorComponent1(stateVector, material2Parameters, direction);
    vector<vector<double> > material2AcousticTensorComponent2 = computeMaterial2AcousticTensorComponent2(stateVector, material2Parameters, direction);

    return MatrixAlgebra::multiplyMatrices(material2AcousticTensorComponent1, material2AcousticTensorComponent2);
}

double HPRIntermediateAcousticTensor::computeMaximumWaveSpeed(HPRIntermediateStateVector stateVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters, int direction)
{
    vector<vector<double> > material1AcousticTensor = computeMaterial1AcousticTensor(stateVector, material1Parameters, direction);
    vector<vector<double> > material2AcousticTensor = computeMaterial2AcousticTensor(stateVector, material2Parameters, direction);

    return max(sqrt(MatrixAlgebra::computeLargestEigenvalue(material1AcousticTensor)), sqrt(MatrixAlgebra::computeLargestEigenvalue(material2AcousticTensor)));
}

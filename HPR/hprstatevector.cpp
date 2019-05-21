#include "hprstatevector.h"

HPRStateVector::HPRStateVector()
{
    density = 1.0;
    pressure = 1.0 / 1.4;
    distortionTensor = MatrixAlgebra::multiplyMatrix(pow(density, (1.0 / 3.0)), MatrixAlgebra::computeIdentityMatrix(3));

    xVelocity = 0.0;
    yVelocity = 0.0;
    zVelocity = 0.0;

    xThermalImpulse = 0.0;
    yThermalImpulse = 0.0;
    zThermalImpulse = 0.0;
}

HPRStateVector::HPRStateVector(double newDensity, double newPressure, vector<vector<double> > newDistortionTensor, double newXVelocity, double newYVelocity, double newZVelocity,
                               double newXThermalImpulse, double newYThermalImpulse, double newZThermalImpulse)
{
    density = newDensity;
    pressure = newPressure;
    distortionTensor = newDistortionTensor;

    xVelocity = newXVelocity;
    yVelocity = newYVelocity;
    zVelocity = newZVelocity;

    xThermalImpulse = newXThermalImpulse;
    yThermalImpulse = newYThermalImpulse;
    zThermalImpulse = newZThermalImpulse;
}

HPRStateVector::HPRStateVector(double newXVelocity, double newYVelocity, double newZVelocity, vector<vector<double> > newDistortionTensor, double newEntropy,
                               HyperelasticMaterialParameters hyperelasticMaterialParameters, HPRMaterialParameters materialParameters)
{
    double referenceMassDensity = hyperelasticMaterialParameters.getReferenceMassDensity();
    density = referenceMassDensity * MatrixAlgebra::computeDeterminant(newDistortionTensor);
    distortionTensor = newDistortionTensor;

    xVelocity = newXVelocity;
    yVelocity = newYVelocity;
    zVelocity = newZVelocity;

    xThermalImpulse = 0.0;
    yThermalImpulse = 0.0;
    zThermalImpulse = 0.0;

    double totalEnergy = ElasticEquationOfState::computeTotalEnergy(distortionTensor, newEntropy, xVelocity, yVelocity, zVelocity, hyperelasticMaterialParameters);
    pressure = HPREquationOfState::computePressure(density, totalEnergy, xVelocity, yVelocity, zVelocity, distortionTensor, xThermalImpulse, yThermalImpulse, zThermalImpulse, materialParameters);
}

void HPRStateVector::setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector, HPRMaterialParameters materialParameters)
{
    density = newPrimitiveVariableVector[0];
    pressure = newPrimitiveVariableVector[1];

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            distortionTensor[i][j] = newPrimitiveVariableVector[2 + (i * 3) + j];
        }
    }

    xVelocity = newPrimitiveVariableVector[11];
    yVelocity = newPrimitiveVariableVector[12];
    zVelocity = newPrimitiveVariableVector[13];

    if (materialParameters.getIsThermal())
    {
        xThermalImpulse = newPrimitiveVariableVector[14];
        yThermalImpulse = newPrimitiveVariableVector[15];
        zThermalImpulse = newPrimitiveVariableVector[16];
    }
}

void HPRStateVector::setConservedVariableVector(vector<double> newConservedVariableVector, HPRMaterialParameters materialParameters)
{
    density = newConservedVariableVector[0];

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            distortionTensor[i][j] = newConservedVariableVector[2 + (i * 3) + j];
        }
    }

    xVelocity = newConservedVariableVector[11] / density;
    yVelocity = newConservedVariableVector[12] / density;
    zVelocity = newConservedVariableVector[13] / density;

    if (materialParameters.getIsThermal())
    {
        xThermalImpulse = newConservedVariableVector[14] / density;
        yThermalImpulse = newConservedVariableVector[15] / density;
        zThermalImpulse = newConservedVariableVector[16] / density;
    }

    double totalEnergy = newConservedVariableVector[1] / density;
    pressure = HPREquationOfState::computePressure(density, totalEnergy, xVelocity, yVelocity, zVelocity, distortionTensor, xThermalImpulse, yThermalImpulse, zThermalImpulse, materialParameters);
}

vector<double> HPRStateVector::computePrimitiveVariableVector(HPRMaterialParameters materialParameters)
{
    vector<double> primitiveVariableVector(17);

    primitiveVariableVector[0] = density;
    primitiveVariableVector[1] = pressure;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            primitiveVariableVector[2 + (i * 3) + j] = distortionTensor[i][j];
        }
    }

    primitiveVariableVector[11] = xVelocity;
    primitiveVariableVector[12] = yVelocity;
    primitiveVariableVector[13] = zVelocity;

    if (materialParameters.getIsThermal())
    {
        primitiveVariableVector[14] = xThermalImpulse;
        primitiveVariableVector[15] = yThermalImpulse;
        primitiveVariableVector[16] = zThermalImpulse;
    }

    return primitiveVariableVector;
}

vector<double> HPRStateVector::computeConservedVariableVector(HPRMaterialParameters materialParameters)
{
    vector<double> conservedVariableVector(17);

    conservedVariableVector[0] = density;
    conservedVariableVector[1] = density * HPREquationOfState::computeTotalEnergy(density, pressure, xVelocity, yVelocity, zVelocity, distortionTensor, xThermalImpulse, yThermalImpulse,
                                                                                  zThermalImpulse, materialParameters);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            conservedVariableVector[2 + (i * 3) + j] = distortionTensor[i][j];
        }
    }

    conservedVariableVector[11] = density * xVelocity;
    conservedVariableVector[12] = density * yVelocity;
    conservedVariableVector[13] = density * zVelocity;

    if (materialParameters.getIsThermal())
    {
        conservedVariableVector[14] = density * xThermalImpulse;
        conservedVariableVector[15] = density * yThermalImpulse;
        conservedVariableVector[16] = density * zThermalImpulse;
    }

    return conservedVariableVector;
}

vector<double> HPRStateVector::computeXFluxVector(vector<double> conservedVariableVector, HPRMaterialParameters materialParameters)
{
    vector<double> fluxVector(17);

    double computedDensity = conservedVariableVector[0];

    vector<vector<double> > computedDistortionTensor(3, vector<double>(3));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            computedDistortionTensor[i][j] = conservedVariableVector[2 + (i * 3) + j];
        }
    }

    vector<vector<double> > computedShearStressTensor = HPREquationOfState::computeShearStressTensor(computedDensity, computedDistortionTensor, materialParameters);

    double computedXVelocity = conservedVariableVector[11] / computedDensity;
    double computedYVelocity = conservedVariableVector[12] / computedDensity;
    double computedZVelocity = conservedVariableVector[13] / computedDensity;

    double computedXThermalImpulse = conservedVariableVector[14] / computedDensity;
    double computedYThermalImpulse = conservedVariableVector[15] / computedDensity;
    double computedZThermalImpulse = conservedVariableVector[16] / computedDensity;

    double computedTotalEnergy = conservedVariableVector[1] / computedDensity;
    double computedPressure = HPREquationOfState::computePressure(computedDensity, computedTotalEnergy, computedXVelocity, computedYVelocity, computedZVelocity, computedDistortionTensor,
                                                                  computedXThermalImpulse, computedYThermalImpulse, computedZThermalImpulse, materialParameters);

    fluxVector[0] = computedDensity * computedXVelocity;
    fluxVector[1] = (computedDensity * (computedXVelocity * computedTotalEnergy)) + (computedPressure * computedXVelocity);

    fluxVector[11] = (computedDensity * (computedXVelocity * computedXVelocity)) + computedPressure;
    fluxVector[12] = computedDensity * (computedXVelocity * computedYVelocity);
    fluxVector[13] = computedDensity * (computedXVelocity * computedZVelocity);

    vector<double> computedVelocityVector(3);
    computedVelocityVector[0] = computedXVelocity;
    computedVelocityVector[1] = computedYVelocity;
    computedVelocityVector[2] = computedZVelocity;
    vector<double> distortionTensorVelocityVectorProduct = MatrixAlgebra::multiplyMatrixByVector(computedDistortionTensor, computedVelocityVector);

    fluxVector[1] -= VectorAlgebra::computeDotProduct(computedShearStressTensor[0], computedVelocityVector);
    for (int i = 0; i < 3; i++)
    {
        fluxVector[11 + i] -= computedShearStressTensor[0][i];
    }

    for (int i = 0; i < 3; i++)
    {
        fluxVector[2 + (i * 3)] = distortionTensorVelocityVectorProduct[i];
    }

    if (materialParameters.getIsThermal())
    {
        double computedTemperature = HPREquationOfState::computeTemperature(computedDensity, computedPressure, materialParameters);
        vector<double> computedHeatFluxVector = HPREquationOfState::computeHeatFluxVector(computedTemperature, computedXThermalImpulse, computedYThermalImpulse, computedZThermalImpulse,
                                                                                          materialParameters);

        fluxVector[1] += computedHeatFluxVector[0];

        fluxVector[14] = (computedDensity * (computedXVelocity * computedXThermalImpulse)) + computedTemperature;
        fluxVector[15] = computedDensity * (computedXVelocity * computedYThermalImpulse);
        fluxVector[16] = computedDensity * (computedXVelocity * computedZThermalImpulse);
    }

    return fluxVector;
}

vector<double> HPRStateVector::computeXFluxVector(HPRMaterialParameters materialParameters)
{
    return computeXFluxVector(computeConservedVariableVector(materialParameters), materialParameters);
}

vector<double> HPRStateVector::computeYFluxVector(vector<double> conservedVariableVector, HPRMaterialParameters materialParameters)
{
    vector<double> fluxVector(17);

    double computedDensity = conservedVariableVector[0];

    vector<vector<double> > computedDistortionTensor(3, vector<double>(3));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            computedDistortionTensor[i][j] = conservedVariableVector[2 + (i * 3) + j];
        }
    }

    vector<vector<double> > computedShearStressTensor = HPREquationOfState::computeShearStressTensor(computedDensity, computedDistortionTensor, materialParameters);

    double computedXVelocity = conservedVariableVector[11] / computedDensity;
    double computedYVelocity = conservedVariableVector[12] / computedDensity;
    double computedZVelocity = conservedVariableVector[13] / computedDensity;

    double computedXThermalImpulse = conservedVariableVector[14] / computedDensity;
    double computedYThermalImpulse = conservedVariableVector[15] / computedDensity;
    double computedZThermalImpulse = conservedVariableVector[16] / computedDensity;

    double computedTotalEnergy = conservedVariableVector[1] / computedDensity;
    double computedPressure = HPREquationOfState::computePressure(computedDensity, computedTotalEnergy, computedXVelocity, computedYVelocity, computedZVelocity, computedDistortionTensor,
                                                                  computedXThermalImpulse, computedYThermalImpulse, computedZThermalImpulse, materialParameters);

    fluxVector[0] = computedDensity * computedYVelocity;
    fluxVector[1] = (computedDensity * (computedYVelocity * computedTotalEnergy)) + (computedPressure * computedYVelocity);

    fluxVector[11] = computedDensity * (computedYVelocity * computedXVelocity);
    fluxVector[12] = (computedDensity * (computedYVelocity * computedYVelocity)) + computedPressure;
    fluxVector[13] = computedDensity * (computedYVelocity * computedZVelocity);

    vector<double> computedVelocityVector(3);
    computedVelocityVector[0] = computedXVelocity;
    computedVelocityVector[1] = computedYVelocity;
    computedVelocityVector[2] = computedZVelocity;
    vector<double> computedDistortionTensorVelocityVectorProduct = MatrixAlgebra::multiplyMatrixByVector(computedDistortionTensor, computedVelocityVector);

    fluxVector[1] -= VectorAlgebra::computeDotProduct(computedShearStressTensor[1], computedVelocityVector);
    for (int i = 0; i < 3; i++)
    {
        fluxVector[11 + i] -= computedShearStressTensor[1][i];
    }

    for (int i = 0; i < 3; i++)
    {
        fluxVector[2 + (i * 3) + 1] = computedDistortionTensorVelocityVectorProduct[i];
    }

    if (materialParameters.getIsThermal())
    {
        double computedTemperature = HPREquationOfState::computeTemperature(computedDensity, computedPressure, materialParameters);
        vector<double> computedHeatFluxVector = HPREquationOfState::computeHeatFluxVector(computedTemperature, computedXThermalImpulse, computedYThermalImpulse, computedZThermalImpulse,
                                                                                          materialParameters);

        fluxVector[1] += computedHeatFluxVector[1];

        fluxVector[14] = computedDensity * (computedYVelocity * computedXThermalImpulse);
        fluxVector[15] = (computedDensity * (computedYVelocity * computedYThermalImpulse)) + computedTemperature;
        fluxVector[16] = computedDensity * (computedYVelocity * computedZThermalImpulse);
    }

    return fluxVector;
}

vector<double> HPRStateVector::computeYFluxVector(HPRMaterialParameters materialParameters)
{
    return computeYFluxVector(computeConservedVariableVector(materialParameters), materialParameters);
}

vector<double> HPRStateVector::computeSourceTermVector(vector<double> conservedVariableVector, HPRMaterialParameters materialParameters)
{
    vector<double> sourceTermVector(17);

    double computedDensity = conservedVariableVector[0];

    vector<vector<double> > computedDistortionTensor(3, vector<double>(3));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            computedDistortionTensor[i][j] = conservedVariableVector[2 + (i * 3) + j];
        }
    }

    vector<vector<double> > totalEnergyDerivativeDistortionTensor = HPRDerivatives::computeTotalEnergyDerivativeDistortionTensor(computedDistortionTensor, materialParameters);
    double theta1Reciprocal = HPRSourceTerms::computeTheta1Reciprocal(computedDensity, computedDistortionTensor, materialParameters);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sourceTermVector[2 + (i * 3) + j] = -totalEnergyDerivativeDistortionTensor[i][j] * theta1Reciprocal;
        }
    }

    if (materialParameters.getIsThermal())
    {
        double computedTotalEnergy = conservedVariableVector[1] / computedDensity;

        double computedXVelocity = conservedVariableVector[11] / computedDensity;
        double computedYVelocity = conservedVariableVector[12] / computedDensity;
        double computedZVelocity = conservedVariableVector[13] / computedDensity;

        double computedXThermalImpulse = conservedVariableVector[14] / computedDensity;
        double computedYThermalImpulse = conservedVariableVector[15] / computedDensity;
        double computedZThermalImpulse = conservedVariableVector[16] / computedDensity;

        double computedPressure = HPREquationOfState::computePressure(computedDensity, computedTotalEnergy, computedXVelocity, computedYVelocity, computedZVelocity, computedDistortionTensor,
                                                                      computedXThermalImpulse, computedYThermalImpulse, computedZThermalImpulse, materialParameters);
        double computedTemperature = HPREquationOfState::computeTemperature(computedDensity, computedPressure, materialParameters);

        vector<double> totalEnergyDerivativeThermalImpulse = HPRDerivatives::computeTotalEnergyDerivativeThermalImpulse(computedXThermalImpulse, computedYThermalImpulse, computedZThermalImpulse,
                                                                                                                        materialParameters);

        double theta2Reciprocal = HPRSourceTerms::computeTheta2Reciprocal(computedDensity, computedTemperature, materialParameters);

        for (int i = 0; i < 3; i++)
        {
            sourceTermVector[14 + i] = -computedDensity * theta2Reciprocal * totalEnergyDerivativeThermalImpulse[i];
        }
    }

    return sourceTermVector;
}

vector<double> HPRStateVector::computeSourceTermVector(HPRMaterialParameters materialParameters)
{
    return computeSourceTermVector(computeConservedVariableVector(materialParameters), materialParameters);
}

double HPRStateVector::computeTotalEnergy(HPRMaterialParameters materialParameters)
{
    return HPREquationOfState::computeTotalEnergy(density, pressure, xVelocity, yVelocity, zVelocity, distortionTensor, xThermalImpulse, yThermalImpulse, zThermalImpulse, materialParameters);
}

double HPRStateVector::computeTemperature(HPRMaterialParameters materialParameters)
{
    return HPREquationOfState::computeTemperature(density, pressure, materialParameters);
}

vector<double> HPRStateVector::computeHeatFluxVector(HPRMaterialParameters materialParameters)
{
    return HPREquationOfState::computeHeatFluxVector(computeTemperature(materialParameters), xThermalImpulse, yThermalImpulse, zThermalImpulse, materialParameters);
}

double HPRStateVector::computeTotalEnergyDerivativeDensity(HPRMaterialParameters materialParameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativeDensity(density, pressure, materialParameters);
}

double HPRStateVector::computeTotalEnergyDerivativePressure(HPRMaterialParameters materialParameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativePressure(density, materialParameters);
}

vector<vector<double> > HPRStateVector::computeTotalEnergyDerivativeDistortionTensor(HPRMaterialParameters materialParameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativeDistortionTensor(distortionTensor, materialParameters);
}

vector<double> HPRStateVector::computeTotalEnergyDerivativeThermalImpulse(HPRMaterialParameters materialParameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativeThermalImpulse(xThermalImpulse, yThermalImpulse, zThermalImpulse, materialParameters);
}

double HPRStateVector::computeTemperatureDerivativeDensity(HPRMaterialParameters materialParameters)
{
    return HPRDerivatives::computeTemperatureDerivativeDensity(density, pressure, materialParameters);
}

double HPRStateVector::computeTemperatureDerivativePressure(HPRMaterialParameters materialParameters)
{
    return HPRDerivatives::computeTemperatureDerivativePressure(density, materialParameters);
}

double HPRStateVector::computeTheta1Reciprocal(HPRMaterialParameters materialParameters)
{
    return HPRSourceTerms::computeTheta1Reciprocal(density, distortionTensor, materialParameters);
}

double HPRStateVector::computeTheta2Reciprocal(HPRMaterialParameters materialParameters)
{
    return HPRSourceTerms::computeTheta2Reciprocal(density, computeTemperature(materialParameters), materialParameters);
}

vector<vector<double> > HPRStateVector::computeShearStressTensor(HPRMaterialParameters materialParameters)
{
    return HPREquationOfState::computeShearStressTensor(density, distortionTensor, materialParameters);
}

vector<vector<double> > HPRStateVector::computeShearStressTensorDerivativeDensity(HPRMaterialParameters materialParameters)
{
    return HPREquationOfState::computeShearStressTensorDerivativeDensity(distortionTensor, materialParameters);
}

vector<vector<vector<vector<double> > > > HPRStateVector::computeShearStressTensorDerivativeDistortionTensor(HPRMaterialParameters materialParameters)
{
    return HPREquationOfState::computeShearStressTensorDerivativeDistortionTensor(density, distortionTensor, materialParameters);
}

void HPRStateVector::setDensity(double newDensity)
{
    density = newDensity;
}

void HPRStateVector::setPressure(double newPressure)
{
    pressure = newPressure;
}

void HPRStateVector::setDistortionTensor(vector<vector<double> > newDistortionTensor)
{
    distortionTensor = newDistortionTensor;
}

void HPRStateVector::setXVelocity(double newXVelocity)
{
    xVelocity = newXVelocity;
}

void HPRStateVector::setYVelocity(double newYVelocity)
{
    yVelocity = newYVelocity;
}

void HPRStateVector::setZVelocity(double newZVelocity)
{
    zVelocity = newZVelocity;
}

void HPRStateVector::setXThermalImpulse(double newXThermalImpulse)
{
    xThermalImpulse = newXThermalImpulse;
}

void HPRStateVector::setYThermalImpulse(double newYThermalImpulse)
{
    yThermalImpulse = newYThermalImpulse;
}

void HPRStateVector::setZThermalImpulse(double newZThermalImpulse)
{
    zThermalImpulse = newZThermalImpulse;
}

double HPRStateVector::getDensity()
{
    return density;
}

double HPRStateVector::getPressure()
{
    return pressure;
}

vector<vector<double> > HPRStateVector::getDistortionTensor()
{
    return distortionTensor;
}

double HPRStateVector::getXVelocity()
{
    return xVelocity;
}

double HPRStateVector::getYVelocity()
{
    return yVelocity;
}

double HPRStateVector::getZVelocity()
{
    return zVelocity;
}

double HPRStateVector::getXThermalImpulse()
{
    return xThermalImpulse;
}

double HPRStateVector::getYThermalImpulse()
{
    return yThermalImpulse;
}

double HPRStateVector::getZThermalImpulse()
{
    return zThermalImpulse;
}

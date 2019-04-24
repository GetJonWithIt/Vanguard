#include "elasticstatevector.h"

ElasticStateVector::ElasticStateVector()
{
    density = 1.0;
    xVelocity = 0.0;
    yVelocity = 0.0;
    zVelocity = 0.0;

    distortionTensor = MatrixAlgebra::computeIdentityMatrix(3);
    entropy = 0.0;
}

ElasticStateVector::ElasticStateVector(double newXVelocity, double newYVelocity, double newZVelocity, vector<vector<double> > newDistortionTensor, double newEntropy,
                                       HyperelasticMaterialParameters materialParameters)
{
    double referenceMassDensity = materialParameters.getReferenceMassDensity();
    density = referenceMassDensity * MatrixAlgebra::computeDeterminant(newDistortionTensor);

    xVelocity = newXVelocity;
    yVelocity = newYVelocity;
    zVelocity = newZVelocity;
    distortionTensor = newDistortionTensor;
    entropy = newEntropy;
}

void ElasticStateVector::setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector)
{
    density = newPrimitiveVariableVector[0];
    xVelocity = newPrimitiveVariableVector[1];
    yVelocity = newPrimitiveVariableVector[2];
    zVelocity = newPrimitiveVariableVector[3];

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            distortionTensor[i][j] = newPrimitiveVariableVector[4 + (i * 3) + j];
        }
    }

    entropy = newPrimitiveVariableVector[13];
}

void ElasticStateVector::setConservedVariableVector(vector<double> newConservedVariableVector, HyperelasticMaterialParameters materialParameters)
{
    density = newConservedVariableVector[0];
    xVelocity = newConservedVariableVector[1] / density;
    yVelocity = newConservedVariableVector[2] / density;
    zVelocity = newConservedVariableVector[3] / density;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            distortionTensor[i][j] = newConservedVariableVector[4 + (i * 3) + j];
        }
    }

    double computedTotalEnergy = newConservedVariableVector[13] / density;
    entropy = ElasticEquationOfState::computeEntropy(computedTotalEnergy, distortionTensor, xVelocity, yVelocity, zVelocity, materialParameters);
}

vector<double> ElasticStateVector::computePrimitiveVariableVector()
{
    vector<double> primitiveVariableVector(14);

    primitiveVariableVector[0] = density;
    primitiveVariableVector[1] = xVelocity;
    primitiveVariableVector[2] = yVelocity;
    primitiveVariableVector[3] = zVelocity;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            primitiveVariableVector[4 + (i * 3) + j] = distortionTensor[i][j];
        }
    }

    primitiveVariableVector[13] = entropy;

    return primitiveVariableVector;
}

vector<double> ElasticStateVector::computeConservedVariableVector(HyperelasticMaterialParameters materialParameters)
{
    vector<double> conservedVariableVector(14);

    conservedVariableVector[0] = density;
    conservedVariableVector[1] = density * xVelocity;
    conservedVariableVector[2] = density * yVelocity;
    conservedVariableVector[3] = density * zVelocity;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            conservedVariableVector[4 + (i * 3) + j] = distortionTensor[i][j];
        }
    }

    double computedTotalEnergy = ElasticEquationOfState::computeTotalEnergy(distortionTensor, entropy, xVelocity, yVelocity, zVelocity, materialParameters);
    conservedVariableVector[13] = density * computedTotalEnergy;

    return conservedVariableVector;
}

vector<double> ElasticStateVector::computeXFluxVector(vector<double> conservedVariableVector, HyperelasticMaterialParameters materialParameters)
{
    vector<double> fluxVector(14);

    double computedDensity = conservedVariableVector[0];
    double computedXVelocity = conservedVariableVector[1] / computedDensity;
    double computedYVelocity = conservedVariableVector[2] / computedDensity;
    double computedZVelocity = conservedVariableVector[3] / computedDensity;

    vector<vector<double> > computedDistortionTensor(3, vector<double>(3));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            computedDistortionTensor[i][j] = conservedVariableVector[4 + (i * 3) + j];
        }
    }

    double computedTotalEnergy = conservedVariableVector[13] / computedDensity;
    double computedEntropy = ElasticEquationOfState::computeEntropy(computedTotalEnergy, computedDistortionTensor, computedXVelocity, computedYVelocity, computedZVelocity, materialParameters);

    vector<vector<double> > computedTotalStressTensor = ElasticEquationOfState::computeTotalStressTensor(computedDensity, computedDistortionTensor, computedEntropy, materialParameters);

    for (int i = 0; i < 14; i++)
    {
        fluxVector[i] = 0.0;
    }

    fluxVector[0] = computedDensity * computedXVelocity;
    fluxVector[1] = (computedDensity * (computedXVelocity * computedXVelocity)) - computedTotalStressTensor[0][0];
    fluxVector[2] = (computedDensity * (computedXVelocity * computedYVelocity)) - computedTotalStressTensor[0][1];
    fluxVector[3] = (computedDensity * (computedXVelocity * computedZVelocity)) - computedTotalStressTensor[0][2];

    vector<double> computedVelocityVector(3);
    computedVelocityVector[0] = computedXVelocity;
    computedVelocityVector[1] = computedYVelocity;
    computedVelocityVector[2] = computedZVelocity;
    vector<double> distortionTensorVelocityVectorProduct = MatrixAlgebra::multiplyMatrixByVector(computedDistortionTensor, computedVelocityVector);

    for (int i = 0; i < 3; i++)
    {
        fluxVector[4 + (i * 3)] = distortionTensorVelocityVectorProduct[i];
    }

    fluxVector[13] = (computedDensity * (computedXVelocity * computedTotalEnergy)) - VectorAlgebra::computeDotProduct(computedTotalStressTensor[0], computedVelocityVector);

    return fluxVector;
}

vector<double> ElasticStateVector::computeXFluxVector(HyperelasticMaterialParameters materialParameters)
{
    return computeXFluxVector(computeConservedVariableVector(materialParameters), materialParameters);
}

vector<double> ElasticStateVector::computeYFluxVector(vector<double> conservedVariableVector, HyperelasticMaterialParameters materialParameters)
{
    vector<double> fluxVector(14);

    double computedDensity = conservedVariableVector[0];
    double computedXVelocity = conservedVariableVector[1] / computedDensity;
    double computedYVelocity = conservedVariableVector[2] / computedDensity;
    double computedZVelocity = conservedVariableVector[3] / computedDensity;

    vector<vector<double> > computedDistortionTensor(3, vector<double>(3));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            computedDistortionTensor[i][j] = conservedVariableVector[4 + (i * 3) + j];
        }
    }

    double computedTotalEnergy = conservedVariableVector[13] / computedDensity;

    double computedEntropy = ElasticEquationOfState::computeEntropy(computedTotalEnergy, computedDistortionTensor, computedXVelocity, computedYVelocity, computedZVelocity, materialParameters);
    vector<vector<double> > computedTotalStressTensor = ElasticEquationOfState::computeTotalStressTensor(computedDensity, computedDistortionTensor, computedEntropy, materialParameters);

    for (int i = 0; i < 14; i++)
    {
        fluxVector[i] = 0.0;
    }

    fluxVector[0] = computedDensity * computedYVelocity;
    fluxVector[1] = (computedDensity * (computedYVelocity * computedXVelocity)) - computedTotalStressTensor[1][0];
    fluxVector[2] = (computedDensity * (computedYVelocity * computedYVelocity)) - computedTotalStressTensor[1][1];
    fluxVector[3] = (computedDensity * (computedYVelocity * computedZVelocity)) - computedTotalStressTensor[1][2];

    vector<double> computedVelocityVector(3);
    computedVelocityVector[0] = computedXVelocity;
    computedVelocityVector[1] = computedYVelocity;
    computedVelocityVector[2] = computedZVelocity;
    vector<double> distortionTensorVelocityVectorProduct = MatrixAlgebra::multiplyMatrixByVector(computedDistortionTensor, computedVelocityVector);

    for (int i = 0; i < 3; i++)
    {
        fluxVector[4 + (i * 3) + 1] = distortionTensorVelocityVectorProduct[i];
    }

    fluxVector[13] = (computedDensity * (computedYVelocity * computedTotalEnergy)) - VectorAlgebra::computeDotProduct(computedTotalStressTensor[1], computedVelocityVector);

    return fluxVector;
}

vector<double> ElasticStateVector::computeYFluxVector(HyperelasticMaterialParameters materialParameters)
{
    return computeYFluxVector(computeConservedVariableVector(materialParameters), materialParameters);
}

double ElasticStateVector::computeSoundSpeed(HyperelasticMaterialParameters materialParameters, int direction)
{
    return ElasticAcousticTensor::computeMaximumWaveSpeed(distortionTensor, entropy, materialParameters, direction);
}

void ElasticStateVector::setDensity(double newDensity)
{
    density = newDensity;
}

void ElasticStateVector::setXVelocity(double newXVelocity)
{
    xVelocity = newXVelocity;
}

void ElasticStateVector::setYVelocity(double newYVelocity)
{
    yVelocity = newYVelocity;
}

void ElasticStateVector::setZVelocity(double newZVelocity)
{
    zVelocity = newZVelocity;
}

void ElasticStateVector::setDistortionTensor(vector<vector<double> > newDistortionTensor)
{
    distortionTensor = newDistortionTensor;
}

void ElasticStateVector::setEntropy(double newEntropy)
{
    entropy = newEntropy;
}

double ElasticStateVector::getDensity()
{
    return density;
}

double ElasticStateVector::getXVelocity()
{
    return xVelocity;
}

double ElasticStateVector::getYVelocity()
{
    return yVelocity;
}

double ElasticStateVector::getZVelocity()
{
    return zVelocity;
}

vector<vector<double> > ElasticStateVector::getDistortionTensor()
{
    return distortionTensor;
}

double ElasticStateVector::getEntropy()
{
    return entropy;
}

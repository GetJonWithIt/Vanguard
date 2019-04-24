#include "elasticreducedstatevector.h"

ElasticReducedStateVector::ElasticReducedStateVector()
{
    material1VolumeFraction = 0.999;
    interfaceXVelocity = 0.0;
    interfaceYVelocity = 0.0;
    interfaceZVelocity = 0.0;

    interfaceDistortionTensor = MatrixAlgebra::computeIdentityMatrix(3);
    interfaceEntropy = 0.0;

    material1Density = 1.0;
    material2Density = 1.0;
}

ElasticReducedStateVector::ElasticReducedStateVector(double newMaterial1VolumeFraction, double newInterfaceXVelocity, double newInterfaceYVelocity, double newInterfaceZVelocity,
                                                     vector<vector<double> > newMaterial1DistortionTensor, double newMaterial1Entropy, vector<vector<double> > newMaterial2DistortionTensor,
                                                     double newMaterial2Entropy, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    double material1ReferenceMassDensity = material1Parameters.getReferenceMassDensity();
    material1Density = material1ReferenceMassDensity * MatrixAlgebra::computeDeterminant(newMaterial1DistortionTensor);

    double material2ReferenceMassDensity = material2Parameters.getReferenceMassDensity();
    material2Density = material2ReferenceMassDensity * MatrixAlgebra::computeDeterminant(newMaterial2DistortionTensor);

    material1VolumeFraction = newMaterial1VolumeFraction;
    interfaceXVelocity = newInterfaceXVelocity;
    interfaceYVelocity = newInterfaceYVelocity;
    interfaceZVelocity = newInterfaceZVelocity;

    double material2VolumeFraction = 1.0 - material1VolumeFraction;
    interfaceDistortionTensor = MatrixAlgebra::addMatrices(MatrixAlgebra::multiplyMatrix(pow(material1VolumeFraction, (1.0 / 3.0)), newMaterial1DistortionTensor),
                                                           MatrixAlgebra::multiplyMatrix(pow(material2VolumeFraction, (1.0 / 3.0)), newMaterial2DistortionTensor));

    interfaceEntropy = (material1VolumeFraction * newMaterial1Entropy) + (material2VolumeFraction * newMaterial2Entropy);
}

void ElasticReducedStateVector::setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector)
{
    material1VolumeFraction = newPrimitiveVariableVector[0];
    interfaceXVelocity = newPrimitiveVariableVector[1];
    interfaceYVelocity = newPrimitiveVariableVector[2];
    interfaceZVelocity = newPrimitiveVariableVector[3];

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            interfaceDistortionTensor[i][j] = newPrimitiveVariableVector[4 + (i * 3) + j];
        }
    }
    interfaceEntropy = newPrimitiveVariableVector[13];

    material1Density = newPrimitiveVariableVector[14];
    material2Density = newPrimitiveVariableVector[15];
}

void ElasticReducedStateVector::setConservedVariableVector(vector<double> newConservedVariableVector, HyperelasticMaterialParameters material1Parameters,
                                                           HyperelasticMaterialParameters material2Parameters)
{
    double totalDensity = newConservedVariableVector[0];

    if ((newConservedVariableVector[1] / totalDensity) < 0.001)
    {
        material1VolumeFraction = 0.001;
    }
    else if ((newConservedVariableVector[1] / totalDensity) > 0.999)
    {
        material1VolumeFraction = 0.999;
    }
    else
    {
        material1VolumeFraction = newConservedVariableVector[1] / totalDensity;
    }
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    interfaceXVelocity = newConservedVariableVector[2] / totalDensity;
    interfaceYVelocity = newConservedVariableVector[3] / totalDensity;
    interfaceZVelocity = newConservedVariableVector[4] / totalDensity;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            interfaceDistortionTensor[i][j] = newConservedVariableVector[5 + (i * 3) + j];
        }
    }

    double computedTotalEnergy = newConservedVariableVector[14] / totalDensity;
    double computedMaterial1Entropy = ElasticEquationOfState::computeEntropy(computedTotalEnergy, interfaceDistortionTensor, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                                             material1Parameters);
    double computedMaterial2Entropy = ElasticEquationOfState::computeEntropy(computedTotalEnergy, interfaceDistortionTensor, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                                             material2Parameters);
    interfaceEntropy = (material1VolumeFraction * computedMaterial1Entropy) + (material2VolumeFraction * computedMaterial2Entropy);

    material1Density = newConservedVariableVector[15] / material1VolumeFraction;
    material2Density = newConservedVariableVector[16] / material2VolumeFraction;
}

vector<double> ElasticReducedStateVector::computePrimitiveVariableVector()
{
    vector<double> primitiveVariableVector(16);

    primitiveVariableVector[0] = material1VolumeFraction;
    primitiveVariableVector[1] = interfaceXVelocity;
    primitiveVariableVector[2] = interfaceYVelocity;
    primitiveVariableVector[3] = interfaceZVelocity;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            primitiveVariableVector[4 + (i * 3) + j] = interfaceDistortionTensor[i][j];
        }
    }
    primitiveVariableVector[13] = interfaceEntropy;

    primitiveVariableVector[14] = material1Density;
    primitiveVariableVector[15] = material2Density;

    return primitiveVariableVector;
}

vector<double> ElasticReducedStateVector::computeConservedVariableVector(HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    vector<double> conservedVariableVector(17);

    if (material1VolumeFraction < 0.001)
    {
        material1VolumeFraction = 0.001;
    }
    else if (material1VolumeFraction > 0.999)
    {
        material1VolumeFraction = 0.999;
    }
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    double totalDensity = (material1VolumeFraction * material1Density) + (material2VolumeFraction * material2Density);

    conservedVariableVector[0] = totalDensity;
    conservedVariableVector[1] = totalDensity * material1VolumeFraction;
    conservedVariableVector[2] = totalDensity * interfaceXVelocity;
    conservedVariableVector[3] = totalDensity * interfaceYVelocity;
    conservedVariableVector[4] = totalDensity * interfaceZVelocity;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            conservedVariableVector[5 + (i * 3) + j] = interfaceDistortionTensor[i][j];
        }
    }

    double computedMaterial1TotalEnergy = ElasticEquationOfState::computeTotalEnergy(interfaceDistortionTensor, interfaceEntropy, interfaceXVelocity, interfaceYVelocity,
                                                                                     interfaceZVelocity, material1Parameters);
    double computedMaterial2TotalEnergy = ElasticEquationOfState::computeTotalEnergy(interfaceDistortionTensor, interfaceEntropy, interfaceXVelocity, interfaceYVelocity,
                                                                                     interfaceZVelocity, material2Parameters);
    double computedTotalEnergy = (material1VolumeFraction * computedMaterial1TotalEnergy) + (material2VolumeFraction * computedMaterial2TotalEnergy);
    conservedVariableVector[14] = totalDensity * computedTotalEnergy;

    conservedVariableVector[15] = material1VolumeFraction * material1Density;
    conservedVariableVector[16] = material2VolumeFraction * material2Density;

    return conservedVariableVector;
}

vector<double> ElasticReducedStateVector::computeXFluxVector(vector<double> conservedVariableVector, HyperelasticMaterialParameters material1Parameters,
                                                             HyperelasticMaterialParameters material2Parameters)
{
    vector<double> fluxVector(17);

    double computedTotalDensity = conservedVariableVector[0];

    double computedMaterial1VolumeFraction;
    if ((conservedVariableVector[1] / computedTotalDensity) < 0.001)
    {
        computedMaterial1VolumeFraction = 0.001;
    }
    else if ((conservedVariableVector[1] / computedTotalDensity) > 0.999)
    {
        computedMaterial1VolumeFraction = 0.999;
    }
    else
    {
        computedMaterial1VolumeFraction = conservedVariableVector[1] / computedTotalDensity;
    }
    double computedMaterial2VolumeFraction = 1.0 - computedMaterial1VolumeFraction;

    double computedInterfaceXVelocity = conservedVariableVector[2] / computedTotalDensity;
    double computedInterfaceYVelocity = conservedVariableVector[3] / computedTotalDensity;
    double computedInterfaceZVelocity = conservedVariableVector[4] / computedTotalDensity;

    vector<vector<double> > computedInterfaceDistortionTensor(3, vector<double>(3));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            computedInterfaceDistortionTensor[i][j] = conservedVariableVector[5 + (i * 3) + j];
        }
    }

    double computedTotalEnergy = conservedVariableVector[14] / computedTotalDensity;
    double computedMaterial1Entropy = ElasticEquationOfState::computeEntropy(computedTotalEnergy, computedInterfaceDistortionTensor, computedInterfaceXVelocity,
                                                                             computedInterfaceYVelocity, computedInterfaceZVelocity, material1Parameters);
    double computedMaterial2Entropy = ElasticEquationOfState::computeEntropy(computedTotalEnergy, computedInterfaceDistortionTensor, computedInterfaceXVelocity,
                                                                             computedInterfaceZVelocity, computedInterfaceZVelocity, material2Parameters);
    double computedInterfaceEntropy = (computedMaterial1VolumeFraction * computedMaterial1Entropy) + (computedMaterial2VolumeFraction * computedMaterial2Entropy);

    vector<vector<double> > computedMaterial1TotalStressTensor = ElasticEquationOfState::computeTotalStressTensor(computedTotalDensity, computedInterfaceDistortionTensor,
                                                                                                                  computedInterfaceEntropy, material1Parameters);
    vector<vector<double> > computedMaterial2TotalStressTensor = ElasticEquationOfState::computeTotalStressTensor(computedTotalDensity, computedInterfaceDistortionTensor,
                                                                                                                  computedInterfaceEntropy, material2Parameters);
    vector<vector<double> > computedInterfaceTotalStressTensor = MatrixAlgebra::addMatrices(MatrixAlgebra::multiplyMatrix(computedMaterial1VolumeFraction, computedMaterial1TotalStressTensor),
                                                                                            MatrixAlgebra::multiplyMatrix(computedMaterial2VolumeFraction, computedMaterial2TotalStressTensor));

    double computedMaterial1Density = conservedVariableVector[15] / computedMaterial1VolumeFraction;
    double computedMaterial2Density = conservedVariableVector[16] / computedMaterial2VolumeFraction;

    vector<double> computedVelocityVector(3);
    computedVelocityVector[0] = computedInterfaceXVelocity;
    computedVelocityVector[1] = computedInterfaceYVelocity;
    computedVelocityVector[2] = computedInterfaceZVelocity;

    for (int i = 0; i < 17; i++)
    {
        fluxVector[i] = 0.0;
    }

    fluxVector[0] = computedTotalDensity * computedInterfaceXVelocity;
    fluxVector[1] = computedTotalDensity * (computedInterfaceXVelocity * computedMaterial1VolumeFraction);
    fluxVector[2] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceXVelocity)) - computedInterfaceTotalStressTensor[0][0];
    fluxVector[3] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceYVelocity)) - computedInterfaceTotalStressTensor[0][1];
    fluxVector[4] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceZVelocity)) - computedInterfaceTotalStressTensor[0][2];

    vector<double> interfaceDistortionTensorVelocityVectorProduct = MatrixAlgebra::multiplyMatrixByVector(computedInterfaceDistortionTensor, computedVelocityVector);
    for (int i = 0; i < 3; i++)
    {
        fluxVector[5 + (i * 3)] = interfaceDistortionTensorVelocityVectorProduct[i];
    }

    fluxVector[14] = (computedTotalDensity * (computedInterfaceXVelocity * computedTotalEnergy)) - VectorAlgebra::computeDotProduct(computedInterfaceTotalStressTensor[0], computedVelocityVector);

    fluxVector[15] = computedMaterial1VolumeFraction * (computedMaterial1Density * computedInterfaceXVelocity);
    fluxVector[16] = computedMaterial2VolumeFraction * (computedMaterial2Density * computedInterfaceXVelocity);

    return fluxVector;
}

vector<double> ElasticReducedStateVector::computeXFluxVector(HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    return computeXFluxVector(computeConservedVariableVector(material1Parameters, material2Parameters), material1Parameters, material2Parameters);
}

vector<double> ElasticReducedStateVector::computeYFluxVector(vector<double> conservedVariableVector, HyperelasticMaterialParameters material1Parameters,
                                                             HyperelasticMaterialParameters material2Parameters)
{
    vector<double> fluxVector(17);

    double computedTotalDensity = conservedVariableVector[0];

    double computedMaterial1VolumeFraction;
    if ((conservedVariableVector[1] / computedTotalDensity) < 0.001)
    {
        computedMaterial1VolumeFraction = 0.001;
    }
    else if ((conservedVariableVector[1] / computedTotalDensity) > 0.999)
    {
        computedMaterial1VolumeFraction = 0.999;
    }
    else
    {
        computedMaterial1VolumeFraction = conservedVariableVector[1] / computedTotalDensity;
    }
    double computedMaterial2VolumeFraction = 1.0 - computedMaterial1VolumeFraction;

    double computedInterfaceXVelocity = conservedVariableVector[2] / computedTotalDensity;
    double computedInterfaceYVelocity = conservedVariableVector[3] / computedTotalDensity;
    double computedInterfaceZVelocity = conservedVariableVector[4] / computedTotalDensity;

    vector<vector<double> > computedInterfaceDistortionTensor(3, vector<double>(3));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            computedInterfaceDistortionTensor[i][j] = conservedVariableVector[5 + (i * 3) + j];
        }
    }

    double computedTotalEnergy = conservedVariableVector[14] / computedTotalDensity;
    double computedMaterial1Entropy = ElasticEquationOfState::computeEntropy(computedTotalEnergy, computedInterfaceDistortionTensor, computedInterfaceXVelocity,
                                                                             computedInterfaceYVelocity, computedInterfaceZVelocity, material1Parameters);
    double computedMaterial2Entropy = ElasticEquationOfState::computeEntropy(computedTotalEnergy, computedInterfaceDistortionTensor, computedInterfaceXVelocity,
                                                                             computedInterfaceYVelocity, computedInterfaceZVelocity, material2Parameters);
    double computedInterfaceEntropy = (computedMaterial1VolumeFraction * computedMaterial1Entropy) + (computedMaterial2VolumeFraction * computedMaterial2Entropy);

    vector<vector<double> > computedMaterial1TotalStressTensor = ElasticEquationOfState::computeTotalStressTensor(computedTotalDensity, computedInterfaceDistortionTensor,
                                                                                                                  computedInterfaceEntropy, material1Parameters);
    vector<vector<double> > computedMaterial2TotalStressTensor = ElasticEquationOfState::computeTotalStressTensor(computedTotalDensity, computedInterfaceDistortionTensor,
                                                                                                                  computedInterfaceEntropy, material2Parameters);
    vector<vector<double> > computedInterfaceTotalStressTensor = MatrixAlgebra::addMatrices(MatrixAlgebra::multiplyMatrix(computedMaterial1VolumeFraction, computedMaterial1TotalStressTensor),
                                                                                            MatrixAlgebra::multiplyMatrix(computedMaterial2VolumeFraction, computedMaterial2TotalStressTensor));

    double computedMaterial1Density = conservedVariableVector[15] / computedMaterial1VolumeFraction;
    double computedMaterial2Density = conservedVariableVector[16] / computedMaterial2VolumeFraction;

    vector<double> computedVelocityVector(3);
    computedVelocityVector[0] = computedInterfaceXVelocity;
    computedVelocityVector[1] = computedInterfaceYVelocity;
    computedVelocityVector[2] = computedInterfaceZVelocity;

    for (int i = 0; i < 17; i++)
    {
        fluxVector[i] = 0.0;
    }

    fluxVector[0] = computedTotalDensity * computedInterfaceYVelocity;
    fluxVector[1] = computedTotalDensity * (computedInterfaceYVelocity * computedMaterial1VolumeFraction);
    fluxVector[2] = (computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceXVelocity)) - computedInterfaceTotalStressTensor[1][0];
    fluxVector[3] = (computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceYVelocity)) - computedInterfaceTotalStressTensor[1][1];
    fluxVector[4] = (computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceZVelocity)) - computedInterfaceTotalStressTensor[1][2];

    vector<double> interfaceDistortionTensorVelocityVectorProduct = MatrixAlgebra::multiplyMatrixByVector(computedInterfaceDistortionTensor, computedVelocityVector);
    for (int i = 0; i < 3; i++)
    {
        fluxVector[5 + (i * 3) + 1] = interfaceDistortionTensorVelocityVectorProduct[i];
    }

    fluxVector[14] = (computedTotalDensity * (computedInterfaceYVelocity * computedTotalEnergy)) - VectorAlgebra::computeDotProduct(computedInterfaceTotalStressTensor[1], computedVelocityVector);

    fluxVector[15] = computedMaterial1VolumeFraction * (computedMaterial1Density * computedInterfaceYVelocity);
    fluxVector[16] = computedMaterial2VolumeFraction * (computedMaterial2Density * computedInterfaceYVelocity);

    return fluxVector;
}

vector<double> ElasticReducedStateVector::computeYFluxVector(HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    return computeYFluxVector(computeConservedVariableVector(material1Parameters, material2Parameters), material1Parameters, material2Parameters);
}

double ElasticReducedStateVector::computeMaterial1SoundSpeed(HyperelasticMaterialParameters material1Parameters, int direction)
{
    return ElasticAcousticTensor::computeMaximumWaveSpeed(interfaceDistortionTensor, interfaceEntropy, material1Parameters, direction);
}

double ElasticReducedStateVector::computeMaterial2SoundSpeed(HyperelasticMaterialParameters material2Parameters, int direction)
{
    return ElasticAcousticTensor::computeMaximumWaveSpeed(interfaceDistortionTensor, interfaceEntropy, material2Parameters, direction);
}

double ElasticReducedStateVector::computeTotalDensity()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1Density) + (material2VolumeFraction * material2Density);
}

void ElasticReducedStateVector::relaxTotalDensity()
{
    double totalDensity = computeTotalDensity();

    material1Density = totalDensity;
    material2Density = totalDensity;
}

void ElasticReducedStateVector::setMaterial1VolumeFraction(double newMaterial1VolumeFraction)
{
    material1VolumeFraction = newMaterial1VolumeFraction;
}

void ElasticReducedStateVector::setInterfaceXVelocity(double newInterfaceXVelocity)
{
    interfaceXVelocity = newInterfaceXVelocity;
}

void ElasticReducedStateVector::setInterfaceYVelocity(double newInterfaceYVelocity)
{
    interfaceYVelocity = newInterfaceYVelocity;
}

void ElasticReducedStateVector::setInterfaceZVelocity(double newInterfaceZVelocity)
{
    interfaceZVelocity = newInterfaceZVelocity;
}

void ElasticReducedStateVector::setInterfaceDistortionTensor(vector<vector<double> > newInterfaceDistortionTensor)
{
    interfaceDistortionTensor = newInterfaceDistortionTensor;
}

void ElasticReducedStateVector::setInterfaceEntropy(double newInterfaceEntropy)
{
    interfaceEntropy = newInterfaceEntropy;
}

void ElasticReducedStateVector::setMaterial1Density(double newMaterial1Density)
{
    material1Density = newMaterial1Density;
}

void ElasticReducedStateVector::setMaterial2Density(double newMaterial2Density)
{
    material2Density = newMaterial2Density;
}

double ElasticReducedStateVector::getMaterial1VolumeFraction()
{
    return material1VolumeFraction;
}

double ElasticReducedStateVector::getInterfaceXVelocity()
{
    return interfaceXVelocity;
}

double ElasticReducedStateVector::getInterfaceYVelocity()
{
    return interfaceYVelocity;
}

double ElasticReducedStateVector::getInterfaceZVelocity()
{
    return interfaceZVelocity;
}

vector<vector<double> > ElasticReducedStateVector::getInterfaceDistortionTensor()
{
    return interfaceDistortionTensor;
}

double ElasticReducedStateVector::getInterfaceEntropy()
{
    return interfaceEntropy;
}

double ElasticReducedStateVector::getMaterial1Density()
{
    return material1Density;
}

double ElasticReducedStateVector::getMaterial2Density()
{
    return material2Density;
}

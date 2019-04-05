#include "elasticmultiphysicsstatevector.h"

ElasticMultiphysicsStateVector::ElasticMultiphysicsStateVector()
{
    material1VolumeFraction = 0.999;
    interfaceXVelocity = 0.0;
    interfaceYVelocity = 0.0;
    interfaceZVelocity = 0.0;

    material1Density = 1.0;
    material1DistortionTensor = MatrixAlgebra::computeIdentityMatrix(3);
    material1Entropy = 0.0;

    material2Density = 1.0;
    material2DistortionTensor = MatrixAlgebra::computeIdentityMatrix(3);
    material2Entropy = 0.0;
}

ElasticMultiphysicsStateVector::ElasticMultiphysicsStateVector(double newMaterial1VolumeFraction, double newInterfaceXVelocity, double newInterfaceYVelocity, double newInterfaceZVelocity,
                                                               vector<vector<double> > newMaterial1DistortionTensor, double newMaterial1Entropy,
                                                               vector<vector<double> > newMaterial2DistortionTensor, double newMaterial2Entropy,
                                                               HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    double material1ReferenceMassDensity = material1Parameters.getReferenceMassDensity();
    material1Density = material1ReferenceMassDensity * MatrixAlgebra::computeDeterminant(newMaterial1DistortionTensor);

    double material2ReferenceMassDensity = material2Parameters.getReferenceMassDensity();
    material2Density = material2ReferenceMassDensity * MatrixAlgebra::computeDeterminant(newMaterial2DistortionTensor);

    material1VolumeFraction = newMaterial1VolumeFraction;
    interfaceXVelocity = newInterfaceXVelocity;
    interfaceYVelocity = newInterfaceYVelocity;
    interfaceZVelocity = newInterfaceZVelocity;

    material1DistortionTensor = newMaterial1DistortionTensor;
    material1Entropy = newMaterial1Entropy;

    material2DistortionTensor = newMaterial2DistortionTensor;
    material2Entropy = newMaterial2Entropy;
}

void ElasticMultiphysicsStateVector::setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector)
{
    material1VolumeFraction = newPrimitiveVariableVector[0];
    interfaceXVelocity = newPrimitiveVariableVector[1];
    interfaceYVelocity = newPrimitiveVariableVector[2];
    interfaceZVelocity = newPrimitiveVariableVector[3];

    material1Density = newPrimitiveVariableVector[4];
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            material1DistortionTensor[i][j] = newPrimitiveVariableVector[5 + (i * 3) + j];
        }
    }

    material1Entropy = newPrimitiveVariableVector[14];

    material2Density = newPrimitiveVariableVector[15];
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            material2DistortionTensor[i][j] = newPrimitiveVariableVector[16 + (i * 3) + j];
        }
    }

    material2Entropy = newPrimitiveVariableVector[25];
}

void ElasticMultiphysicsStateVector::setConservedVariableVector(vector<double> newConservedVariableVector, HyperelasticMaterialParameters material1Parameters,
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

    material1Density = newConservedVariableVector[5] / material1VolumeFraction;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            material1DistortionTensor[i][j] = newConservedVariableVector[6 + (i * 3) + j] / pow(material1VolumeFraction, (1.0 / 3.0));
        }
    }

    double computedMaterial1TotalEnergy = newConservedVariableVector[15] / (material1VolumeFraction * material1Density);
    material1Entropy = ElasticEquationOfState::computeEntropy(computedMaterial1TotalEnergy, material1DistortionTensor, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                              material1Parameters);

    material2Density = newConservedVariableVector[16] / material2VolumeFraction;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            material2DistortionTensor[i][j] = newConservedVariableVector[17 + (i * 3) + j] / pow(material2VolumeFraction, (1.0 / 3.0));
        }
    }

    double computedMaterial2TotalEnergy = newConservedVariableVector[26] / (material2VolumeFraction * material2Density);
    material2Entropy = ElasticEquationOfState::computeEntropy(computedMaterial2TotalEnergy, material2DistortionTensor, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                              material2Parameters);
}

vector<double> ElasticMultiphysicsStateVector::computePrimitiveVariableVector()
{
    vector<double> primitiveVariableVector(26);

    primitiveVariableVector[0] = material1VolumeFraction;
    primitiveVariableVector[1] = interfaceXVelocity;
    primitiveVariableVector[2] = interfaceYVelocity;
    primitiveVariableVector[3] = interfaceZVelocity;

    primitiveVariableVector[4] = material1Density;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            primitiveVariableVector[5 + (i * 3) + j] = material1DistortionTensor[i][j];
        }
    }

    primitiveVariableVector[14] = material1Entropy;

    primitiveVariableVector[15] = material2Density;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            primitiveVariableVector[16 + (i * 3) + j] = material2DistortionTensor[i][j];
        }
    }

    primitiveVariableVector[25] = material2Entropy;

    return primitiveVariableVector;
}

vector<double> ElasticMultiphysicsStateVector::computeConservedVariableVector(HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    vector<double> conservedVariableVector(27);

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

    conservedVariableVector[5] = material1VolumeFraction * material1Density;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            conservedVariableVector[6 + (i * 3) + j] = pow(material1VolumeFraction, (1.0 / 3.0)) * material1DistortionTensor[i][j];
        }
    }

    double computedMaterial1TotalEnergy = ElasticEquationOfState::computeTotalEnergy(material1DistortionTensor, material1Entropy, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                                                     material1Parameters);
    conservedVariableVector[15] = material1VolumeFraction * (material1Density * computedMaterial1TotalEnergy);

    conservedVariableVector[16] = material2VolumeFraction * material2Density;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            conservedVariableVector[17 + (i * 3) + j] = pow(material2VolumeFraction, (1.0 / 3.0)) * material2DistortionTensor[i][j];
        }
    }

    double computedMaterial2TotalEnergy = ElasticEquationOfState::computeTotalEnergy(material2DistortionTensor, material2Entropy, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                                                     material2Parameters);
    conservedVariableVector[26] = material2VolumeFraction * (material2Density * computedMaterial2TotalEnergy);

    return conservedVariableVector;
}

vector<double> ElasticMultiphysicsStateVector::computeXFluxVector(vector<double> conservedVariableVector, HyperelasticMaterialParameters material1Parameters,
                                                                  HyperelasticMaterialParameters material2Parameters)
{
    vector<double> fluxVector(27);

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

    double computedMaterial1Density = conservedVariableVector[5] / computedMaterial1VolumeFraction;

    vector<vector<double> > computedMaterial1DistortionTensor(3, vector<double>(3));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            computedMaterial1DistortionTensor[i][j] = conservedVariableVector[6 + (i * 3) + j] / pow(computedMaterial1VolumeFraction, (1.0 / 3.0));
        }
    }

    double computedMaterial1TotalEnergy = conservedVariableVector[15] / (computedMaterial1VolumeFraction * computedMaterial1Density);
    double computedMaterial1Entropy = ElasticEquationOfState::computeEntropy(computedMaterial1TotalEnergy, computedMaterial1DistortionTensor, computedInterfaceXVelocity,
                                                                             computedInterfaceYVelocity, computedInterfaceZVelocity, material1Parameters);

    vector<vector<double> > computedMaterial1TotalStressTensor = ElasticEquationOfState::computeTotalStressTensor(computedMaterial1Density, computedMaterial1DistortionTensor,
                                                                                                                  computedMaterial1Entropy, material1Parameters);

    double computedMaterial2Density = conservedVariableVector[16] / computedMaterial2VolumeFraction;

    vector<vector<double> > computedMaterial2DistortionTensor(3, vector<double>(3));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            computedMaterial2DistortionTensor[i][j] = conservedVariableVector[17 + (i * 3) + j] / pow(computedMaterial2VolumeFraction, (1.0 / 3.0));
        }
    }

    double computedMaterial2TotalEnergy = conservedVariableVector[26] / (computedMaterial2VolumeFraction * computedMaterial2Density);
    double computedMaterial2Entropy = ElasticEquationOfState::computeEntropy(computedMaterial2TotalEnergy, computedMaterial2DistortionTensor, computedInterfaceXVelocity,
                                                                             computedInterfaceYVelocity, computedInterfaceZVelocity, material2Parameters);

    vector<vector<double> > computedMaterial2TotalStressTensor = ElasticEquationOfState::computeTotalStressTensor(computedMaterial2Density, computedMaterial2DistortionTensor,
                                                                                                                  computedMaterial2Entropy, material2Parameters);

    vector<double> computedVelocityVector(3);
    computedVelocityVector[0] = computedInterfaceXVelocity;
    computedVelocityVector[1] = computedInterfaceYVelocity;
    computedVelocityVector[2] = computedInterfaceZVelocity;

    for (int i = 0; i < 27; i++)
    {
        fluxVector[i] = 0.0;
    }

    fluxVector[0] = computedTotalDensity * computedInterfaceXVelocity;
    fluxVector[1] = computedTotalDensity * (computedInterfaceXVelocity * computedMaterial1VolumeFraction);
    fluxVector[2] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceXVelocity)) -
            ((computedMaterial1VolumeFraction * computedMaterial1TotalStressTensor[0][0]) + (computedMaterial2VolumeFraction * computedMaterial2TotalStressTensor[0][0]));
    fluxVector[3] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceYVelocity)) -
            ((computedMaterial1VolumeFraction * computedMaterial1TotalStressTensor[0][1]) + (computedMaterial2VolumeFraction * computedMaterial2TotalStressTensor[0][1]));
    fluxVector[4] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceZVelocity)) -
            ((computedMaterial1VolumeFraction * computedMaterial1TotalStressTensor[0][2]) + (computedMaterial2VolumeFraction * computedMaterial2TotalStressTensor[0][2]));

    fluxVector[5] = computedMaterial1VolumeFraction * (computedMaterial1Density * computedInterfaceXVelocity);

    vector<double> material1DistortionTensorVelocityVectorProduct = MatrixAlgebra::multiplyMatrixByVector(computedMaterial1DistortionTensor, computedVelocityVector);
    for (int i = 0; i < 3; i++)
    {
        fluxVector[6 + (i * 3)] = pow(computedMaterial1VolumeFraction, (1.0 / 3.0)) * material1DistortionTensorVelocityVectorProduct[i];
    }

    fluxVector[15] = computedMaterial1VolumeFraction * ((computedMaterial1Density * (computedInterfaceXVelocity * computedMaterial1TotalEnergy)) -
                                                        VectorAlgebra::computeDotProduct(computedMaterial1TotalStressTensor[0], computedVelocityVector));

    fluxVector[16] = computedMaterial2VolumeFraction * (computedMaterial2Density * computedInterfaceXVelocity);

    vector<double> material2DistortionTensorVelocityVectorProduct = MatrixAlgebra::multiplyMatrixByVector(computedMaterial2DistortionTensor, computedVelocityVector);
    for (int i = 0; i < 3; i++)
    {
        fluxVector[17 + (i * 3)] = pow(computedMaterial2VolumeFraction, (1.0 / 3.0)) * material2DistortionTensorVelocityVectorProduct[i];
    }

    fluxVector[26] = computedMaterial2VolumeFraction * ((computedMaterial2Density * (computedInterfaceXVelocity * computedMaterial2TotalEnergy)) -
                                                        VectorAlgebra::computeDotProduct(computedMaterial2TotalStressTensor[0], computedVelocityVector));

    return fluxVector;
}

vector<double> ElasticMultiphysicsStateVector::computeXFluxVector(HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    return computeXFluxVector(computeConservedVariableVector(material1Parameters, material2Parameters), material1Parameters, material2Parameters);
}

double ElasticMultiphysicsStateVector::computeMaterial1TotalEnergy(HyperelasticMaterialParameters material1Parameters)
{
    return ElasticEquationOfState::computeTotalEnergy(material1DistortionTensor, material1Entropy, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity, material1Parameters);
}

double ElasticMultiphysicsStateVector::computeMaterial1SoundSpeed(HyperelasticMaterialParameters material1Parameters, int direction)
{
    return ElasticAcousticTensor::computeMaximumWaveSpeed(material1DistortionTensor, material1Entropy, material1Parameters, direction);
}

double ElasticMultiphysicsStateVector::computeMaterial2TotalEnergy(HyperelasticMaterialParameters material2Parameters)
{
    return ElasticEquationOfState::computeTotalEnergy(material2DistortionTensor, material2Entropy, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity, material2Parameters);
}

double ElasticMultiphysicsStateVector::computeMaterial2SoundSpeed(HyperelasticMaterialParameters material2Parameters, int direction)
{
    return ElasticAcousticTensor::computeMaximumWaveSpeed(material2DistortionTensor, material2Entropy, material2Parameters, direction);
}

double ElasticMultiphysicsStateVector::computeTotalDensity()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1Density) + (material2VolumeFraction * material2Density);
}

vector<vector<double> > ElasticMultiphysicsStateVector::computeTotalDistortionTensor()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return MatrixAlgebra::addMatrices(MatrixAlgebra::multiplyMatrix(pow(material1VolumeFraction, (1.0 / 3.0)), material1DistortionTensor),
                                      MatrixAlgebra::multiplyMatrix(pow(material2VolumeFraction, (1.0 / 3.0)), material2DistortionTensor));
}

double ElasticMultiphysicsStateVector::computeTotalEntropy()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1Entropy) + (material2VolumeFraction * material2Entropy);
}

void ElasticMultiphysicsStateVector::relaxTotalDensity()
{
    double totalDensity = computeTotalDensity();

    material1Density = totalDensity;
    material2Density = totalDensity;
}

void ElasticMultiphysicsStateVector::relaxTotalDistortionTensor()
{
    vector<vector<double> > totalDistortionTensor = computeTotalDistortionTensor();

    material1DistortionTensor = totalDistortionTensor;
    material2DistortionTensor = totalDistortionTensor;
}

void ElasticMultiphysicsStateVector::relaxTotalEntropy(HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;
    double totalDensity = computeTotalDensity();
    double totalEnergy = ((material1VolumeFraction * material1Density * computeMaterial1TotalEnergy(material1Parameters)) + (material2VolumeFraction * material2Density *
                                                                                                                             computeMaterial2TotalEnergy(material2Parameters))) / totalDensity;

    material1Entropy = ElasticEquationOfState::computeEntropy(totalEnergy, material1DistortionTensor, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity, material1Parameters);
    material2Entropy = ElasticEquationOfState::computeEntropy(totalEnergy, material2DistortionTensor, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity, material2Parameters);
}

void ElasticMultiphysicsStateVector::setMaterial1VolumeFraction(double newMaterial1VolumeFraction)
{
    material1VolumeFraction = newMaterial1VolumeFraction;
}

void ElasticMultiphysicsStateVector::setInterfaceXVelocity(double newInterfaceXVelocity)
{
    interfaceXVelocity = newInterfaceXVelocity;
}

void ElasticMultiphysicsStateVector::setInterfaceYVelocity(double newInterfaceYVelocity)
{
    interfaceYVelocity = newInterfaceYVelocity;
}

void ElasticMultiphysicsStateVector::setInterfaceZVelocity(double newInterfaceZVelocity)
{
    interfaceZVelocity = newInterfaceZVelocity;
}

void ElasticMultiphysicsStateVector::setMaterial1Density(double newMaterial1Density)
{
    material1Density = newMaterial1Density;
}

void ElasticMultiphysicsStateVector::setMaterial1DistortionTensor(vector<vector<double> > newMaterial1DistortionTensor)
{
    material1DistortionTensor = newMaterial1DistortionTensor;
}

void ElasticMultiphysicsStateVector::setMaterial1Entropy(double newMaterial1Entropy)
{
    material1Entropy = newMaterial1Entropy;
}

void ElasticMultiphysicsStateVector::setMaterial2Density(double newMaterial2Density)
{
    material2Density = newMaterial2Density;
}

void ElasticMultiphysicsStateVector::setMaterial2DistortionTensor(vector<vector<double> > newMaterial2DistortionTensor)
{
    material2DistortionTensor = newMaterial2DistortionTensor;
}

void ElasticMultiphysicsStateVector::setMaterial2Entropy(double newMaterial2Entropy)
{
    material2Entropy = newMaterial2Entropy;
}

double ElasticMultiphysicsStateVector::getMaterial1VolumeFraction()
{
    return material1VolumeFraction;
}

double ElasticMultiphysicsStateVector::getInterfaceXVelocity()
{
    return interfaceXVelocity;
}

double ElasticMultiphysicsStateVector::getInterfaceYVelocity()
{
    return interfaceYVelocity;
}

double ElasticMultiphysicsStateVector::getInterfaceZVelocity()
{
    return interfaceZVelocity;
}

double ElasticMultiphysicsStateVector::getMaterial1Density()
{
    return material1Density;
}

vector<vector<double> > ElasticMultiphysicsStateVector::getMaterial1DistortionTensor()
{
    return material1DistortionTensor;
}

double ElasticMultiphysicsStateVector::getMaterial1Entropy()
{
    return material1Entropy;
}

double ElasticMultiphysicsStateVector::getMaterial2Density()
{
    return material2Density;
}

vector<vector<double> > ElasticMultiphysicsStateVector::getMaterial2DistortionTensor()
{
    return material2DistortionTensor;
}

double ElasticMultiphysicsStateVector::getMaterial2Entropy()
{
    return material2Entropy;
}

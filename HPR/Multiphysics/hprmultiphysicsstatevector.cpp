#include "hprmultiphysicsstatevector.h"

HPRMultiphysicsStateVector::HPRMultiphysicsStateVector()
{
    material1VolumeFraction = 0.999;

    interfaceXVelocity = 0.0;
    interfaceYVelocity = 0.0;
    interfaceZVelocity = 0.0;

    material1Density = 1.0;
    material1Pressure = 1.0 / 1.4;
    material1DistortionTensor = MatrixAlgebra::computeIdentityMatrix(3);

    material1XThermalImpulse = 0.0;
    material1YThermalImpulse = 0.0;
    material1ZThermalImpulse = 0.0;

    material2Density = 1.0;
    material2Pressure = 1.0 / 1.4;
    material2DistortionTensor = MatrixAlgebra::computeIdentityMatrix(3);

    material2XThermalImpulse = 0.0;
    material2YThermalImpulse = 0.0;
    material2ZThermalImpulse = 0.0;
}

HPRMultiphysicsStateVector::HPRMultiphysicsStateVector(double newMaterial1VolumeFraction, double newInterfaceXVelocity, double newInterfaceYVelocity, double newInterfaceZVelocity,
                                                       double newMaterial1Density, double newMaterial1Pressure, vector<vector<double> > newMaterial1DistortionTensor,
                                                       double newMaterial1XThermalImpulse, double newMaterial1YThermalImpulse, double newMaterial1ZThermalImpulse,
                                                       double newMaterial2Density, double newMaterial2Pressure, vector<vector<double> > newMaterial2DistortionTensor,
                                                       double newMaterial2XThermalImpulse, double newMaterial2YThermalImpulse, double newMaterial2ZThermalImpulse)
{
    material1VolumeFraction = newMaterial1VolumeFraction;

    interfaceXVelocity = newInterfaceXVelocity;
    interfaceYVelocity = newInterfaceYVelocity;
    interfaceZVelocity = newInterfaceZVelocity;

    material1Density = newMaterial1Density;
    material1Pressure = newMaterial1Pressure;
    material1DistortionTensor = newMaterial1DistortionTensor;

    material1XThermalImpulse = newMaterial1XThermalImpulse;
    material1YThermalImpulse = newMaterial1YThermalImpulse;
    material1ZThermalImpulse = newMaterial1ZThermalImpulse;

    material2Density = newMaterial2Density;
    material2Pressure = newMaterial2Pressure;
    material2DistortionTensor = newMaterial2DistortionTensor;

    material2XThermalImpulse = newMaterial2XThermalImpulse;
    material2YThermalImpulse = newMaterial2YThermalImpulse;
    material2ZThermalImpulse = newMaterial2ZThermalImpulse;
}

HPRMultiphysicsStateVector::HPRMultiphysicsStateVector(double newMaterial1VolumeFraction, double newInterfaceXVelocity, double newInterfaceYVelocity, double newInterfaceZVelocity,
                                                       vector<vector<double> > newMaterial1DistortionTensor, double newMaterial1Entropy, vector<vector<double> > newMaterial2DistortionTensor,
                                                       double newMaterial2Entropy, HyperelasticMaterialParameters hyperelasticMaterial1Parameters,
                                                       HyperelasticMaterialParameters hyperelasticMaterial2Parameters, HPRMaterialParameters material1Parameters,
                                                       HPRMaterialParameters material2Parameters)
{
    double material1ReferenceMassDensity = hyperelasticMaterial1Parameters.getReferenceMassDensity();
    material1Density = material1ReferenceMassDensity * MatrixAlgebra::computeDeterminant(newMaterial1DistortionTensor);

    double material2ReferenceMassDensity = hyperelasticMaterial2Parameters.getReferenceMassDensity();
    material2Density = material2ReferenceMassDensity * MatrixAlgebra::computeDeterminant(newMaterial2DistortionTensor);

    material1VolumeFraction = newMaterial1VolumeFraction;
    interfaceXVelocity = newInterfaceXVelocity;
    interfaceYVelocity = newInterfaceYVelocity;
    interfaceZVelocity = newInterfaceZVelocity;

    material1DistortionTensor = newMaterial1DistortionTensor;
    material2DistortionTensor = newMaterial2DistortionTensor;

    material1XThermalImpulse = 0.0;
    material1YThermalImpulse = 0.0;
    material1ZThermalImpulse = 0.0;

    material2XThermalImpulse = 0.0;
    material2YThermalImpulse = 0.0;
    material2ZThermalImpulse = 0.0;

    double material1TotalEnergy = ElasticEquationOfState::computeTotalEnergy(material1DistortionTensor, newMaterial1Entropy, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                                             hyperelasticMaterial1Parameters);
    double material2TotalEnergy = ElasticEquationOfState::computeTotalEnergy(material2DistortionTensor, newMaterial2Entropy, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                                             hyperelasticMaterial2Parameters);

    material1Pressure = HPREquationOfState::computePressure(material1Density, material1TotalEnergy, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity, material1DistortionTensor,
                                                            material1XThermalImpulse, material1YThermalImpulse, material1ZThermalImpulse, material1Parameters);
    material2Pressure = HPREquationOfState::computePressure(material2Density, material2TotalEnergy, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity, material2DistortionTensor,
                                                            material2XThermalImpulse, material2YThermalImpulse, material2ZThermalImpulse, material2Parameters);

    double c = 0.0;
}

void HPRMultiphysicsStateVector::setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    material1VolumeFraction = newPrimitiveVariableVector[0];

    interfaceXVelocity = newPrimitiveVariableVector[1];
    interfaceYVelocity = newPrimitiveVariableVector[2];
    interfaceZVelocity = newPrimitiveVariableVector[3];

    material1Density = newPrimitiveVariableVector[4];
    material1Pressure = newPrimitiveVariableVector[5];
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            material1DistortionTensor[i][j] = newPrimitiveVariableVector[6 + (i * 3) + j];
        }
    }

    if (material1Parameters.getIsThermal())
    {
        material1XThermalImpulse = newPrimitiveVariableVector[15];
        material1YThermalImpulse = newPrimitiveVariableVector[16];
        material1ZThermalImpulse = newPrimitiveVariableVector[17];
    }

    material2Density = newPrimitiveVariableVector[18];
    material2Pressure = newPrimitiveVariableVector[19];
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            material2DistortionTensor[i][j] = newPrimitiveVariableVector[20 + (i * 3) + j];
        }
    }

    if (material2Parameters.getIsThermal())
    {
        material2XThermalImpulse = newPrimitiveVariableVector[29];
        material2YThermalImpulse = newPrimitiveVariableVector[30];
        material2ZThermalImpulse = newPrimitiveVariableVector[31];
    }
}

void HPRMultiphysicsStateVector::setConservedVariableVector(vector<double> newConservedVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
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
            material1DistortionTensor[i][j] = newConservedVariableVector[7 + (i * 3) + j] / pow(material1VolumeFraction, (1.0 / 3.0));
        }
    }

    if (material1Parameters.getIsThermal())
    {
        material1XThermalImpulse = newConservedVariableVector[16] / (material1VolumeFraction * material1Density);
        material1YThermalImpulse = newConservedVariableVector[17] / (material1VolumeFraction * material1Density);
        material1ZThermalImpulse = newConservedVariableVector[18] / (material1VolumeFraction * material1Density);
    }

    double computedMaterial1TotalEnergy = newConservedVariableVector[6] / (material1VolumeFraction * material1Density);
    material1Pressure = HPREquationOfState::computePressure(material1Density, computedMaterial1TotalEnergy, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                            material1DistortionTensor, material1XThermalImpulse, material1YThermalImpulse, material1ZThermalImpulse, material1Parameters);

    material2Density = newConservedVariableVector[19] / material2VolumeFraction;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            material2DistortionTensor[i][j] = newConservedVariableVector[21 + (i * 3) + j] / pow(material2VolumeFraction, (1.0 / 3.0));
        }
    }

    if (material2Parameters.getIsThermal())
    {
        material2XThermalImpulse = newConservedVariableVector[30] / (material2VolumeFraction * material2Density);
        material2YThermalImpulse = newConservedVariableVector[31] / (material2VolumeFraction * material2Density);
        material2ZThermalImpulse = newConservedVariableVector[32] / (material2VolumeFraction * material2Density);
    }

    double computedMaterial2TotalEnergy = newConservedVariableVector[20] / (material2VolumeFraction * material2Density);
    material2Pressure = HPREquationOfState::computePressure(material2Density, computedMaterial2TotalEnergy, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                            material2DistortionTensor, material2XThermalImpulse, material2YThermalImpulse, material2ZThermalImpulse, material2Parameters);
}

vector<double> HPRMultiphysicsStateVector::computePrimitiveVariableVector(HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> primitiveVariableVector(32);

    primitiveVariableVector[0] = material1VolumeFraction;

    primitiveVariableVector[1] = interfaceXVelocity;
    primitiveVariableVector[2] = interfaceYVelocity;
    primitiveVariableVector[3] = interfaceZVelocity;

    primitiveVariableVector[4] = material1Density;
    primitiveVariableVector[5] = material1Pressure;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            primitiveVariableVector[6 + (i * 3) + j] = material1DistortionTensor[i][j];
        }
    }

    if (material1Parameters.getIsThermal())
    {
        primitiveVariableVector[15] = material1XThermalImpulse;
        primitiveVariableVector[16] = material1YThermalImpulse;
        primitiveVariableVector[17] = material1ZThermalImpulse;
    }

    primitiveVariableVector[18] = material2Density;
    primitiveVariableVector[19] = material2Pressure;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            primitiveVariableVector[20 + (i * 3) + j] = material2DistortionTensor[i][j];
        }
    }

    if (material2Parameters.getIsThermal())
    {
        primitiveVariableVector[29] = material2XThermalImpulse;
        primitiveVariableVector[30] = material2YThermalImpulse;
        primitiveVariableVector[31] = material2ZThermalImpulse;
    }

    return primitiveVariableVector;
}

vector<double> HPRMultiphysicsStateVector::computeConservedVariableVector(HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> conservedVariableVector(33);

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
            conservedVariableVector[7 + (i * 3) + j] = pow(material1VolumeFraction, (1.0 / 3.0)) * material1DistortionTensor[i][j];
        }
    }

    if (material1Parameters.getIsThermal())
    {
        conservedVariableVector[16] = material1VolumeFraction * (material1Density * material1XThermalImpulse);
        conservedVariableVector[17] = material1VolumeFraction * (material1Density * material1YThermalImpulse);
        conservedVariableVector[18] = material1VolumeFraction * (material1Density * material1ZThermalImpulse);
    }

    double computedMaterial1TotalEnergy = HPREquationOfState::computeTotalEnergy(material1Density, material1Pressure, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                                                 material1DistortionTensor, material1XThermalImpulse, material1YThermalImpulse, material1ZThermalImpulse,
                                                                                 material1Parameters);
    conservedVariableVector[6] = material1VolumeFraction * (material1Density * computedMaterial1TotalEnergy);

    conservedVariableVector[19] = material2VolumeFraction * material2Density;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            conservedVariableVector[21 + (i * 3) + j] = pow(material2VolumeFraction, (1.0 / 3.0)) * material2DistortionTensor[i][j];
        }
    }

    if (material2Parameters.getIsThermal())
    {
        conservedVariableVector[30] = material2VolumeFraction * (material2Density * material2XThermalImpulse);
        conservedVariableVector[31] = material2VolumeFraction * (material2Density * material2YThermalImpulse);
        conservedVariableVector[32] = material2VolumeFraction * (material2Density * material2ZThermalImpulse);
    }

    double computedMaterial2TotalEnergy = HPREquationOfState::computeTotalEnergy(material2Density, material2Pressure, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                                                 material2DistortionTensor, material2XThermalImpulse, material2YThermalImpulse, material2ZThermalImpulse,
                                                                                 material2Parameters);
    conservedVariableVector[20] = material2VolumeFraction * (material2Density * computedMaterial2TotalEnergy);

    return conservedVariableVector;
}

vector<double> HPRMultiphysicsStateVector::computeXFluxVector(vector<double> conservedVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> fluxVector(33);

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
            computedMaterial1DistortionTensor[i][j] = conservedVariableVector[7 + (i * 3) + j] / pow(computedMaterial1VolumeFraction, (1.0 / 3.0));
        }
    }

    double computedMaterial1XThermalImpulse = conservedVariableVector[16] / (computedMaterial1VolumeFraction * computedMaterial1Density);
    double computedMaterial1YThermalImpulse = conservedVariableVector[17] / (computedMaterial1VolumeFraction * computedMaterial1Density);
    double computedMaterial1ZThermalImpulse = conservedVariableVector[18] / (computedMaterial1VolumeFraction * computedMaterial1Density);

    double computedMaterial1TotalEnergy = conservedVariableVector[6] / (computedMaterial1VolumeFraction * computedMaterial1Density);
    double computedMaterial1Pressure = HPREquationOfState::computePressure(computedMaterial1Density, computedMaterial1TotalEnergy, computedInterfaceXVelocity, computedInterfaceYVelocity,
                                                                           computedInterfaceZVelocity, computedMaterial1DistortionTensor, computedMaterial1XThermalImpulse,
                                                                           computedMaterial1YThermalImpulse, computedMaterial1ZThermalImpulse, material1Parameters);

    double computedMaterial2Density = conservedVariableVector[19] / computedMaterial2VolumeFraction;

    vector<vector<double> > computedMaterial2DistortionTensor(3, vector<double>(3));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            computedMaterial2DistortionTensor[i][j] = conservedVariableVector[21 + (i * 3) + j] / pow(computedMaterial2VolumeFraction, (1.0 / 3.0));
        }
    }

    double computedMaterial2XThermalImpulse = conservedVariableVector[30] / (computedMaterial2VolumeFraction * computedMaterial2Density);
    double computedMaterial2YThermalImpulse = conservedVariableVector[31] / (computedMaterial2VolumeFraction * computedMaterial2Density);
    double computedMaterial2ZThermalImpulse = conservedVariableVector[32] / (computedMaterial2VolumeFraction * computedMaterial2Density);

    double computedMaterial2TotalEnergy = conservedVariableVector[20] / (computedMaterial2VolumeFraction * computedMaterial2Density);
    double computedMaterial2Pressure = HPREquationOfState::computePressure(computedMaterial2Density, computedMaterial2TotalEnergy, computedInterfaceXVelocity, computedInterfaceYVelocity,
                                                                           computedInterfaceZVelocity, computedMaterial2DistortionTensor, computedMaterial2XThermalImpulse,
                                                                           computedMaterial2YThermalImpulse, computedMaterial2ZThermalImpulse, material2Parameters);

    double computedInterfacePressure = (computedMaterial1VolumeFraction * computedMaterial1Pressure) + (computedMaterial2VolumeFraction * computedMaterial2Pressure);

    vector<vector<double> > computedMaterial1ShearStressTensor = HPREquationOfState::computeShearStressTensor(computedMaterial1Density, computedMaterial1DistortionTensor, material1Parameters);
    vector<vector<double> > computedMaterial2ShearStressTensor = HPREquationOfState::computeShearStressTensor(computedMaterial2Density, computedMaterial2DistortionTensor, material2Parameters);
    vector<vector<double> > computedInterfaceShearStressTensor = MatrixAlgebra::addMatrices(MatrixAlgebra::multiplyMatrix(computedMaterial1VolumeFraction, computedMaterial1ShearStressTensor),
                                                                                            MatrixAlgebra::multiplyMatrix(computedMaterial2VolumeFraction, computedMaterial2ShearStressTensor));

    for (int i = 0; i < 33; i++)
    {
        fluxVector[i] = 0.0;
    }

    fluxVector[0] = computedTotalDensity * computedInterfaceXVelocity;
    fluxVector[1] = computedTotalDensity * (computedInterfaceXVelocity * computedMaterial1VolumeFraction);

    fluxVector[2] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceXVelocity)) + computedInterfacePressure;
    fluxVector[3] = computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceYVelocity);
    fluxVector[4] = computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceZVelocity);

    fluxVector[5] = computedMaterial1VolumeFraction * (computedMaterial1Density * computedInterfaceXVelocity);
    fluxVector[6] = computedMaterial1VolumeFraction * ((computedMaterial1Density * (computedInterfaceXVelocity * computedMaterial1TotalEnergy)) +
                                                       (computedMaterial1Pressure * computedInterfaceXVelocity));

    fluxVector[19] = computedMaterial2VolumeFraction * (computedMaterial2Density * computedInterfaceXVelocity);
    fluxVector[20] = computedMaterial2VolumeFraction * ((computedMaterial2Density * (computedInterfaceXVelocity * computedMaterial2TotalEnergy)) +
                                                        (computedMaterial2Pressure * computedInterfaceXVelocity));

    vector<double> computedVelocityVector(3);
    computedVelocityVector[0] = computedInterfaceXVelocity;
    computedVelocityVector[1] = computedInterfaceYVelocity;
    computedVelocityVector[2] = computedInterfaceZVelocity;

    vector<double> computedMaterial1DistortionTensorVelocityVectorProduct = MatrixAlgebra::multiplyMatrixByVector(computedMaterial1DistortionTensor, computedVelocityVector);
    vector<double> computedMaterial2DistortionTensorVelocityVectorProduct = MatrixAlgebra::multiplyMatrixByVector(computedMaterial2DistortionTensor, computedVelocityVector);

    fluxVector[6] -= computedMaterial1VolumeFraction * VectorAlgebra::computeDotProduct(computedMaterial1ShearStressTensor[0], computedVelocityVector);
    fluxVector[20] -= computedMaterial2VolumeFraction * VectorAlgebra::computeDotProduct(computedMaterial2ShearStressTensor[0], computedVelocityVector);

    for (int i = 0; i < 3; i++)
    {
        fluxVector[2 + i] -= computedInterfaceShearStressTensor[0][i];
    }

    for (int i = 0; i < 3; i++)
    {
        fluxVector[7 + (i * 3)] = pow(computedMaterial1VolumeFraction, (1.0 / 3.0)) * computedMaterial1DistortionTensorVelocityVectorProduct[i];
    }
    for (int i = 0; i < 3; i++)
    {
        fluxVector[21 + (i * 3)] = pow(computedMaterial2VolumeFraction, (1.0 / 3.0)) * computedMaterial2DistortionTensorVelocityVectorProduct[i];
    }

    if (material1Parameters.getIsThermal())
    {
        double computedMaterial1Temperature = HPREquationOfState::computeTemperature(computedMaterial1Density, computedMaterial1Pressure, material1Parameters);
        vector<double> computedMaterial1HeatFluxVector = HPREquationOfState::computeHeatFluxVector(computedMaterial1Temperature, computedMaterial1XThermalImpulse,
                                                                                                   computedMaterial1YThermalImpulse, computedMaterial1ZThermalImpulse, material1Parameters);

        fluxVector[6] += computedMaterial1VolumeFraction * computedMaterial1HeatFluxVector[0];

        fluxVector[16] = computedMaterial1VolumeFraction * ((computedMaterial1Density * (computedInterfaceXVelocity * computedMaterial1XThermalImpulse)) + computedMaterial1Temperature);
        fluxVector[17] = computedMaterial1VolumeFraction * (computedMaterial1Density * (computedInterfaceXVelocity * computedMaterial1YThermalImpulse));
        fluxVector[18] = computedMaterial1VolumeFraction * (computedMaterial1Density * (computedInterfaceXVelocity * computedMaterial1ZThermalImpulse));
    }

    if (material2Parameters.getIsThermal())
    {
        double computedMaterial2Temperature = HPREquationOfState::computeTemperature(computedMaterial2Density, computedMaterial2Pressure, material2Parameters);
        vector<double> computedMaterial2HeatFluxVector = HPREquationOfState::computeHeatFluxVector(computedMaterial2Temperature, computedMaterial2XThermalImpulse,
                                                                                                   computedMaterial2YThermalImpulse, computedMaterial2ZThermalImpulse, material2Parameters);

        fluxVector[20] += computedMaterial2VolumeFraction * computedMaterial2HeatFluxVector[0];

        fluxVector[30] = computedMaterial2VolumeFraction * ((computedMaterial2Density * (computedInterfaceXVelocity * computedMaterial2XThermalImpulse)) + computedMaterial2Temperature);
        fluxVector[31] = computedMaterial2VolumeFraction * (computedMaterial2Density * (computedInterfaceXVelocity * computedMaterial2YThermalImpulse));
        fluxVector[32] = computedMaterial2VolumeFraction * (computedMaterial2Density * (computedInterfaceXVelocity * computedMaterial2ZThermalImpulse));
    }

    return fluxVector;
}

vector<double> HPRMultiphysicsStateVector::computeXFluxVector(HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    return computeXFluxVector(computeConservedVariableVector(material1Parameters, material2Parameters), material1Parameters, material2Parameters);
}

double HPRMultiphysicsStateVector::computeMaterial1TotalEnergy(HPRMaterialParameters material1Parameters)
{
    return HPREquationOfState::computeTotalEnergy(material1Density, material1Pressure, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity, material1DistortionTensor,
                                                  material1XThermalImpulse, material1YThermalImpulse, material1ZThermalImpulse, material1Parameters);
}

double HPRMultiphysicsStateVector::computeMaterial1Temperature(HPRMaterialParameters material1Parameters)
{
    return HPREquationOfState::computeTemperature(material1Density, material1Pressure, material1Parameters);
}

vector<double> HPRMultiphysicsStateVector::computeMaterial1HeatFluxVector(HPRMaterialParameters material1Parameters)
{
    return HPREquationOfState::computeHeatFluxVector(computeMaterial1Temperature(material1Parameters), material1XThermalImpulse, material1YThermalImpulse, material1ZThermalImpulse,
                                                     material1Parameters);
}

double HPRMultiphysicsStateVector::computeMaterial1TotalEnergyDerivativeDensity(HPRMaterialParameters material1Parameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativeDensity(material1Density, material1Pressure, material1Parameters);
}

double HPRMultiphysicsStateVector::computeMaterial1TotalEnergyDerivativePressure(HPRMaterialParameters material1Parameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativePressure(material1Density, material1Parameters);
}

vector<vector<double> > HPRMultiphysicsStateVector::computeMaterial1TotalEnergyDerivativeDistortionTensor(HPRMaterialParameters material1Parameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativeDistortionTensor(material1DistortionTensor, material1Parameters);
}

vector<double> HPRMultiphysicsStateVector::computeMaterial1TotalEnergyDerivativeThermalImpulse(HPRMaterialParameters material1Parameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativeThermalImpulse(material1XThermalImpulse, material1YThermalImpulse, material1ZThermalImpulse, material1Parameters);
}

double HPRMultiphysicsStateVector::computeMaterial1TemperatureDerivativeDensity(HPRMaterialParameters material1Parameters)
{
    return HPRDerivatives::computeTemperatureDerivativeDensity(material1Density, material1Pressure, material1Parameters);
}

double HPRMultiphysicsStateVector::computeMaterial1TemperatureDerivativePressure(HPRMaterialParameters material1Parameters)
{
    return HPRDerivatives::computeTemperatureDerivativePressure(material1Density, material1Parameters);
}

double HPRMultiphysicsStateVector::computeMaterial1Theta1Reciprocal(HPRMaterialParameters material1Parameters)
{
    return HPRSourceTerms::computeTheta1Reciprocal(material1Density, material1DistortionTensor, material1Parameters);
}

double HPRMultiphysicsStateVector::computeMaterial1Theta2Reciprocal(HPRMaterialParameters material1Parameters)
{
    return HPRSourceTerms::computeTheta2Reciprocal(material1Density, computeMaterial1Temperature(material1Parameters), material1Parameters);
}

vector<vector<double> > HPRMultiphysicsStateVector::computeMaterial1ShearStressTensor(HPRMaterialParameters material1Parameters)
{
    return HPREquationOfState::computeShearStressTensor(material1Density, material1DistortionTensor, material1Parameters);
}

vector<vector<double> > HPRMultiphysicsStateVector::computeMaterial1ShearStressTensorDerivativeDensity(HPRMaterialParameters material1Parameters)
{
    return HPREquationOfState::computeShearStressTensorDerivativeDensity(material1DistortionTensor, material1Parameters);
}

vector<vector<vector<vector<double> > > > HPRMultiphysicsStateVector::computeMaterial1ShearStressTensorDerivativeDistortionTensor(HPRMaterialParameters material1Parameters)
{
    return HPREquationOfState::computeShearStressTensorDerivativeDistortionTensor(material1Density, material1DistortionTensor, material1Parameters);
}

double HPRMultiphysicsStateVector::computeMaterial2TotalEnergy(HPRMaterialParameters material2Parameters)
{
    return HPREquationOfState::computeTotalEnergy(material2Density, material2Pressure, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity, material2DistortionTensor,
                                                  material2XThermalImpulse, material2YThermalImpulse, material2ZThermalImpulse, material2Parameters);
}

double HPRMultiphysicsStateVector::computeMaterial2Temperature(HPRMaterialParameters material2Parameters)
{
    return HPREquationOfState::computeTemperature(material2Density, material2Pressure, material2Parameters);
}

vector<double> HPRMultiphysicsStateVector::computeMaterial2HeatFluxVector(HPRMaterialParameters material2Parameters)
{
    return HPREquationOfState::computeHeatFluxVector(computeMaterial2Temperature(material2Parameters), material2XThermalImpulse, material2YThermalImpulse, material2ZThermalImpulse,
                                                     material2Parameters);
}

double HPRMultiphysicsStateVector::computeMaterial2TotalEnergyDerivativeDensity(HPRMaterialParameters material2Parameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativeDensity(material2Density, material2Pressure, material2Parameters);
}

double HPRMultiphysicsStateVector::computeMaterial2TotalEnergyDerivativePressure(HPRMaterialParameters material2Parameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativePressure(material2Density, material2Parameters);
}

vector<vector<double> > HPRMultiphysicsStateVector::computeMaterial2TotalEnergyDerivativeDistortionTensor(HPRMaterialParameters material2Parameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativeDistortionTensor(material2DistortionTensor, material2Parameters);
}

vector<double> HPRMultiphysicsStateVector::computeMaterial2TotalEnergyDerivativeThermalImpulse(HPRMaterialParameters material2Parameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativeThermalImpulse(material2XThermalImpulse, material2YThermalImpulse, material2ZThermalImpulse, material2Parameters);
}

double HPRMultiphysicsStateVector::computeMaterial2TemperatureDerivativeDensity(HPRMaterialParameters material2Parameters)
{
    return HPRDerivatives::computeTemperatureDerivativeDensity(material2Density, material2Pressure, material2Parameters);
}

double HPRMultiphysicsStateVector::computeMaterial2TemperatureDerivativePressure(HPRMaterialParameters material2Parameters)
{
    return HPRDerivatives::computeTemperatureDerivativePressure(material2Density, material2Parameters);
}

double HPRMultiphysicsStateVector::computeMaterial2Theta1Reciprocal(HPRMaterialParameters material2Parameters)
{
    return HPRSourceTerms::computeTheta1Reciprocal(material2Density, material2DistortionTensor, material2Parameters);
}

double HPRMultiphysicsStateVector::computeMaterial2Theta2Reciprocal(HPRMaterialParameters material2Parameters)
{
    return HPRSourceTerms::computeTheta2Reciprocal(material2Density, computeMaterial2Temperature(material2Parameters), material2Parameters);
}

vector<vector<double> > HPRMultiphysicsStateVector::computeMaterial2ShearStressTensor(HPRMaterialParameters material2Parameters)
{
    return HPREquationOfState::computeShearStressTensor(material2Density, material2DistortionTensor, material2Parameters);
}

vector<vector<double> > HPRMultiphysicsStateVector::computeMaterial2ShearStressTensorDerivativeDensity(HPRMaterialParameters material2Parameters)
{
    return HPREquationOfState::computeShearStressTensorDerivativeDensity(material2DistortionTensor, material2Parameters);
}

vector<vector<vector<vector<double> > > > HPRMultiphysicsStateVector::computeMaterial2ShearStressTensorDerivativeDistortionTensor(HPRMaterialParameters material2Parameters)
{
    return HPREquationOfState::computeShearStressTensorDerivativeDistortionTensor(material2Density, material2DistortionTensor, material2Parameters);
}

double HPRMultiphysicsStateVector::computeTotalDensity()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1Density) + (material2VolumeFraction * material2Density);
}

double HPRMultiphysicsStateVector::computeTotalPressure()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1Pressure) + (material2VolumeFraction * material2Pressure);
}

vector<vector<double> > HPRMultiphysicsStateVector::computeTotalDistortionTensor()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return MatrixAlgebra::addMatrices(MatrixAlgebra::multiplyMatrix(pow(material1VolumeFraction, (1.0 / 3.0)), material1DistortionTensor),
                                      MatrixAlgebra::multiplyMatrix(pow(material2VolumeFraction, (1.0 / 3.0)), material2DistortionTensor));
}

double HPRMultiphysicsStateVector::computeTotalXThermalImpulse()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1XThermalImpulse) + (material2VolumeFraction * material2XThermalImpulse);
}

double HPRMultiphysicsStateVector::computeTotalYThermalImpulse()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1YThermalImpulse) + (material2VolumeFraction * material2YThermalImpulse);
}

double HPRMultiphysicsStateVector::computeTotalZThermalImpulse()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1ZThermalImpulse) + (material2VolumeFraction * material2ZThermalImpulse);
}

void HPRMultiphysicsStateVector::relaxTotalDensity()
{
    double totalDensity = computeTotalDensity();

    material1Density = totalDensity;
    material2Density = totalDensity;
}

void HPRMultiphysicsStateVector::relaxTotalPressure()
{
    double totalPressure = computeTotalPressure();

    material1Pressure = totalPressure;
    material2Pressure = totalPressure;
}

void HPRMultiphysicsStateVector::relaxTotalDistortionTensor()
{
    vector<vector<double> > totalDistortionTensor = computeTotalDistortionTensor();

    material1DistortionTensor = totalDistortionTensor;
    material2DistortionTensor = totalDistortionTensor;
}

void HPRMultiphysicsStateVector::relaxTotalXThermalImpulse()
{
    double totalXThermalImpulse = computeTotalXThermalImpulse();

    material1XThermalImpulse = totalXThermalImpulse;
    material2XThermalImpulse = totalXThermalImpulse;
}

void HPRMultiphysicsStateVector::relaxTotalYThermalImpulse()
{
    double totalYThermalImpulse = computeTotalYThermalImpulse();

    material1YThermalImpulse = totalYThermalImpulse;
    material2YThermalImpulse = totalYThermalImpulse;
}

void HPRMultiphysicsStateVector::relaxTotalZThermalImpulse()
{
    double totalZThermalImpulse = computeTotalZThermalImpulse();

    material1ZThermalImpulse = totalZThermalImpulse;
    material2ZThermalImpulse = totalZThermalImpulse;
}


void HPRMultiphysicsStateVector::setMaterial1VolumeFraction(double newMaterial1VolumeFraction)
{
    material1VolumeFraction = newMaterial1VolumeFraction;
}

void HPRMultiphysicsStateVector::setInterfaceXVelocity(double newInterfaceXVelocity)
{
    interfaceXVelocity = newInterfaceXVelocity;
}

void HPRMultiphysicsStateVector::setInterfaceYVelocity(double newInterfaceYVelocity)
{
    interfaceYVelocity = newInterfaceYVelocity;
}

void HPRMultiphysicsStateVector::setInterfaceZVelocity(double newInterfaceZVelocity)
{
    interfaceZVelocity = newInterfaceZVelocity;
}

void HPRMultiphysicsStateVector::setMaterial1Density(double newMaterial1Density)
{
    material1Density = newMaterial1Density;
}

void HPRMultiphysicsStateVector::setMaterial1Pressure(double newMaterial1Pressure)
{
    material1Pressure = newMaterial1Pressure;
}

void HPRMultiphysicsStateVector::setMaterial1DistortionTensor(vector<vector<double> > newMaterial1DistortionTensor)
{
    material1DistortionTensor = newMaterial1DistortionTensor;
}

void HPRMultiphysicsStateVector::setMaterial1XThermalImpulse(double newMaterial1XThermalImpulse)
{
    material1XThermalImpulse = newMaterial1XThermalImpulse;
}

void HPRMultiphysicsStateVector::setMaterial1YThermalImpulse(double newMaterial1YThermalImpulse)
{
    material1YThermalImpulse = newMaterial1YThermalImpulse;
}

void HPRMultiphysicsStateVector::setMaterial1ZThermalImpulse(double newMaterial1ZThermalImpulse)
{
    material1ZThermalImpulse = newMaterial1ZThermalImpulse;
}

void HPRMultiphysicsStateVector::setMaterial2Density(double newMaterial2Density)
{
    material2Density = newMaterial2Density;
}

void HPRMultiphysicsStateVector::setMaterial2Pressure(double newMaterial2Pressure)
{
    material2Pressure = newMaterial2Pressure;
}

void HPRMultiphysicsStateVector::setMaterial2DistortionTensor(vector<vector<double> > newMaterial2DistortionTensor)
{
    material2DistortionTensor = newMaterial2DistortionTensor;
}

void HPRMultiphysicsStateVector::setMaterial2XThermalImpulse(double newMaterial2XThermalImpulse)
{
    material2XThermalImpulse = newMaterial2XThermalImpulse;
}

void HPRMultiphysicsStateVector::setMaterial2YThermalImpulse(double newMaterial2YThermalImpulse)
{
    material2YThermalImpulse = newMaterial2YThermalImpulse;
}

void HPRMultiphysicsStateVector::setMaterial2ZThermalImpulse(double newMaterial2ZThermalImpulse)
{
    material2ZThermalImpulse = newMaterial2ZThermalImpulse;
}

double HPRMultiphysicsStateVector::getMaterial1VolumeFraction()
{
    return material1VolumeFraction;
}

double HPRMultiphysicsStateVector::getInterfaceXVelocity()
{
    return interfaceXVelocity;
}

double HPRMultiphysicsStateVector::getInterfaceYVelocity()
{
    return interfaceYVelocity;
}

double HPRMultiphysicsStateVector::getInterfaceZVelocity()
{
    return interfaceZVelocity;
}

double HPRMultiphysicsStateVector::getMaterial1Density()
{
    return material1Density;
}

double HPRMultiphysicsStateVector::getMaterial1Pressure()
{
    return material1Pressure;
}

vector<vector<double> > HPRMultiphysicsStateVector::getMaterial1DistortionTensor()
{
    return material1DistortionTensor;
}

double HPRMultiphysicsStateVector::getMaterial1XThermalImpulse()
{
    return material1XThermalImpulse;
}

double HPRMultiphysicsStateVector::getMaterial1YThermalImpulse()
{
    return material1YThermalImpulse;
}

double HPRMultiphysicsStateVector::getMaterial1ZThermalImpulse()
{
    return material1ZThermalImpulse;
}

double HPRMultiphysicsStateVector::getMaterial2Density()
{
    return material2Density;
}

double HPRMultiphysicsStateVector::getMaterial2Pressure()
{
    return material2Pressure;
}

vector<vector<double> > HPRMultiphysicsStateVector::getMaterial2DistortionTensor()
{
    return material2DistortionTensor;
}

double HPRMultiphysicsStateVector::getMaterial2XThermalImpulse()
{
    return material2XThermalImpulse;
}

double HPRMultiphysicsStateVector::getMaterial2YThermalImpulse()
{
    return material2YThermalImpulse;
}

double HPRMultiphysicsStateVector::getMaterial2ZThermalImpulse()
{
    return material2ZThermalImpulse;
}

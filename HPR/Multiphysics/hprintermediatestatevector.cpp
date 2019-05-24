#include "hprintermediatestatevector.h"

HPRIntermediateStateVector::HPRIntermediateStateVector()
{
    material1VolumeFraction = 0.999;
    interfaceDistortionTensor = MatrixAlgebra::computeIdentityMatrix(3);

    interfaceXVelocity = 0.0;
    interfaceYVelocity = 0.0;
    interfaceZVelocity = 0.0;

    interfaceXThermalImpulse = 0.0;
    interfaceYThermalImpulse = 0.0;
    interfaceZThermalImpulse = 0.0;

    material1Density = 1.0;
    material1Pressure = 1.0 / 1.4;

    material2Density = 1.0;
    material2Pressure = 1.0 / 1.4;
}

HPRIntermediateStateVector::HPRIntermediateStateVector(double newMaterial1VolumeFraction, vector<vector<double> > newInterfaceDistortionTensor, double newInterfaceXVelocity,
                                                       double newInterfaceYVelocity, double newInterfaceZVelocity, double newInterfaceXThermalImpulse, double newInterfaceYThermalImpulse,
                                                       double newInterfaceZThermalImpulse, double newMaterial1Density, double newMaterial1Pressure, double newMaterial2Density,
                                                       double newMaterial2Pressure)
{
    material1VolumeFraction = newMaterial1VolumeFraction;
    interfaceDistortionTensor = newInterfaceDistortionTensor;

    interfaceXVelocity = newInterfaceXVelocity;
    interfaceYVelocity = newInterfaceYVelocity;
    interfaceZVelocity = newInterfaceZVelocity;

    interfaceXThermalImpulse = newInterfaceXThermalImpulse;
    interfaceYThermalImpulse = newInterfaceYThermalImpulse;
    interfaceZThermalImpulse = newInterfaceZThermalImpulse;

    material1Density = newMaterial1Density;
    material1Pressure = newMaterial1Pressure;

    material2Density = newMaterial2Density;
    material2Pressure = newMaterial2Pressure;
}

HPRIntermediateStateVector::HPRIntermediateStateVector(double newMaterial1VolumeFraction, double newInterfaceXVelocity, double newInterfaceYVelocity, double newInterfaceZVelocity,
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

    interfaceXThermalImpulse = 0.0;
    interfaceYThermalImpulse = 0.0;
    interfaceZThermalImpulse = 0.0;

    double material2VolumeFraction = 1.0 - material1VolumeFraction;
    interfaceDistortionTensor = MatrixAlgebra::addMatrices(MatrixAlgebra::multiplyMatrix(pow(material1VolumeFraction, (1.0 / 3.0)), newMaterial1DistortionTensor),
                                                           MatrixAlgebra::multiplyMatrix(pow(material2VolumeFraction, (1.0 / 3.0)), newMaterial2DistortionTensor));

    double material1TotalEnergy = ElasticEquationOfState::computeTotalEnergy(newMaterial1DistortionTensor, newMaterial1Entropy, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                                             hyperelasticMaterial1Parameters);
    double material2TotalEnergy = ElasticEquationOfState::computeTotalEnergy(newMaterial2DistortionTensor, newMaterial2Entropy, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                                             hyperelasticMaterial2Parameters);

    material1Pressure = HPREquationOfState::computePressure(material1Density, material1TotalEnergy, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity, newMaterial1DistortionTensor,
                                                            interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse, material1Parameters);
    material2Pressure = HPREquationOfState::computePressure(material2Density, material2TotalEnergy, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity, newMaterial2DistortionTensor,
                                                            interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse, material2Parameters);
}

void HPRIntermediateStateVector::setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    material1VolumeFraction = newPrimitiveVariableVector[0];

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            interfaceDistortionTensor[i][j] = newPrimitiveVariableVector[1 + (i * 3) + j];
        }
    }

    interfaceXVelocity = newPrimitiveVariableVector[10];
    interfaceYVelocity = newPrimitiveVariableVector[11];
    interfaceZVelocity = newPrimitiveVariableVector[12];

    if (material1Parameters.getIsThermal() || material2Parameters.getIsThermal())
    {
        interfaceXThermalImpulse = newPrimitiveVariableVector[13];
        interfaceYThermalImpulse = newPrimitiveVariableVector[14];
        interfaceZThermalImpulse = newPrimitiveVariableVector[15];
    }

    material1Density = newPrimitiveVariableVector[16];
    material1Pressure = newPrimitiveVariableVector[17];

    material2Density = newPrimitiveVariableVector[18];
    material2Pressure = newPrimitiveVariableVector[19];
}

void HPRIntermediateStateVector::setConservedVariableVector(vector<double> newConservedVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
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

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            interfaceDistortionTensor[i][j] = newConservedVariableVector[2 + (i * 3) + j];
        }
    }

    interfaceXVelocity = newConservedVariableVector[11] / totalDensity;
    interfaceYVelocity = newConservedVariableVector[12] / totalDensity;
    interfaceZVelocity = newConservedVariableVector[13] / totalDensity;

    if (material1Parameters.getIsThermal() || material2Parameters.getIsThermal())
    {
        interfaceXThermalImpulse = newConservedVariableVector[14] / totalDensity;
        interfaceYThermalImpulse = newConservedVariableVector[15] / totalDensity;
        interfaceZThermalImpulse = newConservedVariableVector[16] / totalDensity;
    }

    material1Density = newConservedVariableVector[17] / material1VolumeFraction;
    material2Density = newConservedVariableVector[19] / material2VolumeFraction;

    double computedMaterial1TotalEnergy = newConservedVariableVector[18] / (material1VolumeFraction * material1Density);
    material1Pressure = HPREquationOfState::computePressure(material1Density, computedMaterial1TotalEnergy, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                            interfaceDistortionTensor, interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse, material1Parameters);

    double computedMaterial2TotalEnergy = newConservedVariableVector[20] / (material2VolumeFraction * material2Density);
    material2Pressure = HPREquationOfState::computePressure(material2Density, computedMaterial2TotalEnergy, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                            interfaceDistortionTensor, interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse, material2Parameters);
}

vector<double> HPRIntermediateStateVector::computePrimitiveVariableVector(HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> primitiveVariableVector(20);

    primitiveVariableVector[0] = material1VolumeFraction;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            primitiveVariableVector[1 + (i * 3) + j] = interfaceDistortionTensor[i][j];
        }
    }

    primitiveVariableVector[10] = interfaceXVelocity;
    primitiveVariableVector[11] = interfaceYVelocity;
    primitiveVariableVector[12] = interfaceZVelocity;

    if (material1Parameters.getIsThermal() || material2Parameters.getIsThermal())
    {
        primitiveVariableVector[13] = interfaceXThermalImpulse;
        primitiveVariableVector[14] = interfaceYThermalImpulse;
        primitiveVariableVector[15] = interfaceZThermalImpulse;
    }

    primitiveVariableVector[16] = material1Density;
    primitiveVariableVector[17] = material1Pressure;

    primitiveVariableVector[18] = material2Density;
    primitiveVariableVector[19] = material2Pressure;

    return primitiveVariableVector;
}

vector<double> HPRIntermediateStateVector::computeConservedVariableVector(HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> conservedVariableVector(21);

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

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            conservedVariableVector[2 + (i * 3) + j] = interfaceDistortionTensor[i][j];
        }
    }

    conservedVariableVector[11] = totalDensity * interfaceXVelocity;
    conservedVariableVector[12] = totalDensity * interfaceYVelocity;
    conservedVariableVector[13] = totalDensity * interfaceZVelocity;

    if (material1Parameters.getIsThermal() || material2Parameters.getIsThermal())
    {
        conservedVariableVector[14] = totalDensity * interfaceXThermalImpulse;
        conservedVariableVector[15] = totalDensity * interfaceYThermalImpulse;
        conservedVariableVector[16] = totalDensity * interfaceZThermalImpulse;
    }

    conservedVariableVector[17] = material1VolumeFraction * material1Density;
    conservedVariableVector[19] = material2VolumeFraction * material2Density;

    double computedMaterial1TotalEnergy = HPREquationOfState::computeTotalEnergy(material1Density, material1Pressure, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                                                 interfaceDistortionTensor, interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse,
                                                                                 material1Parameters);
    conservedVariableVector[18] = material1VolumeFraction * (material1Density * computedMaterial1TotalEnergy);

    double computedMaterial2TotalEnergy = HPREquationOfState::computeTotalEnergy(material2Density, material2Pressure, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                                                 interfaceDistortionTensor, interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse,
                                                                                 material2Parameters);
    conservedVariableVector[20] = material2VolumeFraction * (material2Density * computedMaterial2TotalEnergy);

    return conservedVariableVector;
}

vector<double> HPRIntermediateStateVector::computeXFluxVector(vector<double> conservedVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> fluxVector(21);

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

    vector<vector<double> > computedInterfaceDistortionTensor(3, vector<double>(3));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            computedInterfaceDistortionTensor[i][j] = conservedVariableVector[2 + (i * 3) + j];
        }
    }

    double computedInterfaceXVelocity = conservedVariableVector[11] / computedTotalDensity;
    double computedInterfaceYVelocity = conservedVariableVector[12] / computedTotalDensity;
    double computedInterfaceZVelocity = conservedVariableVector[13] / computedTotalDensity;

    double computedInterfaceXThermalImpulse = conservedVariableVector[14] / computedTotalDensity;
    double computedInterfaceYThermalImpulse = conservedVariableVector[15] / computedTotalDensity;
    double computedInterfaceZThermalImpulse = conservedVariableVector[16] / computedTotalDensity;

    double computedMaterial1Density = conservedVariableVector[17] / computedMaterial1VolumeFraction;
    double computedMaterial2Density = conservedVariableVector[19] / computedMaterial2VolumeFraction;

    double computedMaterial1TotalEnergy = conservedVariableVector[18] / (computedMaterial1VolumeFraction * computedMaterial1Density);
    double computedMaterial1Pressure = HPREquationOfState::computePressure(computedMaterial1Density, computedMaterial1TotalEnergy, computedInterfaceXVelocity, computedInterfaceYVelocity,
                                                                           computedInterfaceZVelocity, computedInterfaceDistortionTensor, computedInterfaceXThermalImpulse,
                                                                           computedInterfaceYThermalImpulse, computedInterfaceZThermalImpulse, material1Parameters);

    double computedMaterial2TotalEnergy = conservedVariableVector[20] / (computedMaterial2VolumeFraction * computedMaterial2Density);
    double computedMaterial2Pressure = HPREquationOfState::computePressure(computedMaterial2Density, computedMaterial2TotalEnergy, computedInterfaceXVelocity, computedInterfaceYVelocity,
                                                                           computedInterfaceZVelocity, computedInterfaceDistortionTensor, computedInterfaceXThermalImpulse,
                                                                           computedInterfaceYThermalImpulse, computedInterfaceZThermalImpulse, material2Parameters);

    double computedInterfacePressure = (computedMaterial1VolumeFraction * computedMaterial1Pressure) + (computedMaterial2VolumeFraction * computedMaterial2Pressure);

    vector<vector<double> > computedMaterial1ShearStressTensor = HPREquationOfState::computeShearStressTensor(computedMaterial1Density, computedInterfaceDistortionTensor, material1Parameters);
    vector<vector<double> > computedMaterial2ShearStressTensor = HPREquationOfState::computeShearStressTensor(computedMaterial2Density, computedInterfaceDistortionTensor, material2Parameters);
    vector<vector<double> > computedInterfaceShearStressTensor = MatrixAlgebra::addMatrices(MatrixAlgebra::multiplyMatrix(computedMaterial1VolumeFraction, computedMaterial1ShearStressTensor),
                                                                                            MatrixAlgebra::multiplyMatrix(computedMaterial2VolumeFraction, computedMaterial2ShearStressTensor));

    for (int i = 0; i < 21; i++)
    {
        fluxVector[i] = 0.0;
    }

    fluxVector[0] = computedTotalDensity * computedInterfaceXVelocity;
    fluxVector[1] = computedTotalDensity * (computedInterfaceXVelocity * computedMaterial1VolumeFraction);

    fluxVector[18] = computedMaterial1VolumeFraction * ((computedMaterial1Density * (computedInterfaceXVelocity * computedMaterial1TotalEnergy)) +
                                                        (computedMaterial1Pressure * computedInterfaceXVelocity));
    fluxVector[20] = computedMaterial2VolumeFraction * ((computedMaterial2Density * (computedInterfaceXVelocity * computedMaterial2TotalEnergy)) +
                                                        (computedMaterial2Pressure * computedInterfaceXVelocity));

    fluxVector[11] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceXVelocity)) + computedInterfacePressure;
    fluxVector[12] = computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceYVelocity);
    fluxVector[13] = computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceZVelocity);

    vector<double> computedVelocityVector(3);
    computedVelocityVector[0] = computedInterfaceXVelocity;
    computedVelocityVector[1] = computedInterfaceYVelocity;
    computedVelocityVector[2] = computedInterfaceZVelocity;
    vector<double> computedDistortionTensorVelocityVectorProduct = MatrixAlgebra::multiplyMatrixByVector(computedInterfaceDistortionTensor, computedVelocityVector);

    fluxVector[18] -= computedMaterial1VolumeFraction * VectorAlgebra::computeDotProduct(computedMaterial1ShearStressTensor[0], computedVelocityVector);
    fluxVector[20] -= computedMaterial2VolumeFraction * VectorAlgebra::computeDotProduct(computedMaterial2ShearStressTensor[0], computedVelocityVector);

    for (int i = 0; i < 3; i++)
    {
        fluxVector[11 + i] -= computedInterfaceShearStressTensor[0][i];
    }

    for (int i = 0; i < 3; i++)
    {
        fluxVector[2 + (i * 3)] = computedDistortionTensorVelocityVectorProduct[i];
    }

    if (material1Parameters.getIsThermal() || material2Parameters.getIsThermal())
    {
        double computedMaterial1Temperature = HPREquationOfState::computeTemperature(computedMaterial1Density, computedMaterial1Pressure, material1Parameters);
        vector<double> computedMaterial1HeatFluxVector = HPREquationOfState::computeHeatFluxVector(computedMaterial1Temperature, computedInterfaceXThermalImpulse,
                                                                                                   computedInterfaceYThermalImpulse, computedInterfaceZThermalImpulse, material1Parameters);

        double computedMaterial2Temperature = HPREquationOfState::computeTemperature(computedMaterial2Density, computedMaterial2Pressure, material2Parameters);
        vector<double> computedMaterial2HeatFluxVector = HPREquationOfState::computeHeatFluxVector(computedMaterial2Temperature, computedInterfaceXThermalImpulse,
                                                                                                   computedInterfaceYThermalImpulse, computedInterfaceZThermalImpulse, material2Parameters);

        double computedInterfaceTemperature = (computedMaterial1VolumeFraction * computedMaterial1Temperature) + (computedMaterial2VolumeFraction * computedMaterial2Temperature);

        fluxVector[18] += computedMaterial1VolumeFraction * computedMaterial1HeatFluxVector[0];
        fluxVector[20] += computedMaterial2VolumeFraction * computedMaterial2HeatFluxVector[0];

        fluxVector[14] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceXThermalImpulse)) + computedInterfaceTemperature;
        fluxVector[15] = computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceYThermalImpulse);
        fluxVector[16] = computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceZThermalImpulse);
    }

    fluxVector[17] = computedMaterial1VolumeFraction * (computedMaterial1Density * computedInterfaceXVelocity);
    fluxVector[19] = computedMaterial2VolumeFraction * (computedMaterial2Density * computedInterfaceXVelocity);

    return fluxVector;
}

vector<double> HPRIntermediateStateVector::computeXFluxVector(HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    return computeXFluxVector(computeConservedVariableVector(material1Parameters, material2Parameters), material1Parameters, material2Parameters);
}

vector<double> HPRIntermediateStateVector::computeYFluxVector(vector<double> conservedVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> fluxVector(21);

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

    vector<vector<double> > computedInterfaceDistortionTensor(3, vector<double>(3));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            computedInterfaceDistortionTensor[i][j] = conservedVariableVector[2 + (i * 3) + j];
        }
    }

    double computedInterfaceXVelocity = conservedVariableVector[11] / computedTotalDensity;
    double computedInterfaceYVelocity = conservedVariableVector[12] / computedTotalDensity;
    double computedInterfaceZVelocity = conservedVariableVector[13] / computedTotalDensity;

    double computedInterfaceXThermalImpulse = conservedVariableVector[14] / computedTotalDensity;
    double computedInterfaceYThermalImpulse = conservedVariableVector[15] / computedTotalDensity;
    double computedInterfaceZThermalImpulse = conservedVariableVector[16] / computedTotalDensity;

    double computedMaterial1Density = conservedVariableVector[17] / computedMaterial1VolumeFraction;
    double computedMaterial2Density = conservedVariableVector[19] / computedMaterial2VolumeFraction;

    double computedMaterial1TotalEnergy = conservedVariableVector[18] / (computedMaterial1VolumeFraction * computedMaterial1Density);
    double computedMaterial1Pressure = HPREquationOfState::computePressure(computedMaterial1Density, computedMaterial1TotalEnergy, computedInterfaceXVelocity, computedInterfaceYVelocity,
                                                                           computedInterfaceZVelocity, computedInterfaceDistortionTensor, computedInterfaceXThermalImpulse,
                                                                           computedInterfaceYThermalImpulse, computedInterfaceZThermalImpulse, material1Parameters);

    double computedMaterial2TotalEnergy = conservedVariableVector[20] / (computedMaterial2VolumeFraction * computedMaterial2Density);
    double computedMaterial2Pressure = HPREquationOfState::computePressure(computedMaterial2Density, computedMaterial2TotalEnergy, computedInterfaceXVelocity, computedInterfaceYVelocity,
                                                                           computedInterfaceZVelocity, computedInterfaceDistortionTensor, computedInterfaceXThermalImpulse,
                                                                           computedInterfaceYThermalImpulse, computedInterfaceZThermalImpulse, material2Parameters);

    double computedInterfacePressure = (computedMaterial1VolumeFraction * computedMaterial1Pressure) + (computedMaterial2VolumeFraction * computedMaterial2Pressure);

    vector<vector<double> > computedMaterial1ShearStressTensor = HPREquationOfState::computeShearStressTensor(computedMaterial1Density, computedInterfaceDistortionTensor, material1Parameters);
    vector<vector<double> > computedMaterial2ShearStressTensor = HPREquationOfState::computeShearStressTensor(computedMaterial2Density, computedInterfaceDistortionTensor, material2Parameters);
    vector<vector<double> > computedInterfaceShearStressTensor = MatrixAlgebra::addMatrices(MatrixAlgebra::multiplyMatrix(computedMaterial1VolumeFraction, computedMaterial1ShearStressTensor),
                                                                                            MatrixAlgebra::multiplyMatrix(computedMaterial2VolumeFraction, computedMaterial2ShearStressTensor));

    for (int i = 0; i < 21; i++)
    {
        fluxVector[i] = 0.0;
    }

    fluxVector[0] = computedTotalDensity * computedInterfaceYVelocity;
    fluxVector[1] = computedTotalDensity * (computedInterfaceYVelocity * computedMaterial1VolumeFraction);

    fluxVector[18] = computedMaterial1VolumeFraction * ((computedMaterial1Density * (computedInterfaceYVelocity * computedMaterial1TotalEnergy)) +
                                                        (computedMaterial1Pressure * computedInterfaceYVelocity));
    fluxVector[20] = computedMaterial2VolumeFraction * ((computedMaterial2Density * (computedInterfaceYVelocity * computedMaterial2TotalEnergy)) +
                                                        (computedMaterial2Pressure * computedInterfaceYVelocity));

    fluxVector[11] = computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceXVelocity);
    fluxVector[12] = (computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceYVelocity)) + computedInterfacePressure;
    fluxVector[13] = computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceZVelocity);

    vector<double> computedVelocityVector(3);
    computedVelocityVector[0] = computedInterfaceXVelocity;
    computedVelocityVector[1] = computedInterfaceYVelocity;
    computedVelocityVector[2] = computedInterfaceZVelocity;
    vector<double> computedDistortionTensorVelocityVectorProduct = MatrixAlgebra::multiplyMatrixByVector(computedInterfaceDistortionTensor, computedVelocityVector);

    fluxVector[18] -= computedMaterial1VolumeFraction * VectorAlgebra::computeDotProduct(computedMaterial1ShearStressTensor[1], computedVelocityVector);
    fluxVector[20] -= computedMaterial2VolumeFraction * VectorAlgebra::computeDotProduct(computedMaterial2ShearStressTensor[1], computedVelocityVector);

    for (int i = 0; i < 3; i++)
    {
        fluxVector[11 + i] -= computedInterfaceShearStressTensor[1][i];
    }

    for (int i = 0; i < 3; i++)
    {
        fluxVector[2 + (i * 3) + 1] = computedDistortionTensorVelocityVectorProduct[i];
    }

    if (material1Parameters.getIsThermal() || material2Parameters.getIsThermal())
    {
        double computedMaterial1Temperature = HPREquationOfState::computeTemperature(computedMaterial1Density, computedMaterial1Pressure, material1Parameters);
        vector<double> computedMaterial1HeatFluxVector = HPREquationOfState::computeHeatFluxVector(computedMaterial1Temperature, computedInterfaceXThermalImpulse,
                                                                                                   computedInterfaceYThermalImpulse, computedInterfaceZThermalImpulse, material1Parameters);

        double computedMaterial2Temperature = HPREquationOfState::computeTemperature(computedMaterial2Density, computedMaterial2Pressure, material2Parameters);
        vector<double> computedMaterial2HeatFluxVector = HPREquationOfState::computeHeatFluxVector(computedMaterial2Temperature, computedInterfaceXThermalImpulse,
                                                                                                   computedInterfaceYThermalImpulse, computedInterfaceZThermalImpulse, material2Parameters);

        double computedInterfaceTemperature = (computedMaterial1VolumeFraction * computedMaterial1Temperature) + (computedMaterial2VolumeFraction * computedMaterial2Temperature);

        fluxVector[18] += computedMaterial1VolumeFraction * computedMaterial1HeatFluxVector[1];
        fluxVector[20] += computedMaterial2VolumeFraction * computedMaterial2HeatFluxVector[1];

        fluxVector[14] = computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceXThermalImpulse);
        fluxVector[15] = (computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceYThermalImpulse)) + computedInterfaceTemperature;
        fluxVector[16] = computedTotalDensity * (computedInterfaceYVelocity * computedInterfaceZThermalImpulse);
    }

    fluxVector[17] = computedMaterial1VolumeFraction * (computedMaterial1Density * computedInterfaceYVelocity);
    fluxVector[19] = computedMaterial2VolumeFraction * (computedMaterial2Density * computedInterfaceYVelocity);

    return fluxVector;
}

vector<double> HPRIntermediateStateVector::computeYFluxVector(HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    return computeYFluxVector(computeConservedVariableVector(material1Parameters, material2Parameters), material1Parameters, material2Parameters);
}

vector<double> HPRIntermediateStateVector::computeSourceTermVector(vector<double> conservedVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> sourceTermVector(21);

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

    vector<vector<double> > computedInterfaceDistortionTensor(3, vector<double>(3));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            computedInterfaceDistortionTensor[i][j] = conservedVariableVector[2 + (i * 3) + j];
        }
    }

    double computedMaterial1Density = conservedVariableVector[17] / computedMaterial1VolumeFraction;
    double computedMaterial2Density = conservedVariableVector[19] / computedMaterial2VolumeFraction;

    vector<vector<double> > material1TotalEnergyDerivativeDistortionTensor = HPRDerivatives::computeTotalEnergyDerivativeDistortionTensor(computedInterfaceDistortionTensor, material1Parameters);
    vector<vector<double> > material2TotalEnergyDerivativeDistortionTensor = HPRDerivatives::computeTotalEnergyDerivativeDistortionTensor(computedInterfaceDistortionTensor, material2Parameters);
    vector<vector<double> > totalEnergyDerivativeDistortionTensor = MatrixAlgebra::addMatrices(MatrixAlgebra::multiplyMatrix(computedMaterial1VolumeFraction,
                                                                                                                             material1TotalEnergyDerivativeDistortionTensor),
                                                                                               MatrixAlgebra::multiplyMatrix(computedMaterial2VolumeFraction,
                                                                                                                             material2TotalEnergyDerivativeDistortionTensor));

    double material1Theta1Reciprocal = HPRSourceTerms::computeTheta1Reciprocal(computedMaterial1Density, computedInterfaceDistortionTensor, material1Parameters);
    double material2Theta1Reciprocal = HPRSourceTerms::computeTheta1Reciprocal(computedMaterial2Density, computedInterfaceDistortionTensor, material2Parameters);
    double theta1Reciprocal = (computedMaterial1VolumeFraction * material1Theta1Reciprocal) + (computedMaterial2VolumeFraction * material2Theta1Reciprocal);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sourceTermVector[2 + (i * 3) + j] = -totalEnergyDerivativeDistortionTensor[i][j] * theta1Reciprocal;
        }
    }

    if (material1Parameters.getIsThermal() || material2Parameters.getIsThermal())
    {
        double computedMaterial1TotalEnergy = conservedVariableVector[18] / (computedMaterial1VolumeFraction * computedMaterial1Density);
        double computedMaterial2TotalEnergy = conservedVariableVector[20] / (computedMaterial2VolumeFraction * computedMaterial2Density);

        double computedInterfaceXVelocity = conservedVariableVector[11] / computedTotalDensity;
        double computedInterfaceYVelocity = conservedVariableVector[12] / computedTotalDensity;
        double computedInterfaceZVelocity = conservedVariableVector[13] / computedTotalDensity;

        double computedInterfaceXThermalImpulse = conservedVariableVector[14] / computedTotalDensity;
        double computedInterfaceYThermalImpulse = conservedVariableVector[15] / computedTotalDensity;
        double computedInterfaceZThermalImpulse = conservedVariableVector[16] / computedTotalDensity;

        double computedMaterial1Pressure = HPREquationOfState::computePressure(computedMaterial1Density, computedMaterial1TotalEnergy, computedInterfaceXVelocity, computedInterfaceYVelocity,
                                                                               computedInterfaceZVelocity, computedInterfaceDistortionTensor, computedInterfaceXThermalImpulse,
                                                                               computedInterfaceYThermalImpulse, computedInterfaceZThermalImpulse, material1Parameters);
        double computedMaterial2Pressure = HPREquationOfState::computePressure(computedMaterial2Density, computedMaterial2TotalEnergy, computedInterfaceXVelocity, computedInterfaceYVelocity,
                                                                               computedInterfaceZVelocity, computedInterfaceDistortionTensor, computedInterfaceXThermalImpulse,
                                                                               computedInterfaceYThermalImpulse, computedInterfaceZThermalImpulse, material2Parameters);

        double computedMaterial1Temperature = HPREquationOfState::computeTemperature(computedMaterial1Density, computedMaterial1Pressure, material1Parameters);
        double computedMaterial2Temperature = HPREquationOfState::computeTemperature(computedMaterial2Density, computedMaterial2Pressure, material2Parameters);

        vector<double> material1TotalEnergyDerivativeThermalImpulse = HPRDerivatives::computeTotalEnergyDerivativeThermalImpulse(computedInterfaceXThermalImpulse, computedInterfaceYThermalImpulse,
                                                                                                                                 computedInterfaceZThermalImpulse, material1Parameters);
        vector<double> material2TotalEnergyDerivativeThermalImpulse = HPRDerivatives::computeTotalEnergyDerivativeThermalImpulse(computedInterfaceXThermalImpulse, computedInterfaceYThermalImpulse,
                                                                                                                                 computedInterfaceZThermalImpulse, material2Parameters);
        vector<double> totalEnergyDerivativeThermalImpulse = VectorAlgebra::addVectors(VectorAlgebra::multiplyVector(computedMaterial1VolumeFraction, material1TotalEnergyDerivativeThermalImpulse),
                                                                                       VectorAlgebra::multiplyVector(computedMaterial2VolumeFraction, material2TotalEnergyDerivativeThermalImpulse));

        double material1Theta2Reciprocal = HPRSourceTerms::computeTheta2Reciprocal(computedMaterial1Density, computedMaterial1Temperature, material1Parameters);
        double material2Theta2Reciprocal = HPRSourceTerms::computeTheta2Reciprocal(computedMaterial2Density, computedMaterial2Temperature, material2Parameters);
        double theta2Reciprocal = (computedMaterial1VolumeFraction * material1Theta2Reciprocal) + (computedMaterial2VolumeFraction * material2Theta2Reciprocal);

        for (int i = 0; i < 3; i++)
        {
            sourceTermVector[14 + i] -= computedTotalDensity * theta2Reciprocal * totalEnergyDerivativeThermalImpulse[i];
        }
    }

    return sourceTermVector;
}

vector<double> HPRIntermediateStateVector::computeSourceTermVector(HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    return computeSourceTermVector(computeConservedVariableVector(material1Parameters, material2Parameters), material1Parameters, material2Parameters);
}

double HPRIntermediateStateVector::computeMaterial1TotalEnergy(HPRMaterialParameters material1Parameters)
{
    return HPREquationOfState::computeTotalEnergy(material1Density, material1Pressure, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity, interfaceDistortionTensor,
                                                  interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse, material1Parameters);
}

double HPRIntermediateStateVector::computeMaterial1Temperature(HPRMaterialParameters material1Parameters)
{
    return HPREquationOfState::computeTemperature(material1Density, material1Pressure, material1Parameters);
}

vector<double> HPRIntermediateStateVector::computeMaterial1HeatFluxVector(HPRMaterialParameters material1Parameters)
{
    return HPREquationOfState::computeHeatFluxVector(computeMaterial1Temperature(material1Parameters), interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse,
                                                     material1Parameters);
}

double HPRIntermediateStateVector::computeMaterial1TotalEnergyDerivativeDensity(HPRMaterialParameters material1Parameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativeDensity(material1Density, material1Pressure, material1Parameters);
}

double HPRIntermediateStateVector::computeMaterial1TotalEnergyDerivativePressure(HPRMaterialParameters material1Parameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativePressure(material1Density, material1Parameters);
}

double HPRIntermediateStateVector::computeMaterial1TemperatureDerivativeDensity(HPRMaterialParameters material1Parameters)
{
    return HPRDerivatives::computeTemperatureDerivativeDensity(material1Density, material1Pressure, material1Parameters);
}

double HPRIntermediateStateVector::computeMaterial1TemperatureDerivativePressure(HPRMaterialParameters material1Parameters)
{
    return HPRDerivatives::computeTemperatureDerivativePressure(material1Density, material1Parameters);
}

double HPRIntermediateStateVector::computeMaterial1Theta1Reciprocal(HPRMaterialParameters material1Parameters)
{
    return HPRSourceTerms::computeTheta1Reciprocal(material1Density, interfaceDistortionTensor, material1Parameters);
}

double HPRIntermediateStateVector::computeMaterial1Theta2Reciprocal(HPRMaterialParameters material1Parameters)
{
    return HPRSourceTerms::computeTheta2Reciprocal(material1Density, computeMaterial1Temperature(material1Parameters), material1Parameters);
}

vector<vector<double> > HPRIntermediateStateVector::computeMaterial1ShearStressTensor(HPRMaterialParameters material1Parameters)
{
    return HPREquationOfState::computeShearStressTensor(material1Density, interfaceDistortionTensor, material1Parameters);
}

vector<vector<vector<vector<double> > > >  HPRIntermediateStateVector::computeMaterial1ShearStressTensorDerivativeDistortionTensor(HPRMaterialParameters material1Parameters)
{
    return HPREquationOfState::computeShearStressTensorDerivativeDistortionTensor(material1Density, interfaceDistortionTensor, material1Parameters);
}

double HPRIntermediateStateVector::computeMaterial2TotalEnergy(HPRMaterialParameters material2Parameters)
{
    return HPREquationOfState::computeTotalEnergy(material2Density, material2Pressure, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity, interfaceDistortionTensor,
                                                  interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse, material2Parameters);
}

double HPRIntermediateStateVector::computeMaterial2Temperature(HPRMaterialParameters material2Parameters)
{
    return HPREquationOfState::computeTemperature(material2Density, material2Pressure, material2Parameters);
}

vector<double> HPRIntermediateStateVector::computeMaterial2HeatFluxVector(HPRMaterialParameters material2Parameters)
{
    return HPREquationOfState::computeHeatFluxVector(computeMaterial2Temperature(material2Parameters), interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse,
                                                     material2Parameters);
}

double HPRIntermediateStateVector::computeMaterial2TotalEnergyDerivativeDensity(HPRMaterialParameters material2Parameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativeDensity(material2Density, material2Pressure, material2Parameters);
}

double HPRIntermediateStateVector::computeMaterial2TotalEnergyDerivativePressure(HPRMaterialParameters material2Parameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativePressure(material2Density, material2Parameters);
}

double HPRIntermediateStateVector::computeMaterial2TemperatureDerivativeDensity(HPRMaterialParameters material2Parameters)
{
    return HPRDerivatives::computeTemperatureDerivativeDensity(material2Density, material2Pressure, material2Parameters);
}

double HPRIntermediateStateVector::computeMaterial2TemperatureDerivativePressure(HPRMaterialParameters material2Parameters)
{
    return HPRDerivatives::computeTemperatureDerivativePressure(material2Density, material2Parameters);
}

double HPRIntermediateStateVector::computeMaterial2Theta1Reciprocal(HPRMaterialParameters material2Parameters)
{
    return HPRSourceTerms::computeTheta1Reciprocal(material2Density, interfaceDistortionTensor, material2Parameters);
}

double HPRIntermediateStateVector::computeMaterial2Theta2Reciprocal(HPRMaterialParameters material2Parameters)
{
    return HPRSourceTerms::computeTheta2Reciprocal(material2Density, computeMaterial2Temperature(material2Parameters), material2Parameters);
}

vector<vector<double> > HPRIntermediateStateVector::computeMaterial2ShearStressTensor(HPRMaterialParameters material2Parameters)
{
    return HPREquationOfState::computeShearStressTensor(material2Density, interfaceDistortionTensor, material2Parameters);
}

vector<vector<vector<vector<double> > > > HPRIntermediateStateVector::computeMaterial2ShearStressTensorDerivativeDistortionTensor(HPRMaterialParameters material2Parameters)
{
    return HPREquationOfState::computeShearStressTensorDerivativeDistortionTensor(material2Density, interfaceDistortionTensor, material2Parameters);
}

vector<vector<double> > HPRIntermediateStateVector::computeTotalEnergyDerivativeDistortionTensor(HPRMaterialParameters materialParameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativeDistortionTensor(interfaceDistortionTensor, materialParameters);
}

vector<double> HPRIntermediateStateVector::computeTotalEnergyDerivativeThermalImpulse(HPRMaterialParameters materialParameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativeThermalImpulse(interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse, materialParameters);
}

vector<vector<double> > HPRIntermediateStateVector::computeShearStressTensorDerivativeDensity(HPRMaterialParameters materialParameters)
{
    return HPREquationOfState::computeShearStressTensorDerivativeDensity(interfaceDistortionTensor, materialParameters);
}

double HPRIntermediateStateVector::computeTotalDensity()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1Density) + (material2VolumeFraction * material2Density);
}

double HPRIntermediateStateVector::computeTotalPressure()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1Pressure) + (material2VolumeFraction * material2Pressure);
}

void HPRIntermediateStateVector::relaxTotalDensity()
{
    double totalDensity = computeTotalDensity();

    material1Density = totalDensity;
    material2Density = totalDensity;
}

void HPRIntermediateStateVector::relaxTotalPressure()
{
    double totalPressure = computeTotalPressure();

    material1Pressure = totalPressure;
    material2Pressure = totalPressure;
}

void HPRIntermediateStateVector::setMaterial1VolumeFraction(double newMaterial1VolumeFraction)
{
    material1VolumeFraction = newMaterial1VolumeFraction;
}

void HPRIntermediateStateVector::setInterfaceDistortionTensor(vector<vector<double> > newInterfaceDistortionTensor)
{
    interfaceDistortionTensor = newInterfaceDistortionTensor;
}

void HPRIntermediateStateVector::setInterfaceXVelocity(double newInterfaceXVelocity)
{
    interfaceXVelocity = newInterfaceXVelocity;
}

void HPRIntermediateStateVector::setInterfaceYVelocity(double newInterfaceYVelocity)
{
    interfaceYVelocity = newInterfaceYVelocity;
}

void HPRIntermediateStateVector::setInterfaceZVelocity(double newInterfaceZVelocity)
{
    interfaceZVelocity = newInterfaceZVelocity;
}

void HPRIntermediateStateVector::setInterfaceXThermalImpulse(double newInterfaceXThermalImpulse)
{
    interfaceXThermalImpulse = newInterfaceXThermalImpulse;
}

void HPRIntermediateStateVector::setInterfaceYThermalImpulse(double newInterfaceYThermalImpulse)
{
    interfaceYThermalImpulse = newInterfaceYThermalImpulse;
}

void HPRIntermediateStateVector::setInterfaceZThermalImpulse(double newInterfaceZThermalImpulse)
{
    interfaceZThermalImpulse = newInterfaceZThermalImpulse;
}

void HPRIntermediateStateVector::setMaterial1Density(double newMaterial1Density)
{
    material1Density = newMaterial1Density;
}

void HPRIntermediateStateVector::setMaterial1Pressure(double newMaterial1Pressure)
{
    material1Pressure = newMaterial1Pressure;
}

void HPRIntermediateStateVector::setMaterial2Density(double newMaterial2Density)
{
    material2Density = newMaterial2Density;
}

void HPRIntermediateStateVector::setMaterial2Pressure(double newMaterial2Pressure)
{
    material2Pressure = newMaterial2Pressure;
}

double HPRIntermediateStateVector::getMaterial1VolumeFraction()
{
    return material1VolumeFraction;
}

vector<vector<double> > HPRIntermediateStateVector::getInterfaceDistortionTensor()
{
    return interfaceDistortionTensor;
}

double HPRIntermediateStateVector::getInterfaceXVelocity()
{
    return interfaceXVelocity;
}

double HPRIntermediateStateVector::getInterfaceYVelocity()
{
    return interfaceYVelocity;
}

double HPRIntermediateStateVector::getInterfaceZVelocity()
{
    return interfaceZVelocity;
}

double HPRIntermediateStateVector::getInterfaceXThermalImpulse()
{
    return interfaceXThermalImpulse;
}

double HPRIntermediateStateVector::getInterfaceYThermalImpulse()
{
    return interfaceYThermalImpulse;
}

double HPRIntermediateStateVector::getInterfaceZThermalImpulse()
{
    return interfaceZThermalImpulse;
}

double HPRIntermediateStateVector::getMaterial1Density()
{
    return material1Density;
}

double HPRIntermediateStateVector::getMaterial1Pressure()
{
    return material1Pressure;
}

double HPRIntermediateStateVector::getMaterial2Density()
{
    return material2Density;
}

double HPRIntermediateStateVector::getMaterial2Pressure()
{
    return material2Pressure;
}

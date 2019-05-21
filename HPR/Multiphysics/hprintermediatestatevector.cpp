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

    //vector<vector<double> > computedMaterial1ShearStressTensor = HPREquationOfState::computeShearStressTensor(computedMaterial1Density, computedInterfaceDistortionTensor, material1Parameters);
    //vector<vector<double> > computedMaterial2ShearStressTensor = HPREquationOfState::computeShearStressTensor(computedMaterial2Density, computedInterfaceDistortionTensor, material2Parameters);

    return fluxVector;
}

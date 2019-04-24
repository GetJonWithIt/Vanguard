#include "hprreducedstatevector.h"

HPRReducedStateVector::HPRReducedStateVector()
{
    material1VolumeFraction = 0.999;
    interfacePressure = 1.0 / 1.4;
    interfaceDistortionTensor = MatrixAlgebra::computeIdentityMatrix(3);

    interfaceXVelocity = 0.0;
    interfaceYVelocity = 0.0;
    interfaceZVelocity = 0.0;

    interfaceXThermalImpulse = 0.0;
    interfaceYThermalImpulse = 0.0;
    interfaceZThermalImpulse = 0.0;

    material1Density = 1.0;
    material2Density = 1.0;
}

HPRReducedStateVector::HPRReducedStateVector(double newMaterial1VolumeFraction, double newInterfacePressure, vector<vector<double> > newInterfaceDistortionTensor, double newInterfaceXVelocity,
                                             double newInterfaceYVelocity, double newInterfaceZVelocity, double newInterfaceXThermalImpulse, double newInterfaceYThermalImpulse,
                                             double newInterfaceZThermalImpulse, double newMaterial1Density, double newMaterial2Density)
{
    material1VolumeFraction = newMaterial1VolumeFraction;
    interfacePressure = newInterfacePressure;
    interfaceDistortionTensor = newInterfaceDistortionTensor;

    interfaceXVelocity = newInterfaceXVelocity;
    interfaceYVelocity = newInterfaceYVelocity;
    interfaceZVelocity = newInterfaceZVelocity;

    interfaceXThermalImpulse = newInterfaceXThermalImpulse;
    interfaceYThermalImpulse = newInterfaceYThermalImpulse;
    interfaceZThermalImpulse = newInterfaceZThermalImpulse;

    material1Density = newMaterial1Density;
    material2Density = newMaterial2Density;
}

HPRReducedStateVector::HPRReducedStateVector(double newMaterial1VolumeFraction, double newInterfaceXVelocity, double newInterfaceYVelocity, double newInterfaceZVelocity,
                                             vector<vector<double> > newMaterial1DistortionTensor, double newMaterial1Entropy, vector<vector<double> > newMaterial2DistortionTensor,
                                             double newMaterial2Entropy, HyperelasticMaterialParameters hyperelasticMaterial1Parameters,
                                             HyperelasticMaterialParameters hyperelasticMaterial2Parameters, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
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

    double material1Pressure = HPREquationOfState::computePressure(material1Density, material1TotalEnergy, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                                   newMaterial1DistortionTensor, interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse, material1Parameters);
    double material2Pressure = HPREquationOfState::computePressure(material2Density, material2TotalEnergy, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                                   newMaterial2DistortionTensor, interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse, material2Parameters);

    interfacePressure = (material1VolumeFraction * material1Pressure) + (material2VolumeFraction * material2Pressure);
}

void HPRReducedStateVector::setPrimitiveVariableVector(vector<double> newPrimitiveVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    material1VolumeFraction = newPrimitiveVariableVector[0];
    interfacePressure = newPrimitiveVariableVector[1];

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            interfaceDistortionTensor[i][j] = newPrimitiveVariableVector[2 + (i * 3) + j];
        }
    }

    interfaceXVelocity = newPrimitiveVariableVector[11];
    interfaceYVelocity = newPrimitiveVariableVector[12];
    interfaceZVelocity = newPrimitiveVariableVector[13];

    if (material1Parameters.getIsThermal() || material2Parameters.getIsThermal())
    {
        interfaceXThermalImpulse = newPrimitiveVariableVector[14];
        interfaceYThermalImpulse = newPrimitiveVariableVector[15];
        interfaceZThermalImpulse = newPrimitiveVariableVector[16];
    }

    material1Density = newPrimitiveVariableVector[17];
    material2Density = newPrimitiveVariableVector[18];
}

void HPRReducedStateVector::setConservedVariableVector(vector<double> newConservedVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
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
            interfaceDistortionTensor[i][j] = newConservedVariableVector[3 + (i * 3) + j];
        }
    }

    interfaceXVelocity = newConservedVariableVector[12] / totalDensity;
    interfaceYVelocity = newConservedVariableVector[13] / totalDensity;
    interfaceZVelocity = newConservedVariableVector[14] / totalDensity;

    if (material1Parameters.getIsThermal() || material2Parameters.getIsThermal())
    {
        interfaceXThermalImpulse = newConservedVariableVector[15] / totalDensity;
        interfaceYThermalImpulse = newConservedVariableVector[16] / totalDensity;
        interfaceZThermalImpulse = newConservedVariableVector[17] / totalDensity;
    }

    material1Density = newConservedVariableVector[18] / material1VolumeFraction;
    material2Density = newConservedVariableVector[19] / material2VolumeFraction;

    double computedTotalEnergy = newConservedVariableVector[2] / totalDensity;
    double computedMaterial1Pressure = HPREquationOfState::computePressure(material1Density, computedTotalEnergy, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                                           interfaceDistortionTensor, interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse,
                                                                           material1Parameters);
    double computedMaterial2Pressure = HPREquationOfState::computePressure(material2Density, computedTotalEnergy, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                                           interfaceDistortionTensor, interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse,
                                                                           material2Parameters);
    interfacePressure = (material1VolumeFraction * computedMaterial1Pressure) + (material2VolumeFraction * computedMaterial2Pressure);
}

vector<double> HPRReducedStateVector::computePrimitiveVariableVector(HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> primitiveVariableVector(19);

    primitiveVariableVector[0] = material1VolumeFraction;
    primitiveVariableVector[1] = interfacePressure;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            primitiveVariableVector[2 + (i * 3) + j] = interfaceDistortionTensor[i][j];
        }
    }

    primitiveVariableVector[11] = interfaceXVelocity;
    primitiveVariableVector[12] = interfaceYVelocity;
    primitiveVariableVector[13] = interfaceZVelocity;

    if (material1Parameters.getIsThermal() || material2Parameters.getIsThermal())
    {
        primitiveVariableVector[14] = interfaceXThermalImpulse;
        primitiveVariableVector[15] = interfaceYThermalImpulse;
        primitiveVariableVector[16] = interfaceZThermalImpulse;
    }

    primitiveVariableVector[17] = material1Density;
    primitiveVariableVector[18] = material2Density;

    return primitiveVariableVector;
}

vector<double> HPRReducedStateVector::computeConservedVariableVector(HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> conservedVariableVector(20);

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

    double computedMaterial1TotalEnergy = HPREquationOfState::computeTotalEnergy(material1Density, interfacePressure, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                                                 interfaceDistortionTensor, interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse,
                                                                                 material1Parameters);
    double computedMaterial2TotalEnergy = HPREquationOfState::computeTotalEnergy(material2Density, interfacePressure, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity,
                                                                                 interfaceDistortionTensor, interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse,
                                                                                 material2Parameters);
    double computedTotalEnergy = (material1VolumeFraction * computedMaterial1TotalEnergy) + (material2VolumeFraction * computedMaterial2TotalEnergy);
    conservedVariableVector[2] = totalDensity * computedTotalEnergy;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            conservedVariableVector[3 + (i * 3) + j] = interfaceDistortionTensor[i][j];
        }
    }

    conservedVariableVector[12] = totalDensity * interfaceXVelocity;
    conservedVariableVector[13] = totalDensity * interfaceYVelocity;
    conservedVariableVector[14] = totalDensity * interfaceZVelocity;

    if (material1Parameters.getIsThermal() || material2Parameters.getIsThermal())
    {
        conservedVariableVector[15] = totalDensity * interfaceXThermalImpulse;
        conservedVariableVector[16] = totalDensity * interfaceYThermalImpulse;
        conservedVariableVector[17] = totalDensity * interfaceZThermalImpulse;
    }

    conservedVariableVector[18] = material1VolumeFraction * material1Density;
    conservedVariableVector[19] = material2VolumeFraction * material2Density;

    return conservedVariableVector;
}

vector<double> HPRReducedStateVector::computeXFluxVector(vector<double> conservedVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> fluxVector(20);

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
            computedInterfaceDistortionTensor[i][j] = conservedVariableVector[3 + (i * 3) + j];
        }
    }

    double computedInterfaceXVelocity = conservedVariableVector[12] / computedTotalDensity;
    double computedInterfaceYVelocity = conservedVariableVector[13] / computedTotalDensity;
    double computedInterfaceZVelocity = conservedVariableVector[14] / computedTotalDensity;

    double computedInterfaceXThermalImpulse = conservedVariableVector[15] / computedTotalDensity;
    double computedInterfaceYThermalImpulse = conservedVariableVector[16] / computedTotalDensity;
    double computedInterfaceZThermalImpulse = conservedVariableVector[17] / computedTotalDensity;

    double computedMaterial1Density = conservedVariableVector[18] / computedMaterial1VolumeFraction;
    double computedMaterial2Density = conservedVariableVector[19] / computedMaterial2VolumeFraction;

    double computedTotalEnergy = conservedVariableVector[2] / computedTotalDensity;
    double computedMaterial1Pressure = HPREquationOfState::computePressure(computedMaterial1Density, computedTotalEnergy, computedInterfaceXVelocity, computedInterfaceYVelocity,
                                                                           computedInterfaceZVelocity, computedInterfaceDistortionTensor, computedInterfaceXThermalImpulse,
                                                                           computedInterfaceYThermalImpulse, computedInterfaceZThermalImpulse, material1Parameters);
    double computedMaterial2Pressure = HPREquationOfState::computePressure(computedMaterial2Density, computedTotalEnergy, computedInterfaceXVelocity, computedInterfaceYVelocity,
                                                                           computedInterfaceZVelocity, computedInterfaceDistortionTensor, computedInterfaceXThermalImpulse,
                                                                           computedInterfaceYThermalImpulse, computedInterfaceZThermalImpulse, material2Parameters);
    double computedInterfacePressure = (computedMaterial1VolumeFraction * computedMaterial1Pressure) + (computedMaterial2VolumeFraction * computedMaterial2Pressure);

    vector<vector<double> > computedMaterial1ShearStressTensor = HPREquationOfState::computeShearStressTensor(computedMaterial1Density, computedInterfaceDistortionTensor, material1Parameters);
    vector<vector<double> > computedMaterial2ShearStressTensor = HPREquationOfState::computeShearStressTensor(computedMaterial2Density, computedInterfaceDistortionTensor, material2Parameters);
    vector<vector<double> > computedInterfaceShearStressTensor = MatrixAlgebra::addMatrices(MatrixAlgebra::multiplyMatrix(computedMaterial1VolumeFraction, computedMaterial1ShearStressTensor),
                                                                                            MatrixAlgebra::multiplyMatrix(computedMaterial2VolumeFraction, computedMaterial2ShearStressTensor));

    for (int i = 0; i < 20; i++)
    {
        fluxVector[i] = 0.0;
    }

    fluxVector[0] = computedTotalDensity * computedInterfaceXVelocity;
    fluxVector[1] = computedTotalDensity * (computedInterfaceXVelocity * computedMaterial1VolumeFraction);
    fluxVector[2] = (computedTotalDensity * (computedInterfaceXVelocity * computedTotalEnergy)) + (computedInterfacePressure * computedInterfaceXVelocity);

    fluxVector[12] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceXVelocity)) + computedInterfacePressure;
    fluxVector[13] = computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceYVelocity);
    fluxVector[14] = computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceZVelocity);

    vector<double> computedVelocityVector(3);
    computedVelocityVector[0] = computedInterfaceXVelocity;
    computedVelocityVector[1] = computedInterfaceYVelocity;
    computedVelocityVector[2] = computedInterfaceZVelocity;
    vector<double> distortionTensorVelocityVectorProduct = MatrixAlgebra::multiplyMatrixByVector(computedInterfaceDistortionTensor, computedVelocityVector);

    fluxVector[2] -= VectorAlgebra::computeDotProduct(computedInterfaceShearStressTensor[0], computedVelocityVector);
    for (int i = 0; i < 3; i++)
    {
        fluxVector[12 + i] -= computedInterfaceShearStressTensor[0][i];
    }

    for (int i = 0; i < 3; i++)
    {
        fluxVector[3 + (i * 3)] = distortionTensorVelocityVectorProduct[i];
    }

    if (material1Parameters.getIsThermal() || material2Parameters.getIsThermal())
    {
        double computedMaterial1Temperature = HPREquationOfState::computeTemperature(computedMaterial1Density, computedInterfacePressure, material1Parameters);
        vector<double> computedMaterial1HeatFluxVector = HPREquationOfState::computeHeatFluxVector(computedMaterial1Temperature, computedInterfaceXThermalImpulse,
                                                                                                   computedInterfaceYThermalImpulse, computedInterfaceZThermalImpulse, material1Parameters);

        double computedMaterial2Temperature = HPREquationOfState::computeTemperature(computedMaterial2Density, computedInterfacePressure, material2Parameters);
        vector<double> computedMaterial2HeatFluxVector = HPREquationOfState::computeHeatFluxVector(computedMaterial2Temperature, computedInterfaceXThermalImpulse,
                                                                                                   computedInterfaceYThermalImpulse, computedInterfaceZThermalImpulse, material2Parameters);

        double computedInterfaceTemperature = (computedMaterial1VolumeFraction * computedMaterial1Temperature) + (computedMaterial1VolumeFraction * computedMaterial2Temperature);
        vector<double> computedInterfaceHeatFluxVector = VectorAlgebra::addVectors(VectorAlgebra::multiplyVector(computedMaterial1VolumeFraction, computedMaterial1HeatFluxVector),
                                                                                   VectorAlgebra::multiplyVector(computedMaterial2VolumeFraction, computedMaterial2HeatFluxVector));

        fluxVector[2] += computedInterfaceHeatFluxVector[0];

        fluxVector[15] = (computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceXThermalImpulse)) + computedInterfaceTemperature;
        fluxVector[16] = computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceYThermalImpulse);
        fluxVector[17] = computedTotalDensity * (computedInterfaceXVelocity * computedInterfaceZThermalImpulse);
    }

    fluxVector[18] = computedMaterial1VolumeFraction * (computedMaterial1Density * computedInterfaceXVelocity);
    fluxVector[19] = computedMaterial2VolumeFraction * (computedMaterial2Density * computedInterfaceXVelocity);

    return fluxVector;
}

vector<double> HPRReducedStateVector::computeXFluxVector(HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    return computeXFluxVector(computeConservedVariableVector(material1Parameters, material2Parameters), material1Parameters, material2Parameters);
}

vector<double> HPRReducedStateVector::computeSourceTermVector(vector<double> conservedVariableVector, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> sourceTermVector(20);

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
            computedInterfaceDistortionTensor[i][j] = conservedVariableVector[3 + (i * 3) + j];
        }
    }

    double computedMaterial1Density = conservedVariableVector[18] / computedMaterial1VolumeFraction;
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
            sourceTermVector[3 + (i * 3) + j] = -totalEnergyDerivativeDistortionTensor[i][j] * theta1Reciprocal;
        }
    }

    if (material1Parameters.getIsThermal() || material2Parameters.getIsThermal())
    {
        double computedTotalEnergy = conservedVariableVector[2] / computedTotalDensity;

        double computedInterfaceXVelocity = conservedVariableVector[12] / computedTotalDensity;
        double computedInterfaceYVelocity = conservedVariableVector[13] / computedTotalDensity;
        double computedInterfaceZVelocity = conservedVariableVector[14] / computedTotalDensity;

        double computedInterfaceXThermalImpulse = conservedVariableVector[15] / computedTotalDensity;
        double computedInterfaceYThermalImpulse = conservedVariableVector[16] / computedTotalDensity;
        double computedInterfaceZThermalImpulse = conservedVariableVector[17] / computedTotalDensity;

        double computedMaterial1Pressure = HPREquationOfState::computePressure(computedMaterial1Density, computedTotalEnergy, computedInterfaceXVelocity, computedInterfaceYVelocity,
                                                                               computedInterfaceZVelocity, computedInterfaceDistortionTensor, computedInterfaceXThermalImpulse,
                                                                               computedInterfaceYThermalImpulse, computedInterfaceZThermalImpulse, material1Parameters);
        double computedMaterial2Pressure = HPREquationOfState::computePressure(computedMaterial2Density, computedTotalEnergy, computedInterfaceXVelocity, computedInterfaceYVelocity,
                                                                               computedInterfaceZVelocity, computedInterfaceDistortionTensor, computedInterfaceXThermalImpulse,
                                                                               computedInterfaceYThermalImpulse, computedInterfaceZThermalImpulse, material2Parameters);
        double computedInterfacePressure = (computedMaterial1VolumeFraction * computedMaterial1Pressure) + (computedMaterial2VolumeFraction * computedMaterial2Pressure);

        double computedMaterial1Temperature = HPREquationOfState::computeTemperature(computedMaterial1Density, computedInterfacePressure, material1Parameters);
        double computedMaterial2Temperature = HPREquationOfState::computeTemperature(computedMaterial2Density, computedInterfacePressure, material2Parameters);
        double computedInterfaceTemperature = (computedMaterial1VolumeFraction * computedMaterial1Temperature) + (computedMaterial2VolumeFraction * computedMaterial2Temperature);

        vector<double> material1TotalEnergyDerivativeThermalImpulse = HPRDerivatives::computeTotalEnergyDerivativeThermalImpulse(computedInterfaceXThermalImpulse, computedInterfaceYThermalImpulse,
                                                                                                                                 computedInterfaceZThermalImpulse, material1Parameters);
        vector<double> material2TotalEnergyDerivativeThermalImpulse = HPRDerivatives::computeTotalEnergyDerivativeThermalImpulse(computedInterfaceXThermalImpulse, computedInterfaceYThermalImpulse,
                                                                                                                                 computedInterfaceZThermalImpulse, material2Parameters);
        vector<double> totalEnergyDerivativeThermalImpulse = VectorAlgebra::addVectors(VectorAlgebra::multiplyVector(computedMaterial1VolumeFraction, material1TotalEnergyDerivativeThermalImpulse),
                                                                                       VectorAlgebra::multiplyVector(computedMaterial2VolumeFraction, material2TotalEnergyDerivativeThermalImpulse));

        double material1Theta2Reciprocal = HPRSourceTerms::computeTheta2Reciprocal(computedMaterial1Density, computedInterfaceTemperature, material1Parameters);
        double material2Theta2Reciprocal = HPRSourceTerms::computeTheta2Reciprocal(computedMaterial2Density, computedInterfaceTemperature, material2Parameters);
        double theta2Reciprocal = (computedMaterial1VolumeFraction * material1Theta2Reciprocal) + (computedMaterial2VolumeFraction * material2Theta2Reciprocal);

        for (int i = 0; i < 3; i++)
        {
            sourceTermVector[15 + i] = -computedTotalDensity * theta2Reciprocal * totalEnergyDerivativeThermalImpulse[i];
        }
    }

    return sourceTermVector;
}

double HPRReducedStateVector::computeMaterial1TotalEnergy(HPRMaterialParameters material1Parameters)
{
    return HPREquationOfState::computeTotalEnergy(material1Density, interfacePressure, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity, interfaceDistortionTensor,
                                                  interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse, material1Parameters);
}

double HPRReducedStateVector::computeMaterial1Temperature(HPRMaterialParameters material1Parameters)
{
    return HPREquationOfState::computeTemperature(material1Density, interfacePressure, material1Parameters);
}

vector<double> HPRReducedStateVector::computeMaterial1HeatFluxVector(HPRMaterialParameters material1Parameters)
{
    return HPREquationOfState::computeHeatFluxVector(computeMaterial1Temperature(material1Parameters), interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse,
                                                     material1Parameters);
}

double HPRReducedStateVector::computeMaterial1TotalEnergyDerivativeDensity(HPRMaterialParameters material1Parameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativeDensity(material1Density, interfacePressure, material1Parameters);
}

double HPRReducedStateVector::computeMaterial1TotalEnergyDerivativePressure(HPRMaterialParameters material1Parameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativePressure(material1Density, material1Parameters);
}

double HPRReducedStateVector::computeMaterial1TemperatureDerivativeDensity(HPRMaterialParameters material1Parameters)
{
    return HPRDerivatives::computeTemperatureDerivativeDensity(material1Density, interfacePressure, material1Parameters);
}

double HPRReducedStateVector::computeMaterial1TemperatureDerivativePressure(HPRMaterialParameters material1Parameters)
{
    return HPRDerivatives::computeTemperatureDerivativePressure(material1Density, material1Parameters);
}

double HPRReducedStateVector::computeMaterial1Theta1Reciprocal(HPRMaterialParameters material1Parameters)
{
    return HPRSourceTerms::computeTheta1Reciprocal(material1Density, interfaceDistortionTensor, material1Parameters);
}

double HPRReducedStateVector::computeMaterial1Theta2Reciprocal(HPRMaterialParameters material1Parameters)
{
    return HPRSourceTerms::computeTheta2Reciprocal(material1Density, computeMaterial1Temperature(material1Parameters), material1Parameters);
}

vector<vector<double> > HPRReducedStateVector::computeMaterial1ShearStressTensor(HPRMaterialParameters material1Parameters)
{
    return HPREquationOfState::computeShearStressTensor(material1Density, interfaceDistortionTensor, material1Parameters);
}

vector<vector<vector<vector<double> > > > HPRReducedStateVector::computeMaterial1ShearStressTensorDerivativeDistortionTensor(HPRMaterialParameters material1Parameters)
{
    return HPREquationOfState::computeShearStressTensorDerivativeDistortionTensor(material1Density, interfaceDistortionTensor, material1Parameters);
}

double HPRReducedStateVector::computeMaterial2TotalEnergy(HPRMaterialParameters material2Parameters)
{
    return HPREquationOfState::computeTotalEnergy(material2Density, interfacePressure, interfaceXVelocity, interfaceYVelocity, interfaceZVelocity, interfaceDistortionTensor,
                                                  interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse, material2Parameters);
}

double HPRReducedStateVector::computeMaterial2Temperature(HPRMaterialParameters material2Parameters)
{
    return HPREquationOfState::computeTemperature(material2Density, interfacePressure, material2Parameters);
}

vector<double> HPRReducedStateVector::computeMaterial2HeatFluxVector(HPRMaterialParameters material2Parameters)
{
    return HPREquationOfState::computeHeatFluxVector(computeMaterial2Temperature(material2Parameters), interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse,
                                                     material2Parameters);
}

double HPRReducedStateVector::computeMaterial2TotalEnergyDerivativeDensity(HPRMaterialParameters material2Parameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativeDensity(material2Density, interfacePressure, material2Parameters);
}

double HPRReducedStateVector::computeMaterial2TotalEnergyDerivativePressure(HPRMaterialParameters material2Parameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativePressure(material2Density, material2Parameters);
}

double HPRReducedStateVector::computeMaterial2TemperatureDerivativeDensity(HPRMaterialParameters material2Parameters)
{
    return HPRDerivatives::computeTemperatureDerivativeDensity(material2Density, interfacePressure, material2Parameters);
}

double HPRReducedStateVector::computeMaterial2TemperatureDerivativePressure(HPRMaterialParameters material2Parameters)
{
    return HPRDerivatives::computeTemperatureDerivativePressure(material2Density, material2Parameters);
}

double HPRReducedStateVector::computeMaterial2Theta1Reciprocal(HPRMaterialParameters material2Parameters)
{
    return HPRSourceTerms::computeTheta1Reciprocal(material2Density, interfaceDistortionTensor, material2Parameters);
}

double HPRReducedStateVector::computeMaterial2Theta2Reciprocal(HPRMaterialParameters material2Parameters)
{
    return HPRSourceTerms::computeTheta2Reciprocal(material2Density, computeMaterial2Temperature(material2Parameters), material2Parameters);
}

vector<vector<double> > HPRReducedStateVector::computeMaterial2ShearStressTensor(HPRMaterialParameters material2Parameters)
{
    return HPREquationOfState::computeShearStressTensor(material2Density, interfaceDistortionTensor, material2Parameters);
}

vector<vector<vector<vector<double> > > > HPRReducedStateVector::computeMaterial2ShearStressTensorDerivativeDistortionTensor(HPRMaterialParameters material2Parameters)
{
    return HPREquationOfState::computeShearStressTensorDerivativeDistortionTensor(material2Density, interfaceDistortionTensor, material2Parameters);
}

vector<vector<double> > HPRReducedStateVector::computeTotalEnergyDerivativeDistortionTensor(HPRMaterialParameters materialParameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativeDistortionTensor(interfaceDistortionTensor, materialParameters);
}

vector<double> HPRReducedStateVector::computeTotalEnergyDerivativeThermalImpulse(HPRMaterialParameters materialParameters)
{
    return HPRDerivatives::computeTotalEnergyDerivativeThermalImpulse(interfaceXThermalImpulse, interfaceYThermalImpulse, interfaceZThermalImpulse, materialParameters);
}

vector<vector<double> > HPRReducedStateVector::computeShearStressTensorDerivativeDensity(HPRMaterialParameters materialParameters)
{
    return HPREquationOfState::computeShearStressTensorDerivativeDensity(interfaceDistortionTensor, materialParameters);
}

double HPRReducedStateVector::computeTotalDensity()
{
    double material2VolumeFraction = 1.0 - material1VolumeFraction;

    return (material1VolumeFraction * material1Density) + (material2VolumeFraction * material2Density);
}

void HPRReducedStateVector::relaxTotalDensity()
{
    double totalDensity = computeTotalDensity();

    material1Density = totalDensity;
    material2Density = totalDensity;
}

void HPRReducedStateVector::setMaterial1VolumeFraction(double newMaterial1VolumeFraction)
{
    material1VolumeFraction = newMaterial1VolumeFraction;
}

void HPRReducedStateVector::setInterfacePressure(double newInterfacePressure)
{
    interfacePressure = newInterfacePressure;
}

void HPRReducedStateVector::setInterfaceDistortionTensor(vector<vector<double> > newInterfaceDistortionTensor)
{
    interfaceDistortionTensor = newInterfaceDistortionTensor;
}

void HPRReducedStateVector::setInterfaceXVelocity(double newInterfaceXVelocity)
{
    interfaceXVelocity = newInterfaceXVelocity;
}

void HPRReducedStateVector::setInterfaceYVelocity(double newInterfaceYVelocity)
{
    interfaceYVelocity = newInterfaceYVelocity;
}

void HPRReducedStateVector::setInterfaceZVelocity(double newInterfaceZVelocity)
{
    interfaceZVelocity = newInterfaceZVelocity;
}

void HPRReducedStateVector::setInterfaceXThermalImpulse(double newInterfaceXThermalImpulse)
{
    interfaceXThermalImpulse = newInterfaceXThermalImpulse;
}

void HPRReducedStateVector::setInterfaceYThermalImpulse(double newInterfaceYThermalImpulse)
{
    interfaceYThermalImpulse = newInterfaceYThermalImpulse;
}

void HPRReducedStateVector::setInterfaceZThermalImpulse(double newInterfaceZThermalImpulse)
{
    interfaceZThermalImpulse = newInterfaceZThermalImpulse;
}

void HPRReducedStateVector::setMaterial1Density(double newMaterial1Density)
{
    material1Density = newMaterial1Density;
}

void HPRReducedStateVector::setMaterial2Density(double newMaterial2Density)
{
    material2Density = newMaterial2Density;
}

double HPRReducedStateVector::getMaterial1VolumeFraction()
{
    return material1VolumeFraction;
}

double HPRReducedStateVector::getInterfacePressure()
{
    return interfacePressure;
}

vector<vector<double> > HPRReducedStateVector::getInterfaceDistortionTensor()
{
    return interfaceDistortionTensor;
}

double HPRReducedStateVector::getInterfaceXVelocity()
{
    return interfaceXVelocity;
}

double HPRReducedStateVector::getInterfaceYVelocity()
{
    return interfaceYVelocity;
}

double HPRReducedStateVector::getInterfaceZVelocity()
{
    return interfaceZVelocity;
}

double HPRReducedStateVector::getInterfaceXThermalImpulse()
{
    return interfaceXThermalImpulse;
}

double HPRReducedStateVector::getInterfaceYThermalImpulse()
{
    return interfaceYThermalImpulse;
}

double HPRReducedStateVector::getInterfaceZThermalImpulse()
{
    return interfaceZThermalImpulse;
}

double HPRReducedStateVector::getMaterial1Density()
{
    return material1Density;
}

double HPRReducedStateVector::getMaterial2Density()
{
    return material2Density;
}

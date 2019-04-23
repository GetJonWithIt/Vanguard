#include "hprequationofstate.h"

HPREquationOfState::HPREquationOfState()
{
}

double HPREquationOfState::computeMicroscaleEnergy(double density, double pressure, HPRMaterialParameters materialParameters)
{
    return HPRMieGruneisen::computeInternalEnergy(density, pressure, materialParameters);
}

double HPREquationOfState::computeMesoscaleDistortionEnergy(vector<vector<double> > distortionTensor, HPRMaterialParameters materialParameters)
{
    vector<vector<double> > strainTensor = TensorAlgebra::computeGramianMatrix(distortionTensor);
    vector<vector<double> > deviatoricStrainTensor = TensorAlgebra::computeDeviator(strainTensor);

    double transverseWaveSpeedSquared = materialParameters.computeTransverseWaveSpeedSquared();

    return (0.25 * transverseWaveSpeedSquared) * MatrixAlgebra::computeComponentSum(MatrixAlgebra::multiplyMatricesElementWise(deviatoricStrainTensor, deviatoricStrainTensor));
}

double HPREquationOfState::computeMesoscaleThermalEnergy(double xThermalImpulse, double yThermalImpulse, double zThermalImpulse, HPRMaterialParameters materialParameters)
{
    double heatWaveSpeedSquared = materialParameters.computeHeatWaveSpeedSquared();
    double thermalImpulseSquared = (xThermalImpulse * xThermalImpulse) + (yThermalImpulse * yThermalImpulse) + (zThermalImpulse * zThermalImpulse);

    return (0.5 * heatWaveSpeedSquared) * thermalImpulseSquared;
}

double HPREquationOfState::computeMacroscaleEnergy(double xVelocity, double yVelocity, double zVelocity)
{
    double velocitySquared = (xVelocity * xVelocity) + (yVelocity * yVelocity) + (zVelocity * zVelocity);

    return (0.5 * velocitySquared);
}

double HPREquationOfState::computeTotalEnergy(double density, double pressure, double xVelocity, double yVelocity, double zVelocity, vector<vector<double> > distortionTensor,
                                              double xThermalImpulse, double yThermalImpulse, double zThermalImpulse, HPRMaterialParameters materialParameters)
{
    double totalEnergy = computeMicroscaleEnergy(density, pressure, materialParameters) + computeMesoscaleDistortionEnergy(distortionTensor, materialParameters) +
            computeMacroscaleEnergy(xVelocity, yVelocity, zVelocity);

    if (materialParameters.getIsThermal())
    {
        totalEnergy += computeMesoscaleThermalEnergy(xThermalImpulse, yThermalImpulse, zThermalImpulse, materialParameters);
    }

    return totalEnergy;
}

double HPREquationOfState::computePressure(double density, double totalEnergy, double xVelocity, double yVelocity, double zVelocity, vector<vector<double> > distortionTensor,
                                           double xThermalImpulse, double yThermalImpulse, double zThermalImpulse, HPRMaterialParameters materialParameters)
{
    double microscaleEnergy = totalEnergy - (computeMacroscaleEnergy(xVelocity, yVelocity, zVelocity) + computeMesoscaleDistortionEnergy(distortionTensor, materialParameters));

    if (materialParameters.getIsThermal())
    {
        microscaleEnergy -= computeMesoscaleThermalEnergy(xThermalImpulse, yThermalImpulse, zThermalImpulse, materialParameters);
    }

    return HPRMieGruneisen::computePressure(density, microscaleEnergy, materialParameters);
}

double HPREquationOfState::computeTemperature(double density, double pressure, HPRMaterialParameters materialParameters)
{
    return HPRMieGruneisen::computeTemperature(density, pressure, materialParameters);
}

vector<double> HPREquationOfState::computeHeatFluxVector(double temperature, double xThermalImpulse, double yThermalImpulse, double zThermalImpulse, HPRMaterialParameters materialParameters)
{
    vector<double> totalEnergyDerivativeThermalImpulse = HPRDerivatives::computeTotalEnergyDerivativeThermalImpulse(xThermalImpulse, yThermalImpulse, zThermalImpulse, materialParameters);

    return VectorAlgebra::multiplyVector(temperature, totalEnergyDerivativeThermalImpulse);
}

vector<vector<double> > HPREquationOfState::computeShearStressTensor(double density, vector<vector<double> > distortionTensor, HPRMaterialParameters materialParameters)
{
    vector<vector<double> > totalEnergyDerivativeDistortionTensor = HPRDerivatives::computeTotalEnergyDerivativeDistortionTensor(distortionTensor, materialParameters);

    return MatrixAlgebra::multiplyMatrix(-density, MatrixAlgebra::multiplyMatrices(MatrixAlgebra::computeTranspose(distortionTensor), totalEnergyDerivativeDistortionTensor));
}

vector<vector<vector<vector<double> > > > HPREquationOfState::computeShearStressTensorDerivativeDistortionTensor(double density, vector<vector<double> > distortionTensor,
                                                                                                                 HPRMaterialParameters materialParameters)
{
    double transverseWaveSpeedSquared = materialParameters.computeTransverseWaveSpeedSquared();

    vector<vector<double> > strainTensor = TensorAlgebra::computeGramianMatrix(distortionTensor);
    vector<vector<vector<vector<double> > > > strainTensorDistortionTensorProduct = TensorAlgebra::multiplyThreeTensors(TensorAlgebra::reshapeMatrixAtLevelOne(strainTensor),
                                                                                                                        TensorAlgebra::reshapeMatrixAtLevelTwo(distortionTensor));

    vector<vector<vector<vector<double> > > > shearStressTensorDerivativeDistortionTensor = TensorAlgebra::subtractFourTensors(
                TensorAlgebra::addFourTensors(TensorAlgebra::swapFourTensorAxes(strainTensorDistortionTensorProduct, 0, 3),
                                              TensorAlgebra::swapFourTensorAxes(strainTensorDistortionTensorProduct, 1, 3)),
                TensorAlgebra::multiplyFourTensor((2.0 / 3.0), strainTensorDistortionTensorProduct));

    vector<vector<double> > distortionTensorDeviatorStrainTensorTranspose = MatrixAlgebra::computeTranspose(
                MatrixAlgebra::multiplyMatrices(distortionTensor,TensorAlgebra::computeDeviator(strainTensor)));

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                shearStressTensorDerivativeDistortionTensor[i][j][k][i] += distortionTensorDeviatorStrainTensorTranspose[j][k];
                shearStressTensorDerivativeDistortionTensor[j][i][k][i] += distortionTensorDeviatorStrainTensorTranspose[j][k];
            }
        }
    }

    return TensorAlgebra::multiplyFourTensor(-(density * transverseWaveSpeedSquared), shearStressTensorDerivativeDistortionTensor);
}

vector<vector<double> > HPREquationOfState::computeShearStressTensorDerivativeDensity(vector<vector<double> > distortionTensor, HPRMaterialParameters materialParameters)
{
    vector<vector<double> > totalEnergyDerivativeDistortionTensor = HPRDerivatives::computeTotalEnergyDerivativeDistortionTensor(distortionTensor, materialParameters);

    return MatrixAlgebra::multiplyMatrix(-1.0, MatrixAlgebra::multiplyMatrices(MatrixAlgebra::computeTranspose(distortionTensor), totalEnergyDerivativeDistortionTensor));
}

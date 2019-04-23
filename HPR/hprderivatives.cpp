#include "hprderivatives.h"

HPRDerivatives::HPRDerivatives()
{
}

double HPRDerivatives::computeTotalEnergyDerivativeDensity(double density, double pressure, HPRMaterialParameters materialParameters)
{
    return HPRMieGruneisen::computeInternalEnergyDerivativeDensity(density, pressure, materialParameters);
}

double HPRDerivatives::computeTotalEnergyDerivativePressure(double density, HPRMaterialParameters materialParameters)
{
    return HPRMieGruneisen::computeInternalEnergyDerivativePressure(density, materialParameters);
}

vector<vector<double> > HPRDerivatives::computeTotalEnergyDerivativeDistortionTensor(vector<vector<double> > distortionTensor, HPRMaterialParameters materialParameters)
{
    vector<vector<double> > strainTensor = TensorAlgebra::computeGramianMatrix(distortionTensor);
    vector<vector<double> > deviatoricStrainTensor = TensorAlgebra::computeDeviator(strainTensor);

    double transverseWaveSpeedSquared = materialParameters.computeTransverseWaveSpeedSquared();

    return MatrixAlgebra::multiplyMatrix(transverseWaveSpeedSquared, MatrixAlgebra::multiplyMatrices(distortionTensor, deviatoricStrainTensor));
}

vector<double> HPRDerivatives::computeTotalEnergyDerivativeThermalImpulse(double xThermalImpulse, double yThermalImpulse, double zThermalImpulse, HPRMaterialParameters materialParameters)
{
    double heatWaveSpeedSquared = materialParameters.computeHeatWaveSpeedSquared();

    vector<double> totalEnergyDerivativeThermalImpulse(3);
    totalEnergyDerivativeThermalImpulse[0] = heatWaveSpeedSquared * xThermalImpulse;
    totalEnergyDerivativeThermalImpulse[1] = heatWaveSpeedSquared * yThermalImpulse;
    totalEnergyDerivativeThermalImpulse[2] = heatWaveSpeedSquared * zThermalImpulse;

    return totalEnergyDerivativeThermalImpulse;
}

double HPRDerivatives::computeTemperatureDerivativeDensity(double density, double pressure, HPRMaterialParameters materialParameters)
{
    return HPRMieGruneisen::computeTemperatureDerivativeDensity(density, pressure, materialParameters);
}

double HPRDerivatives::computeTemperatureDerivativePressure(double density, HPRMaterialParameters materialParameters)
{
    return HPRMieGruneisen::computeTemperatureDerivativePressure(density, materialParameters);
}

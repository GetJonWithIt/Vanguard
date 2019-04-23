#include "hprsourceterms.h"

HPRSourceTerms::HPRSourceTerms()
{
}

double HPRSourceTerms::computeTheta1Reciprocal(double density, vector<vector<double> > distortionTensor, HPRMaterialParameters materialParameters)
{
    double strainDissipationTime = materialParameters.getStrainDissipationTime();
    double transverseWaveSpeedSquared = materialParameters.computeTransverseWaveSpeedSquared();

    if (materialParameters.getIsPlastic())
    {
        double elasticPlasticTransitionParameter = materialParameters.getElasticPlasticTransitionParameter();
        double powerLawIndex = materialParameters.getPowerLawIndex();

        vector<vector<double> > shearStressTensor = HPREquationOfState::computeShearStressTensor(density, distortionTensor, materialParameters);
        double shearStressTensorNorm = TensorAlgebra::computeSigmaNorm(shearStressTensor);
        shearStressTensorNorm = min(shearStressTensorNorm, pow(10.0, 8.0));

        return 3.0 * (pow(MatrixAlgebra::computeDeterminant(distortionTensor), (5.0 / 3.0)) / (transverseWaveSpeedSquared * strainDissipationTime)) *
                pow(shearStressTensorNorm / elasticPlasticTransitionParameter, powerLawIndex);
    }
    else
    {
        return 3.0 * (pow(MatrixAlgebra::computeDeterminant(distortionTensor), (5.0 / 3.0)) / (transverseWaveSpeedSquared * strainDissipationTime));
    }
}

double HPRSourceTerms::computeTheta2Reciprocal(double density, double temperature, HPRMaterialParameters materialParameters)
{
    double referenceDensity = materialParameters.getReferenceDensity();
    double initialTemperature = materialParameters.getInitialTemperature();

    double heatWaveSpeedSquared = materialParameters.computeHeatWaveSpeedSquared();
    double thermalImpulseRelaxationTime = materialParameters.getThermalImpulseRelaxationTime();

    return 1.0 / (heatWaveSpeedSquared * thermalImpulseRelaxationTime * ((density / referenceDensity) * (initialTemperature / temperature)));
}

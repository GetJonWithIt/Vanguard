#include "elasticequationofstate.h"

ElasticEquationOfState::ElasticEquationOfState()
{
}

double ElasticEquationOfState::computeHydrodynamicInternalEnergyFirstComponent(double fingerTensorThirdInvariant, double bulkSoundSpeedSquared, double alphaParameter)
{
    return (bulkSoundSpeedSquared / (2.0 * (alphaParameter * alphaParameter))) * (pow(fingerTensorThirdInvariant, (0.5 * alphaParameter)) - 1.0) * (pow(fingerTensorThirdInvariant,
                                                                                                                                                        (0.5 * alphaParameter)) - 1.0);
}

double ElasticEquationOfState::computeHydrodynamicInternalEnergySecondComponent(double fingerTensorThirdInvariant, double entropy, double specificHeatCapacity, double initialTemperature,
                                                                                double gammaParameter)
{
    return specificHeatCapacity * initialTemperature * pow(fingerTensorThirdInvariant, (0.5 * gammaParameter)) * (exp(entropy / specificHeatCapacity) - 1.0);
}

double ElasticEquationOfState::computeShearDeformationInternalEnergy(double fingerTensorFirstInvariant, double fingerTensorSecondInvariant, double fingerTensorThirdInvariant,
                                                                     double shearWaveSpeedSquared, double betaParameter)
{
    return (0.5 * shearWaveSpeedSquared) * pow(fingerTensorThirdInvariant, (0.5 * betaParameter)) * (((fingerTensorFirstInvariant * fingerTensorFirstInvariant) / 3.0) -
                                                                                                     fingerTensorSecondInvariant);
}

vector<vector<double> > ElasticEquationOfState::computeHydrodynamicInternalEnergyFirstComponentDerivative(double fingerTensorThirdInvariant, double bulkSoundSpeedSquared, double alphaParameter)
{
    return MatrixAlgebra::multiplyMatrix((bulkSoundSpeedSquared / (2.0 * alphaParameter)) * (pow(fingerTensorThirdInvariant, alphaParameter) -
                                                                                             pow(fingerTensorThirdInvariant, (0.5 * alphaParameter))), MatrixAlgebra::computeIdentityMatrix(3));
}

vector<vector<double> > ElasticEquationOfState::computeHydrodynamicInternalEnergySecondComponentDerivative(double fingerTensorThirdInvariant, double entropy, double specificHeatCapacity,
                                                                                                           double initialTemperature, double gammaParameter)
{
    return MatrixAlgebra::multiplyMatrix(specificHeatCapacity * initialTemperature * (0.5 * gammaParameter) * (exp(entropy / specificHeatCapacity) - 1.0) *
                                         pow(fingerTensorThirdInvariant, (0.5 * gammaParameter)), MatrixAlgebra::computeIdentityMatrix(3));
}

vector<vector<double> > ElasticEquationOfState::computeShearDeformationInternalEnergyDerivative(vector<vector<double> > fingerTensor, double fingerTensorFirstInvariant,
                                                                                                double fingerTensorSecondInvariant, double fingerTensorThirdInvariant,
                                                                                                double shearWaveSpeedSquared, double betaParameter)
{
    return MatrixAlgebra::multiplyMatrix(
                (0.5 * shearWaveSpeedSquared) * pow(fingerTensorThirdInvariant, (0.5 * betaParameter)), MatrixAlgebra::addMatrices(
                    MatrixAlgebra::subtractMatrices(MatrixAlgebra::multiplyMatrix(
                                                        (0.5 * betaParameter) * (((fingerTensorFirstInvariant * fingerTensorFirstInvariant) / 3.0) - fingerTensorSecondInvariant),
                                                        MatrixAlgebra::computeIdentityMatrix(3)), MatrixAlgebra::multiplyMatrix((fingerTensorFirstInvariant / 3.0), fingerTensor)),
                    MatrixAlgebra::multiplyMatrices(fingerTensor, fingerTensor)));
}

double ElasticEquationOfState::computeTotalEnergy(vector<vector<double> > distortionTensor, double entropy, double xVelocity, double yVelocity, double zVelocity,
                                                  HyperelasticMaterialParameters materialParameters)
{
    vector<vector<double> > fingerTensor = TensorAlgebra::computeGramianMatrix(distortionTensor);

    double fingerTensorFirstInvariant = TensorAlgebra::computeFirstInvariant(fingerTensor);
    double fingerTensorSecondInvariant = TensorAlgebra::computeSecondInvariant(fingerTensor);
    double fingerTensorThirdInvariant = TensorAlgebra::computeThirdInvariant(fingerTensor);

    double bulkSoundSpeedSquared = materialParameters.computeBulkSoundSpeedSquared();
    double shearWaveSpeedSquared = materialParameters.computeShearWaveSpeedSquared();
    double specificHeatCapacity = materialParameters.getSpecificHeatCapacity();
    double initialTemperature = materialParameters.getInitialTemperature();

    double alphaParameter = materialParameters.getAlphaParameter();
    double betaParameter = materialParameters.getBetaParameter();
    double gammaParameter = materialParameters.getGammaParameter();

    double internalEnergy = computeHydrodynamicInternalEnergyFirstComponent(fingerTensorThirdInvariant, bulkSoundSpeedSquared, alphaParameter) +
            computeHydrodynamicInternalEnergySecondComponent(fingerTensorThirdInvariant, entropy, specificHeatCapacity, initialTemperature, gammaParameter) +
            computeShearDeformationInternalEnergy(fingerTensorFirstInvariant, fingerTensorSecondInvariant, fingerTensorThirdInvariant, shearWaveSpeedSquared, betaParameter);
    double velocitySquared = (xVelocity * xVelocity) + (yVelocity * yVelocity) + (zVelocity * zVelocity);

    return internalEnergy + (0.5 * velocitySquared);
}

double ElasticEquationOfState::computeEntropy(double totalEnergy, vector<vector<double> > distortionTensor, double xVelocity, double yVelocity, double zVelocity,
                                              HyperelasticMaterialParameters materialParameters)
{
    vector<vector<double> > fingerTensor = TensorAlgebra::computeGramianMatrix(distortionTensor);

    double fingerTensorFirstInvariant = TensorAlgebra::computeFirstInvariant(fingerTensor);
    double fingerTensorSecondInvariant = TensorAlgebra::computeSecondInvariant(fingerTensor);
    double fingerTensorThirdInvariant = TensorAlgebra::computeThirdInvariant(fingerTensor);

    double bulkSoundSpeedSquared = materialParameters.computeBulkSoundSpeedSquared();
    double shearWaveSpeedSquared = materialParameters.computeShearWaveSpeedSquared();
    double specificHeatCapacity = materialParameters.getSpecificHeatCapacity();
    double initialTemperature = materialParameters.getInitialTemperature();

    double alphaParameter = materialParameters.getAlphaParameter();
    double betaParameter = materialParameters.getBetaParameter();
    double gammaParameter = materialParameters.getGammaParameter();

    double velocitySquared = (xVelocity * xVelocity) + (yVelocity * yVelocity) + (zVelocity * zVelocity);
    double hydrodynamicInternalEnergySecondComponent = totalEnergy -
            (computeHydrodynamicInternalEnergyFirstComponent(fingerTensorThirdInvariant, bulkSoundSpeedSquared, alphaParameter) +
             computeShearDeformationInternalEnergy(fingerTensorFirstInvariant, fingerTensorSecondInvariant, fingerTensorThirdInvariant, shearWaveSpeedSquared, betaParameter) +
             (0.5 * velocitySquared));

    if (hydrodynamicInternalEnergySecondComponent > 0.0)
    {
        return specificHeatCapacity * log(1.0 + (hydrodynamicInternalEnergySecondComponent / (specificHeatCapacity * initialTemperature * pow(fingerTensorThirdInvariant,
                                                                                                                                              (0.5 * gammaParameter)))));
    }
    else
    {
        return 0.0;
    }
}

double ElasticEquationOfState::computePressure(double density, vector<vector<double> > distortionTensor, double entropy, HyperelasticMaterialParameters materialParameters)
{
    vector<vector<double> > fingerTensor = TensorAlgebra::computeGramianMatrix(distortionTensor);

    double fingerTensorFirstInvariant = TensorAlgebra::computeFirstInvariant(fingerTensor);
    double fingerTensorSecondInvariant = TensorAlgebra::computeSecondInvariant(fingerTensor);
    double fingerTensorThirdInvariant = TensorAlgebra::computeThirdInvariant(fingerTensor);

    double bulkSoundSpeedSquared = materialParameters.computeBulkSoundSpeedSquared();
    double shearWaveSpeedSquared = materialParameters.computeShearWaveSpeedSquared();
    double specificHeatCapacity = materialParameters.getSpecificHeatCapacity();
    double initialTemperature = materialParameters.getInitialTemperature();

    double alphaParameter = materialParameters.getAlphaParameter();
    double betaParameter = materialParameters.getBetaParameter();
    double gammaParameter = materialParameters.getGammaParameter();

    return 2.0 * density * (((bulkSoundSpeedSquared / (2.0 * alphaParameter)) * (pow(fingerTensorThirdInvariant, alphaParameter) - pow(fingerTensorThirdInvariant, (0.5 * alphaParameter)))) +
                            (specificHeatCapacity * initialTemperature * (0.5 * gammaParameter) * (exp(entropy / specificHeatCapacity) - 1.0) *
                             pow(fingerTensorThirdInvariant, (0.5 * gammaParameter))) +
                            (((0.5 * shearWaveSpeedSquared) * pow(fingerTensorThirdInvariant, (0.5 * betaParameter))) *
                             ((0.5 * betaParameter) * (((fingerTensorFirstInvariant * fingerTensorFirstInvariant) / 3.0) - fingerTensorSecondInvariant))));
}

vector<vector<double> > ElasticEquationOfState::computeTotalStressTensor(double density, vector<vector<double> > distortionTensor, double entropy, HyperelasticMaterialParameters materialParameters)
{
    vector<vector<double> > fingerTensor = TensorAlgebra::computeGramianMatrix(distortionTensor);

    double fingerTensorFirstInvariant = TensorAlgebra::computeFirstInvariant(fingerTensor);
    double fingerTensorSecondInvariant = TensorAlgebra::computeSecondInvariant(fingerTensor);
    double fingerTensorThirdInvariant = TensorAlgebra::computeThirdInvariant(fingerTensor);

    double bulkSoundSpeedSquared = materialParameters.computeBulkSoundSpeedSquared();
    double shearWaveSpeedSquared = materialParameters.computeShearWaveSpeedSquared();
    double specificHeatCapacity = materialParameters.getSpecificHeatCapacity();
    double initialTemperature = materialParameters.getInitialTemperature();

    double alphaParameter = materialParameters.getAlphaParameter();
    double betaParameter = materialParameters.getBetaParameter();
    double gammaParameter = materialParameters.getGammaParameter();

    return MatrixAlgebra::multiplyMatrix(
                -2.0 * density, MatrixAlgebra::addMatrices(
                    MatrixAlgebra::addMatrices(computeHydrodynamicInternalEnergyFirstComponentDerivative(fingerTensorThirdInvariant, bulkSoundSpeedSquared, alphaParameter),
                                               computeHydrodynamicInternalEnergySecondComponentDerivative(fingerTensorThirdInvariant, entropy, specificHeatCapacity, initialTemperature,
                                                                                                          gammaParameter)),
                    computeShearDeformationInternalEnergyDerivative(fingerTensor, fingerTensorFirstInvariant, fingerTensorSecondInvariant, fingerTensorThirdInvariant, shearWaveSpeedSquared,
                                                                    betaParameter)));
}

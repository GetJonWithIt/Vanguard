#include "elasticacoustictensor.h"

ElasticAcousticTensor::ElasticAcousticTensor()
{
}

vector<vector<double> > ElasticAcousticTensor::computeStressTensorCentredDifference(vector<vector<double> > distortionTensor, double entropy, HyperelasticMaterialParameters materialParameters,
                                                                                    double epsilon, int direction)
{
    vector<vector<double> > stressTensorCentredDifference(3, vector<double>(3));

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            vector<vector<double> > distortionTensorForwardStep = distortionTensor;
            vector<vector<double> > distortionTensorBackwardStep = distortionTensor;

            distortionTensorForwardStep[direction][j] += epsilon;
            distortionTensorBackwardStep[direction][j] -= epsilon;

            double referenceMassDensity = materialParameters.getReferenceMassDensity();
            double densityForwardStep = referenceMassDensity * MatrixAlgebra::computeDeterminant(distortionTensorForwardStep);
            double densityBackwardStep = referenceMassDensity * MatrixAlgebra::computeDeterminant(distortionTensorBackwardStep);

            vector<vector<double> > stressTensorForwardStep = ElasticEquationOfState::computeTotalStressTensor(densityForwardStep, distortionTensorForwardStep, entropy, materialParameters);
            vector<vector<double> > stressTensorBackwardStep = ElasticEquationOfState::computeTotalStressTensor(densityBackwardStep, distortionTensorBackwardStep, entropy, materialParameters);

            stressTensorCentredDifference[i][j] = stressTensorForwardStep[direction][i] - stressTensorBackwardStep[direction][i];
        }
    }

    return stressTensorCentredDifference;
}

vector<vector<double> > ElasticAcousticTensor::computeAcousticTensor(vector<vector<double> > distortionTensor, double entropy, HyperelasticMaterialParameters materialParameters, double epsilon,
                                                                     int direction)
{
    vector<vector<double> > stressTensorDerivative = MatrixAlgebra::multiplyMatrix((1.0 / (2.0 * epsilon)), computeStressTensorCentredDifference(distortionTensor, entropy, materialParameters,
                                                                                                                                                 epsilon, direction));

    return MatrixAlgebra::multiplyMatrix(-1.0, MatrixAlgebra::multiplyMatrices(stressTensorDerivative, distortionTensor));
}

double ElasticAcousticTensor::computeMaximumWaveSpeed(vector<vector<double> > distortionTensor, double entropy, HyperelasticMaterialParameters materialParameters, int direction)
{
    vector<vector<double> > acousticTensor = computeAcousticTensor(distortionTensor, entropy, materialParameters, pow(10.0, -8.0), direction);

    double referenceMassDensity = materialParameters.getReferenceMassDensity();
    double density = referenceMassDensity * MatrixAlgebra::computeDeterminant(distortionTensor);

    return sqrt(MatrixAlgebra::computeLargestEigenvalue(acousticTensor) / density);
}

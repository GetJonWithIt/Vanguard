#include "slopelimiters.h"

SlopeLimiters::SlopeLimiters()
{
}

double SlopeLimiters::computeRCoefficient(double steepness, double bias)
{
    return 2.0 / ((1.0 - bias) + ((1.0 + bias) * steepness));
}

double SlopeLimiters::computeSuperBeeSlopeLimiter(double steepness, double bias)
{
    if (steepness < 0.0)
    {
        return 0.0;
    }
    else if (steepness < 0.5)
    {
        return 2.0 * steepness;
    }
    else if (steepness < 1.0)
    {
        return 1.0;
    }
    else
    {
        return min(min(steepness, computeRCoefficient(steepness, bias)), 2.0);
    }
}

double SlopeLimiters::computeVanLeerSlopeLimiter(double steepness, double bias)
{
    if (steepness < 0.0)
    {
        return 0.0;
    }
    else
    {
        return min((2.0 * steepness) / (1.0 + steepness), computeRCoefficient(steepness, bias));
    }
}

double SlopeLimiters::computeMinBeeSlopeLimiter(double steepness, double bias)
{
    if (steepness < 0.0)
    {
        return 0.0;
    }
    else if (steepness < 1.0)
    {
        return steepness;
    }
    else
    {
        return min(1.0, computeRCoefficient(steepness, bias));
    }
}

double SlopeLimiters::computeSlopeLimiter(double steepness, double bias, int slopeLimiter)
{
    if (slopeLimiter == 0)
    {
        return computeSuperBeeSlopeLimiter(steepness, bias);
    }
    else if (slopeLimiter == 1)
    {
        return computeVanLeerSlopeLimiter(steepness, bias);
    }
    else
    {
        return computeMinBeeSlopeLimiter(steepness, bias);
    }
}

vector<double> SlopeLimiters::computeSlopeVector(EulerStateVector leftStateVector, EulerStateVector middleStateVector, EulerStateVector rightStateVector, double bias, int slopeLimiter,
                                                 EulerMaterialParameters materialParameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(materialParameters);
    vector<double> middleConservedVariableVector= middleStateVector.computeConservedVariableVector(materialParameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(materialParameters);

    return computeSlopeVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector, bias, slopeLimiter);
}

vector<double> SlopeLimiters::computeSlopeVector(EulerMultiphysicsStateVector leftStateVector, EulerMultiphysicsStateVector middleStateVector, EulerMultiphysicsStateVector rightStateVector,
                                                 double bias, int slopeLimiter, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    return computeSlopeVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector, bias, slopeLimiter);
}

vector<double> SlopeLimiters::computeSlopeVector(EulerReducedStateVector leftStateVector, EulerReducedStateVector middleStateVector, EulerReducedStateVector rightStateVector, double bias,
                                                 int slopeLimiter, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    return computeSlopeVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector, bias, slopeLimiter);
}

vector<double> SlopeLimiters::computeSlopeVector(ElasticStateVector leftStateVector, ElasticStateVector middleStateVector, ElasticStateVector rightStateVector, double bias, int slopeLimiter,
                                                 HyperelasticMaterialParameters materialParameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(materialParameters);
    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(materialParameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(materialParameters);

    return computeSlopeVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector, bias, slopeLimiter);
}

vector<double> SlopeLimiters::computeSlopeVector(ElasticMultiphysicsStateVector leftStateVector, ElasticMultiphysicsStateVector middleStateVector,
                                                 ElasticMultiphysicsStateVector rightStateVector, double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters,
                                                 HyperelasticMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    return computeSlopeVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector, bias, slopeLimiter);
}

vector<double> SlopeLimiters::computeSlopeVector(ElasticReducedStateVector leftStateVector, ElasticReducedStateVector middleStateVector, ElasticReducedStateVector rightStateVector,
                                                 double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    return computeSlopeVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector, bias, slopeLimiter);
}

vector<double> SlopeLimiters::computeSlopeVector(MHDStateVector leftStateVector, MHDStateVector middleStateVector, MHDStateVector rightStateVector, double bias, int slopeLimiter,
                                                 MHDMaterialParameters materialParameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(materialParameters);
    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(materialParameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(materialParameters);

    return computeSlopeVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector, bias, slopeLimiter);
}

vector<double> SlopeLimiters::computeSlopeVector(MHDMultiphysicsStateVector leftStateVector, MHDMultiphysicsStateVector middleStateVector, MHDMultiphysicsStateVector rightStateVector,
                                                 double bias, int slopeLimiter, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    return computeSlopeVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector, bias, slopeLimiter);
}

vector<double> SlopeLimiters::computeSlopeVector(MHDIntermediateStateVector leftStateVector, MHDIntermediateStateVector middleStateVector, MHDIntermediateStateVector rightStateVector,
                                                 double bias, int slopeLimiter, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    return computeSlopeVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector, bias, slopeLimiter);
}

vector<double> SlopeLimiters::computeSlopeVector(MHDReducedStateVector leftStateVector, MHDReducedStateVector middleStateVector, MHDReducedStateVector rightStateVector, double bias,
                                                 int slopeLimiter, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    return computeSlopeVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector, bias, slopeLimiter);
}

vector<double> SlopeLimiters::computeSlopeVector(HPRStateVector leftStateVector, HPRStateVector middleStateVector, HPRStateVector rightStateVector, double bias, int slopeLimiter,
                                                 HPRMaterialParameters materialParameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(materialParameters);
    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(materialParameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(materialParameters);

    return computeSlopeVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector, bias, slopeLimiter);
}

vector<double> SlopeLimiters::computeSlopeVector(HPRMultiphysicsStateVector leftStateVector, HPRMultiphysicsStateVector middleStateVector, HPRMultiphysicsStateVector rightStateVector,
                                                 double bias, int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    return computeSlopeVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector, bias, slopeLimiter);
}

vector<double> SlopeLimiters::computeSlopeVector(HPRIntermediateStateVector leftStateVector, HPRIntermediateStateVector middleStateVector, HPRIntermediateStateVector rightStateVector,
                                                 double bias, int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    return computeSlopeVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector, bias, slopeLimiter);
}

vector<double> SlopeLimiters::computeSlopeVector(HPRReducedStateVector leftStateVector, HPRReducedStateVector middleStateVector, HPRReducedStateVector rightStateVector, double bias,
                                                 int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> leftConservedVariableVector = leftStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> rightConservedVariableVector = rightStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);

    return computeSlopeVector(leftConservedVariableVector, middleConservedVariableVector, rightConservedVariableVector, bias, slopeLimiter);
}

vector<double> SlopeLimiters::computeSlopeVector(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector, vector<double> rightConservedVariableVector,
                                                 double bias, int slopeLimiter)
{
    vector<double> leftConservedVariableVectorDifference = VectorAlgebra::subtractVectors(middleConservedVariableVector, leftConservedVariableVector);
    vector<double> rightConservedVariableVectorDifference = VectorAlgebra::subtractVectors(rightConservedVariableVector, middleConservedVariableVector);

    vector<double> slopeVector = VectorAlgebra::addVectors(VectorAlgebra::multiplyVector(0.5 * (1.0 + bias), leftConservedVariableVectorDifference),
                                                           VectorAlgebra::multiplyVector(0.5 * (1.0 - bias), rightConservedVariableVectorDifference));
    int componentCount = slopeVector.size();

    for (int i = 0; i < componentCount; i++)
    {
        double numerator = leftConservedVariableVectorDifference[i];
        double denominator = rightConservedVariableVectorDifference[i];

        if (abs(numerator) < pow(10.0, -5.0))
        {
            numerator = pow(10.0, -5.0);
        }
        if (abs(denominator) < pow(10.0, -5.0))
        {
            denominator = pow(10.0, -5.0);
        }

        double steepness = numerator / denominator;
        slopeVector[i] *= computeSlopeLimiter(steepness, bias, slopeLimiter);
    }

    return slopeVector;
}

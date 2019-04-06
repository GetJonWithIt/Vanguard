#ifndef SLOPELIMITERS_H
#define SLOPELIMITERS_H

#include "Euler/eulerstatevector.h"
#include "Euler/Multiphysics/eulermultiphysicsstatevector.h"
#include "Elasticity/elasticstatevector.h"
#include "Elasticity/Multiphysics/elasticmultiphysicsstatevector.h"
#include "Elasticity/Multiphysics/elasticreducedstatevector.h"
#include "Mathematics/vectoralgebra.h"
using namespace std;

class SlopeLimiters
{
public:
    SlopeLimiters();

    static double computeRCoefficient(double steepness, double bias);

    static double computeSuperBeeSlopeLimiter(double steepness, double bias);
    static double computeVanLeerSlopeLimiter(double steepness, double bias);
    static double computeMinBeeSlopeLimiter(double steepness, double bias);

    static double computeSlopeLimiter(double steepness, double bias, int slopeLimiter);

    static vector<double> computeSlopeVector(EulerStateVector leftStateVector, EulerStateVector middleStateVector, EulerStateVector rightStateVector, double bias, int slopeLimiter,
                                             EulerMaterialParameters materialParameters);
    static vector<double> computeSlopeVector(EulerMultiphysicsStateVector leftStateVector, EulerMultiphysicsStateVector middleStateVector, EulerMultiphysicsStateVector rightStateVector,
                                             double bias, int slopeLimiter, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> computeSlopeVector(ElasticStateVector leftStateVector, ElasticStateVector middleStateVector, ElasticStateVector rightStateVector, double bias, int slopeLimiter,
                                             HyperelasticMaterialParameters materialParameters);
    static vector<double> computeSlopeVector(ElasticMultiphysicsStateVector leftStateVector, ElasticMultiphysicsStateVector middleStateVector, ElasticMultiphysicsStateVector rightStateVector,
                                             double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static vector<double> computeSlopeVector(ElasticReducedStateVector leftStateVector, ElasticReducedStateVector middleStateVector, ElasticReducedStateVector rightStateVector, double bias,
                                             int slopeLimiter, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeSlopeVector(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector, vector<double> rightConservedVariableVector,
                                             double bias, int slopeLimiter);
};

#endif // SLOPELIMITERS_H

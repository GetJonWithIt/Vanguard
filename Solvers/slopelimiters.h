#ifndef SLOPELIMITERS_H
#define SLOPELIMITERS_H

#include "Euler/eulerstatevector.h"
#include "Euler/Multiphysics/eulermultiphysicsstatevector.h"
#include "Euler/Multiphysics/eulerreducedstatevector.h"
#include "Elasticity/elasticstatevector.h"
#include "Elasticity/Multiphysics/elasticmultiphysicsstatevector.h"
#include "Elasticity/Multiphysics/elasticreducedstatevector.h"
#include "MHD/mhdstatevector.h"
#include "MHD/Multiphysics/mhdmultiphysicsstatevector.h"
#include "MHD/Multiphysics/mhdintermediatestatevector.h"
#include "MHD/Multiphysics/mhdreducedstatevector.h"
#include "HPR/hprstatevector.h"
#include "HPR/Multiphysics/hprmultiphysicsstatevector.h"
#include "HPR/Multiphysics/hprintermediatestatevector.h"
#include "HPR/Multiphysics/hprreducedstatevector.h"
#include "Mathematics/vectoralgebra.h"

#include <omp.h>
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
    static vector<double> computeSlopeVector(EulerReducedStateVector leftStateVector, EulerReducedStateVector middleStateVector, EulerReducedStateVector rightStateVector, double bias,
                                             int slopeLimiter, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> computeSlopeVector(ElasticStateVector leftStateVector, ElasticStateVector middleStateVector, ElasticStateVector rightStateVector, double bias, int slopeLimiter,
                                             HyperelasticMaterialParameters materialParameters);
    static vector<double> computeSlopeVector(ElasticMultiphysicsStateVector leftStateVector, ElasticMultiphysicsStateVector middleStateVector, ElasticMultiphysicsStateVector rightStateVector,
                                             double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static vector<double> computeSlopeVector(ElasticReducedStateVector leftStateVector, ElasticReducedStateVector middleStateVector, ElasticReducedStateVector rightStateVector, double bias,
                                             int slopeLimiter, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeSlopeVector(MHDStateVector leftStateVector, MHDStateVector middleStateVector, MHDStateVector rightStateVector, double bias, int slopeLimiter,
                                             MHDMaterialParameters materialParameters);
    static vector<double> computeSlopeVector(MHDMultiphysicsStateVector leftStateVector, MHDMultiphysicsStateVector middleStateVector, MHDMultiphysicsStateVector rightStateVector,
                                             double bias, int slopeLimiter, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static vector<double> computeSlopeVector(MHDIntermediateStateVector leftStateVector, MHDIntermediateStateVector middleStateVector, MHDIntermediateStateVector rightStateVector,
                                             double bias, int slopeLimiter, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static vector<double> computeSlopeVector(MHDReducedStateVector leftStateVector, MHDReducedStateVector middleStateVector, MHDReducedStateVector rightStateVector, double bias,
                                             int slopeLimiter, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static vector<double> computeSlopeVector(HPRStateVector leftStateVector, HPRStateVector middleStateVector, HPRStateVector rightStateVector, double bias, int slopeLimiter,
                                             HPRMaterialParameters materialParameters);
    static vector<double> computeSlopeVector(HPRMultiphysicsStateVector leftStateVector, HPRMultiphysicsStateVector middleStateVector, HPRMultiphysicsStateVector rightStateVector, double bias,
                                             int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static vector<double> computeSlopeVector(HPRIntermediateStateVector leftStateVector, HPRIntermediateStateVector middleStateVector, HPRIntermediateStateVector rightStateVector, double bias,
                                             int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static vector<double> computeSlopeVector(HPRReducedStateVector leftStateVector, HPRReducedStateVector middleStateVector, HPRReducedStateVector rightStateVector, double bias,
                                             int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static vector<double> computeSlopeVector(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector, vector<double> rightConservedVariableVector,
                                             double bias, int slopeLimiter);
};

#endif // SLOPELIMITERS_H

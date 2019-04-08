#ifndef ELASTICSOLVERS_H
#define ELASTICSOLVERS_H

#include "secondordersolver.h"
using namespace std;

class ElasticSolvers
{
public:
    ElasticSolvers();

    static vector<ElasticStateVector> insertBoundaryCells(vector<ElasticStateVector> & currentCells, int boundarySize);
    static vector<ElasticMultiphysicsStateVector> insertBoundaryCells(vector<ElasticMultiphysicsStateVector> & currentCells, int boundarySize);
    static vector<ElasticReducedStateVector> insertBoundaryCells(vector<ElasticReducedStateVector> & currentCells, int boundarySize);

    static vector<vector<ElasticStateVector> > insertBoundaryCells2D(vector<vector<ElasticStateVector> > & currentCells, int boundarySize);
    static vector<vector<ElasticReducedStateVector> > insertBoundaryCells2D(vector<vector<ElasticReducedStateVector> > & currentCells, int boundarySize);

    static double computeMaximumWaveSpeed(vector<ElasticStateVector> & currentCells, HyperelasticMaterialParameters materialParameters);
    static double computeMaximumWaveSpeed(vector<ElasticMultiphysicsStateVector> & currentCells, HyperelasticMaterialParameters material1Parameters,
                                          HyperelasticMaterialParameters material2Parameters);
    static double computeMaximumWaveSpeed(vector<ElasticReducedStateVector> & currentCells, HyperelasticMaterialParameters material1Parameters,
                                          HyperelasticMaterialParameters material2Parameters);

    static double computeMaximumWaveSpeed2D(vector<vector<ElasticStateVector> > & currentCells, HyperelasticMaterialParameters materialParameters);
    static double computeMaximumWaveSpeed2D(vector<vector<ElasticReducedStateVector> > & currentCells, HyperelasticMaterialParameters material1Parameters,
                                            HyperelasticMaterialParameters material2Parameters);

    static double computeStableTimeStep(vector<ElasticStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                        HyperelasticMaterialParameters materialParameters);
    static double computeStableTimeStep(vector<ElasticMultiphysicsStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                        int currentIteration, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static double computeStableTimeStep(vector<ElasticReducedStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                        HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static double computeStableTimeStep2D(vector<vector<ElasticStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                          int currentIteration, HyperelasticMaterialParameters materialParameters);
    static double computeStableTimeStep2D(vector<vector<ElasticReducedStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                          int currentIteration, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static ElasticStateVector evolveStateByHalfXTimeStep(ElasticStateVector leftStateVector, ElasticStateVector middleStateVector, ElasticStateVector rightStateVector, double cellSpacing,
                                                         double timeStep, double bias, int slopeLimiter, int side, HyperelasticMaterialParameters materialParameters);
    static ElasticMultiphysicsStateVector evolveStateByHalfXTimeStep(ElasticMultiphysicsStateVector leftStateVector, ElasticMultiphysicsStateVector middleStateVector,
                                                                     ElasticMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                     int side, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static ElasticReducedStateVector evolveStateByHalfXTimeStep(ElasticReducedStateVector leftStateVector, ElasticReducedStateVector middleStateVector, ElasticReducedStateVector rightStateVector,
                                                                double cellSpacing, double timeStep, double bias, int slopeLimiter, int side, HyperelasticMaterialParameters material1Parameters,
                                                                HyperelasticMaterialParameters material2Parameters);

    static ElasticStateVector evolveStateByHalfYTimeStep(ElasticStateVector topStateVector, ElasticStateVector middleStateVector, ElasticStateVector bottomStateVector, double cellSpacing,
                                                         double timeStep, double bias, int slopeLimiter, int side, HyperelasticMaterialParameters materialParameters);
    static ElasticReducedStateVector evolveStateByHalfYTimeStep(ElasticReducedStateVector topStateVector, ElasticReducedStateVector middleStateVector,
                                                                ElasticReducedStateVector bottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter, int side,
                                                                HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static ElasticStateVector evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                        HyperelasticMaterialParameters materialParameters);
    static ElasticMultiphysicsStateVector evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                                    HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static ElasticReducedStateVector evolveStateByHalfTimeStepReduced(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                                      HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
};

#endif // ELASTICSOLVERS_H

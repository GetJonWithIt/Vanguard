#ifndef ELASTICSECONDORDERSOLVER_H
#define ELASTICSECONDORDERSOLVER_H

#include "elasticfirstordersolver.h"

class ElasticSecondOrderSolver
{
public:
    ElasticSecondOrderSolver();

    static vector<double> computeXSLICFlux(ElasticStateVector leftLeftStateVector, ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                           ElasticStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           HyperelasticMaterialParameters materialParameters);
    static vector<double> computeXSLICFlux(ElasticMultiphysicsStateVector leftLeftStateVector, ElasticMultiphysicsStateVector leftStateVector,
                                           ElasticMultiphysicsStateVector rightStateVector, ElasticMultiphysicsStateVector rightRightStateVector, double cellSpacing, double timeStep,
                                           double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static vector<double> computeXSLICFlux(ElasticReducedStateVector leftLeftStateVector, ElasticReducedStateVector leftStateVector, ElasticReducedStateVector rightStateVector,
                                           ElasticReducedStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeYSLICFlux(ElasticStateVector topTopStateVector, ElasticStateVector topStateVector, ElasticStateVector bottomStateVector,
                                           ElasticStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           HyperelasticMaterialParameters materialParameters);
    static vector<double> computeYSLICFlux(ElasticReducedStateVector topTopStateVector, ElasticReducedStateVector topStateVector, ElasticReducedStateVector bottomStateVector,
                                           ElasticReducedStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static void computeSLICTimeStep(vector<ElasticStateVector> & currentCells, vector<ElasticStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                    double bias, int slopeLimiter, HyperelasticMaterialParameters materialParameters);
    static void computeSLICTimeStep(vector<ElasticMultiphysicsStateVector> & currentCells, vector<ElasticMultiphysicsStateVector> & currentCellsWithBoundary, double cellSpacing,
                                    double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters,
                                    HyperelasticMaterialParameters material2Parameters);
    static void computeSLICTimeStep(vector<ElasticReducedStateVector> & currentCells, vector<ElasticReducedStateVector> & currentCellsWithBoundary, double cellSpacing,
                                    double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters,
                                    HyperelasticMaterialParameters material2Parameters);

    static void computeXSLICTimeStep2D(vector<vector<ElasticStateVector> > & currentCells, vector<vector<ElasticStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                       double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters materialParameters);
    static void computeXSLICTimeStep2D(vector<vector<ElasticReducedStateVector> > & currentCells, vector<vector<ElasticReducedStateVector> > & currentCellsWithBoundary,
                                       double cellSpacing, double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters,
                                       HyperelasticMaterialParameters material2Parameters);

    static void computeYSLICTimeStep2D(vector<vector<ElasticStateVector> > & currentCells, vector<vector<ElasticStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                       double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters materialParameters);
    static void computeYSLICTimeStep2D(vector<vector<ElasticReducedStateVector> > & currentCells, vector<vector<ElasticReducedStateVector> > & currentCellsWithBoundary,
                                       double cellSpacing, double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters,
                                       HyperelasticMaterialParameters material2Parameters);

    static vector<ElasticStateVector> solve(vector<ElasticStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                            int subcyclingIterations, HyperelasticMaterialParameters materialParameters);
    static vector<ElasticMultiphysicsStateVector> solve(vector<ElasticMultiphysicsStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                        double bias, int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency,
                                                        HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static vector<ElasticReducedStateVector> solve(vector<ElasticReducedStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                   int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, HyperelasticMaterialParameters material1Parameters,
                                                   HyperelasticMaterialParameters material2Parameters);

    static vector<vector<ElasticStateVector> > solve2D(vector<vector<ElasticStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                       int slopeLimiter, int subcyclingIterations, HyperelasticMaterialParameters materialParameters);
    static vector<vector<ElasticReducedStateVector> > solve2D(vector<vector<ElasticReducedStateVector> > & initialCells, double cellSpacing, double CFLCoefficient,
                                                              double finalTime, double bias, int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency,
                                                              HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
};

#endif // ELASTICSECONDORDERSOLVER_H

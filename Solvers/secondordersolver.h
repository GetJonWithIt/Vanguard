#ifndef SECONDORDERSOLVER_H
#define SECONDORDERSOLVER_H

#include "firstordersolver.h"
using namespace std;

class SecondOrderSolver
{
public:
    SecondOrderSolver();

    static vector<double> computeXSLICFlux(EulerStateVector leftLeftStateVector, EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerStateVector rightRightStateVector,
                                           double cellSpacing, double timeStep, double bias, int slopeLimiter, EulerMaterialParameters materialParameters);
    static vector<double> computeXSLICFlux(EulerMultiphysicsStateVector leftLeftStateVector, EulerMultiphysicsStateVector leftStateVector, EulerMultiphysicsStateVector rightStateVector,
                                           EulerMultiphysicsStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> computeXSLICFlux(ElasticStateVector leftLeftStateVector, ElasticStateVector leftStateVector, ElasticStateVector rightStateVector,
                                           ElasticStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           HyperelasticMaterialParameters materialParameters);
    static vector<double> computeXSLICFlux(ElasticMultiphysicsStateVector leftLeftStateVector, ElasticMultiphysicsStateVector leftStateVector, ElasticMultiphysicsStateVector rightStateVector,
                                           ElasticMultiphysicsStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static vector<double> computeXSLICFlux(ElasticReducedStateVector leftLeftStateVector, ElasticReducedStateVector leftStateVector, ElasticReducedStateVector rightStateVector,
                                           ElasticReducedStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeYSLICFlux(EulerStateVector topTopStateVector, EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerStateVector bottomBottomStateVector,
                                           double cellSpacing, double timeStep, double bias, int slopeLimiter, EulerMaterialParameters materialParameters);
    static vector<double> computeYSLICFlux(EulerMultiphysicsStateVector topTopStateVector, EulerMultiphysicsStateVector topStateVector, EulerMultiphysicsStateVector bottomStateVector,
                                           EulerMultiphysicsStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> computeYSLICFlux(ElasticStateVector topTopStateVector, ElasticStateVector topStateVector, ElasticStateVector bottomStateVector,
                                           ElasticStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           HyperelasticMaterialParameters materialParameters);
    static vector<double> computeYSLICFlux(ElasticReducedStateVector topTopStateVector, ElasticReducedStateVector topStateVector, ElasticReducedStateVector bottomStateVector,
                                           ElasticReducedStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static void computeSLICTimeStep(vector<EulerStateVector> & currentCells, vector<EulerStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep, double bias,
                                    int slopeLimiter, EulerMaterialParameters materialParameters);
    static void computeSLICTimeStep(vector<EulerMultiphysicsStateVector> & currentCells, vector<EulerMultiphysicsStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                    double bias, int slopeLimiter, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static void computeSLICTimeStep(vector<ElasticStateVector> & currentCells, vector<ElasticStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep, double bias,
                                    int slopeLimiter, HyperelasticMaterialParameters materialParameters);
    static void computeSLICTimeStep(vector<ElasticMultiphysicsStateVector> & currentCells, vector<ElasticMultiphysicsStateVector> & currentCellsWithBoundary, double cellSpacing,
                                    double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static void computeSLICTimeStep(vector<ElasticReducedStateVector> & currentCells, vector<ElasticReducedStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                    double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static void computeXSLICTimeStep2D(vector<vector<EulerStateVector> > & currentCells, vector<vector<EulerStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                       double bias, int slopeLimiter, EulerMaterialParameters materialParameters);
    static void computeXSLICTimeStep2D(vector<vector<EulerMultiphysicsStateVector> > & currentCells, vector<vector<EulerMultiphysicsStateVector> > & currentCellsWithBoundary,
                                       double cellSpacing, double timeStep, double bias, int slopeLimiter, EulerMaterialParameters material1Parameters,
                                       EulerMaterialParameters material2Parameters);

    static void computeXSLICTimeStep2D(vector<vector<ElasticStateVector> > & currentCells, vector<vector<ElasticStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                       double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters materialParameters);
    static void computeXSLICTimeStep2D(vector<vector<ElasticReducedStateVector> > & currentCells, vector<vector<ElasticReducedStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                       double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static void computeYSLICTimeStep2D(vector<vector<EulerStateVector> > & currentCells, vector<vector<EulerStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                       double bias, int slopeLimiter, EulerMaterialParameters materialParameters);
    static void computeYSLICTimeStep2D(vector<vector<EulerMultiphysicsStateVector> > & currentCells, vector<vector<EulerMultiphysicsStateVector> > & currentCellsWithBoundary,
                                       double cellSpacing, double timeStep, double bias, int slopeLimiter, EulerMaterialParameters material1Parameters,
                                       EulerMaterialParameters material2Parameters);

    static void computeYSLICTimeStep2D(vector<vector<ElasticStateVector> > & currentCells, vector<vector<ElasticStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                       double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters materialParameters);
    static void computeYSLICTimeStep2D(vector<vector<ElasticReducedStateVector> > & currentCells, vector<vector<ElasticReducedStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                       double timeStep, double bias, int slopeLimiter, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<EulerStateVector> solve(vector<EulerStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                          int subcyclingIterations, EulerMaterialParameters materialParameters);
    static vector<EulerMultiphysicsStateVector> solve(vector<EulerMultiphysicsStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                      int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, EulerMaterialParameters material1Parameters,
                                                      EulerMaterialParameters material2Parameters);

    static vector<ElasticStateVector> solve(vector<ElasticStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                            int subcyclingIterations, HyperelasticMaterialParameters materialParameters);
    static vector<ElasticMultiphysicsStateVector> solve(vector<ElasticMultiphysicsStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                        int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, HyperelasticMaterialParameters material1Parameters,
                                                        HyperelasticMaterialParameters material2Parameters);
    static vector<ElasticReducedStateVector> solve(vector<ElasticReducedStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                                   int subcyclingIterations, int reinitialisationFrequency, HyperelasticMaterialParameters material1Parameters,
                                                   HyperelasticMaterialParameters material2Parameters);

    static vector<vector<EulerStateVector> > solve2D(vector<vector<EulerStateVector> > & initialCells, double cellSpacing, double CFLCofficient, double finalTime, double bias, int slopeLimiter,
                                                     int subcyclingIterations, EulerMaterialParameters materialParameters);
    static vector<vector<EulerMultiphysicsStateVector> > solve2D(vector<vector<EulerMultiphysicsStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                                 double bias, int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency,
                                                                 EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<vector<ElasticStateVector> > solve2D(vector<vector<ElasticStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                       int slopeLimiter, int subcyclingIterations, HyperelasticMaterialParameters materialParameters);
    static vector<vector<ElasticReducedStateVector> > solve2D(vector<vector<ElasticReducedStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                              double bias, int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency,
                                                              HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
};

#endif // SECONDORDERSOLVER_H

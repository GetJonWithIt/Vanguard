#ifndef MHDSECONDORDERSOLVER_H
#define MHDSECONDORDERSOLVER_H

#include "mhdfirstordersolver.h"
using namespace std;

class MHDSecondOrderSolver
{
public:
    MHDSecondOrderSolver();

    static vector<double> computeXSLICFlux(MHDStateVector leftLeftStateVector, MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDStateVector rightRightStateVector,
                                           double cellSpacing, double timeStep, double bias, int slopeLimiter, MHDMaterialParameters materialParameters);
    static vector<double> computeXSLICFlux(MHDMultiphysicsStateVector leftLeftStateVector, MHDMultiphysicsStateVector leftStateVector, MHDMultiphysicsStateVector rightStateVector,
                                           MHDMultiphysicsStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static vector<double> computeXSLICFlux(MHDIntermediateStateVector leftLeftStateVector, MHDIntermediateStateVector leftStateVector, MHDIntermediateStateVector rightStateVector,
                                           MHDIntermediateStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static vector<double> computeXSLICFlux(MHDReducedStateVector leftLeftStateVector, MHDReducedStateVector leftStateVector, MHDReducedStateVector rightStateVector,
                                           MHDReducedStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static vector<double> computeYSLICFlux(MHDStateVector topTopStateVector, MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDStateVector bottomBottomStateVector,
                                           double cellSpacing, double timeStep, double bias, int slopeLimiter, MHDMaterialParameters materialParameters);
    static vector<double> computeYSLICFlux(MHDIntermediateStateVector topTopStateVector, MHDIntermediateStateVector topStateVector, MHDIntermediateStateVector bottomStateVector,
                                           MHDIntermediateStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static vector<double> computeYSLICFlux(MHDReducedStateVector topTopStateVector, MHDReducedStateVector topStateVector, MHDReducedStateVector bottomStateVector,
                                           MHDReducedStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static void computeSLICTimeStep(vector<MHDStateVector> & currentCells, vector<MHDStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep, double bias,
                                    int slopeLimiter, MHDMaterialParameters materialParameters);
    static void computeSLICTimeStep(vector<MHDMultiphysicsStateVector> & currentCells, vector<MHDMultiphysicsStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                    double bias, int slopeLimiter, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static void computeSLICTimeStep(vector<MHDIntermediateStateVector> & currentCells, vector<MHDIntermediateStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                    double bias, int slopeLimiter, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static void computeSLICTimeStep(vector<MHDReducedStateVector> & currentCells, vector<MHDReducedStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                    double bias, int slopeLimiter, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static void computeXSLICTimeStep2D(vector<vector<MHDStateVector> > & currentCells, vector<vector<MHDStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                       double bias, int slopeLimiter, MHDMaterialParameters materialParameters);
    static void computeXSLICTimeStep2D(vector<vector<MHDIntermediateStateVector> > & currentCells, vector<vector<MHDIntermediateStateVector> > & currentCellsWithBoundary,
                                       double cellSpacing, double timeStep, double bias, int slopeLimiter, MHDMaterialParameters material1Parameters,
                                       MHDMaterialParameters material2Parameters);
    static void computeXSLICTimeStep2D(vector<vector<MHDReducedStateVector> > & currentCells, vector<vector<MHDReducedStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                       double timeStep, double bias, int slopeLimiter, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static void computeYSLICTimeStep2D(vector<vector<MHDStateVector> > & currentCells, vector<vector<MHDStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                       double bias, int slopeLimiter, MHDMaterialParameters materialParameters);
    static void computeYSLICTimeStep2D(vector<vector<MHDIntermediateStateVector> > & currentCells, vector<vector<MHDIntermediateStateVector> > & currentCellsWithBoundary,
                                       double cellSpacing, double timeStep, double bias, int slopeLimiter, MHDMaterialParameters material1Parameters,
                                       MHDMaterialParameters material2Parameters);
    static void computeYSLICTimeStep2D(vector<vector<MHDReducedStateVector> > & currentCells, vector<vector<MHDReducedStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                       double timeStep, double bias, int slopeLimiter, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static vector<MHDStateVector> solve(vector<MHDStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                        int subcyclingIterations, MHDMaterialParameters materialParameters);
    static vector<MHDMultiphysicsStateVector> solve(vector<MHDMultiphysicsStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                    int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, MHDMaterialParameters material1Parameters,
                                                    MHDMaterialParameters material2Parameters);
    static vector<MHDIntermediateStateVector> solve(vector<MHDIntermediateStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                    int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, MHDMaterialParameters material1Parameters,
                                                    MHDMaterialParameters material2Parameters);
    static vector<MHDReducedStateVector> solve(vector<MHDReducedStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                               int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, MHDMaterialParameters material1Parameters,
                                               MHDMaterialParameters material2Parameters);

    static vector<vector<MHDStateVector> > solve2D(vector<vector<MHDStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                   int slopeLimiter, int subcyclingIterations, MHDMaterialParameters materialParameters);
    static vector<vector<MHDIntermediateStateVector> > solve2D(vector<vector<MHDIntermediateStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                               double bias, int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency,
                                                               MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static vector<vector<MHDReducedStateVector> > solve2D(vector<vector<MHDReducedStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                          double bias, int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency,
                                                          MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
};

#endif // MHDSECONDORDERSOLVER_H

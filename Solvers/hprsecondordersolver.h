#ifndef HPRSECONDORDERSOLVER_H
#define HPRSECONDORDERSOLVER_H

#include "hprfirstordersolver.h"
using namespace std;

class HPRSecondOrderSolver
{
public:
    HPRSecondOrderSolver();

    static vector<double> computeXSLICFlux(HPRStateVector leftLeftStateVector, HPRStateVector leftStateVector, HPRStateVector rightStateVector, HPRStateVector rightRightStateVector,
                                           double cellSpacing, double timeStep, double bias, int slopeLimiter, HPRMaterialParameters materialParameters);
    static vector<double> computeXSLICFlux(HPRMultiphysicsStateVector leftLeftStateVector, HPRMultiphysicsStateVector leftStateVector, HPRMultiphysicsStateVector rightStateVector,
                                           HPRMultiphysicsStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static vector<double> computeXSLICFlux(HPRIntermediateStateVector leftLeftStateVector, HPRIntermediateStateVector leftStateVector, HPRIntermediateStateVector rightStateVector,
                                           HPRIntermediateStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static vector<double> computeXSLICFlux(HPRReducedStateVector leftLeftStateVector, HPRReducedStateVector leftStateVector, HPRReducedStateVector rightStateVector,
                                           HPRReducedStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static vector<double> computeYSLICFlux(HPRStateVector topTopStateVector, HPRStateVector topStateVector, HPRStateVector bottomStateVector, HPRStateVector bottomBottomStateVector,
                                           double cellSpacing, double timeStep, double bias, int slopeLimiter, HPRMaterialParameters materialParameters);
    static vector<double> computeYSLICFlux(HPRIntermediateStateVector topTopStateVector, HPRIntermediateStateVector topStateVector, HPRIntermediateStateVector bottomStateVector,
                                           HPRIntermediateStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static vector<double> computeYSLICFlux(HPRReducedStateVector topTopStateVector, HPRReducedStateVector topStateVector, HPRReducedStateVector bottomStateVector,
                                           HPRReducedStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static void computeSLICTimeStep(vector<HPRStateVector> & currentCells, vector<HPRStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                    HPRMaterialParameters materialParameters);
    static void computeSLICTimeStep(vector<HPRMultiphysicsStateVector> & currentCells, vector<HPRMultiphysicsStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                    double bias, int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static void computeSLICTimeStep(vector<HPRIntermediateStateVector> & currentCells, vector<HPRIntermediateStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                    double bias, int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static void computeSLICTimeStep(vector<HPRReducedStateVector> & currentCells, vector<HPRReducedStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep, double bias,
                                    int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static void computeXSLICTimeStep2D(vector<vector<HPRStateVector> > & currentCells, vector<vector<HPRStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep, double bias,
                                       int slopeLimiter, HPRMaterialParameters materialParameters);
    static void computeXSLICTimeStep2D(vector<vector<HPRIntermediateStateVector> > & currentCells, vector<vector<HPRIntermediateStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                       double timeStep, double bias, int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static void computeXSLICTimeStep2D(vector<vector<HPRReducedStateVector> > & currentCells, vector<vector<HPRReducedStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                       double timeStep, double bias, int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static void computeYSLICTimeStep2D(vector<vector<HPRStateVector> > & currentCells, vector<vector<HPRStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep, double bias,
                                       int slopeLimiter, HPRMaterialParameters materialParameters);
    static void computeYSLICTimeStep2D(vector<vector<HPRIntermediateStateVector> > & currentCells, vector<vector<HPRIntermediateStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                       double timeStep, double bias, int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static void computeYSLICTimeStep2D(vector<vector<HPRReducedStateVector> > & currentCells, vector<vector<HPRReducedStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                       double timeStep, double bias, int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static vector<HPRStateVector> solve(vector<HPRStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                        int subcyclingIterations, HPRMaterialParameters materialParameters);
    static vector<HPRMultiphysicsStateVector> solve(vector<HPRMultiphysicsStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                    int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, HPRMaterialParameters material1Parameters,
                                                    HPRMaterialParameters material2Parameters);
    static vector<HPRIntermediateStateVector> solve(vector<HPRIntermediateStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                    int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, HPRMaterialParameters material1Parameters,
                                                    HPRMaterialParameters material2Parameters);
    static vector<HPRReducedStateVector> solve(vector<HPRReducedStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                               int subcyclingIterations, int reinitialisationFrequency, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static vector<vector<HPRStateVector> > solve2D(vector<vector<HPRStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                                   int subcyclingIterations, HPRMaterialParameters materialParameters);
    static vector<vector<HPRIntermediateStateVector> > solve2D(vector<vector<HPRIntermediateStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                               double bias, int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency,
                                                               HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static vector<vector<HPRReducedStateVector> > solve2D(vector<vector<HPRReducedStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                          int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, HPRMaterialParameters material1Parameters,
                                                          HPRMaterialParameters material2Parameters);
};

#endif // HPRSECONDORDERSOLVER_H

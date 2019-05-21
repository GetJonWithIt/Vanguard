#ifndef SECONDORDERSOLVER_H
#define SECONDORDERSOLVER_H

#include "firstordersolver.h"
using namespace std;

class SecondOrderSolver
{
public:
    SecondOrderSolver();

    static vector<double> computeXSLICFlux(EulerStateVector leftLeftStateVector, EulerStateVector leftStateVector, EulerStateVector rightStateVector,
                                           EulerStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           EulerMaterialParameters materialParameters);
    static vector<double> computeXSLICFlux(EulerMultiphysicsStateVector leftLeftStateVector, EulerMultiphysicsStateVector leftStateVector,
                                           EulerMultiphysicsStateVector rightStateVector, EulerMultiphysicsStateVector rightRightStateVector, double cellSpacing, double timeStep,
                                           double bias, int slopeLimiter, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static vector<double> computeXSLICFlux(EulerReducedStateVector leftLeftStateVector, EulerReducedStateVector leftStateVector, EulerReducedStateVector rightStateVector,
                                           EulerReducedStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> computeYSLICFlux(EulerStateVector topTopStateVector, EulerStateVector topStateVector, EulerStateVector bottomStateVector,
                                           EulerStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           EulerMaterialParameters materialParameters);
    static vector<double> computeYSLICFlux(EulerMultiphysicsStateVector topTopStateVector, EulerMultiphysicsStateVector topStateVector,
                                           EulerMultiphysicsStateVector bottomStateVector, EulerMultiphysicsStateVector bottomBottomStateVector, double cellSpacing, double timeStep,
                                           double bias, int slopeLimiter, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static vector<double> computeYSLICFlux(EulerReducedStateVector topTopStateVector, EulerReducedStateVector topStateVector, EulerReducedStateVector bottomStateVector,
                                           EulerReducedStateVector bottomBottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static void computeSLICTimeStep(vector<EulerStateVector> & currentCells, vector<EulerStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep, double bias,
                                    int slopeLimiter, EulerMaterialParameters materialParameters);
    static void computeSLICTimeStep(vector<EulerMultiphysicsStateVector> & currentCells, vector<EulerMultiphysicsStateVector> & currentCellsWithBoundary, double cellSpacing,
                                    double timeStep, double bias, int slopeLimiter, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static void computeSLICTimeStep(vector<EulerReducedStateVector> & currentCells, vector<EulerReducedStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                    double bias, int slopeLimiter, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static void computeXSLICTimeStep2D(vector<vector<EulerStateVector> > & currentCells, vector<vector<EulerStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                       double timeStep, double bias, int slopeLimiter, EulerMaterialParameters materialParameters);
    static void computeXSLICTimeStep2D(vector<vector<EulerMultiphysicsStateVector> > & currentCells, vector<vector<EulerMultiphysicsStateVector> > & currentCellsWithBoundary,
                                       double cellSpacing, double timeStep, double bias, int slopeLimiter, EulerMaterialParameters material1Parameters,
                                       EulerMaterialParameters material2Parameters);
    static void computeXSLICTimeStep2D(vector<vector<EulerReducedStateVector> > & currentCells, vector<vector<EulerReducedStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                       double timeStep, double bias, int slopeLimiter, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static void computeYSLICTimeStep2D(vector<vector<EulerStateVector> > & currentCells, vector<vector<EulerStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                       double timeStep, double bias, int slopeLimiter, EulerMaterialParameters materialParameters);
    static void computeYSLICTimeStep2D(vector<vector<EulerMultiphysicsStateVector> > & currentCells, vector<vector<EulerMultiphysicsStateVector> > & currentCellsWithBoundary,
                                       double cellSpacing, double timeStep, double bias, int slopeLimiter, EulerMaterialParameters material1Parameters,
                                       EulerMaterialParameters material2Parameters);
    static void computeYSLICTimeStep2D(vector<vector<EulerReducedStateVector> > & currentCells, vector<vector<EulerReducedStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                       double timeStep, double bias, int slopeLimiter, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameter);

    static vector<EulerStateVector> solve(vector<EulerStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                          int subcyclingIterations, EulerMaterialParameters materialParameters);
    static vector<EulerMultiphysicsStateVector> solve(vector<EulerMultiphysicsStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                      int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, EulerMaterialParameters material1Parameters,
                                                      EulerMaterialParameters material2Parameters);
    static vector<EulerReducedStateVector> solve(vector<EulerReducedStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                 int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, EulerMaterialParameters material1Parameters,
                                                 EulerMaterialParameters material2Parameters);

    static vector<vector<EulerStateVector> > solve2D(vector<vector<EulerStateVector> > & initialCells, double cellSpacing, double CFLCofficient, double finalTime, double bias,
                                                     int slopeLimiter, int subcyclingIterations, EulerMaterialParameters materialParameters);
    static vector<vector<EulerMultiphysicsStateVector> > solve2D(vector<vector<EulerMultiphysicsStateVector> > & initialCells, double cellSpacing, double CFLCoefficient,
                                                                 double finalTime, double bias, int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency,
                                                                 EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static vector<vector<EulerReducedStateVector> > solve2D(vector<vector<EulerReducedStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                            double bias, int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency,
                                                            EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
};

#endif // SECONDORDERSOLVER_H

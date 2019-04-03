#ifndef SOLVERS_H
#define SOLVERS_H

#include "slopelimiters.h"
using namespace std;

class Solvers
{
public:
    Solvers();

    static vector<EulerStateVector> insertBoundaryCells(vector<EulerStateVector> & currentCells, int boundarySize);
    static vector<EulerMultiphysicsStateVector> insertBoundaryCells(vector<EulerMultiphysicsStateVector> & currentCells, int boundarySize);

    static vector<vector<EulerStateVector> > insertBoundaryCells2D(vector<vector<EulerStateVector> > & currentCells, int boundarySize);
    static vector<vector<EulerMultiphysicsStateVector> > insertBoundaryCells2D(vector<vector<EulerMultiphysicsStateVector> > & currentCells, int boundarySize);

    static double computeMaximumWaveSpeed(vector<EulerStateVector> & currentCells, EulerMaterialParameters materialParameters);
    static double computeMaximumWaveSpeed(vector<EulerMultiphysicsStateVector> & currentCells, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static double computeMaximumWaveSpeed2D(vector<vector<EulerStateVector> > & currentCells, EulerMaterialParameters materialParameters);
    static double computeMaximumWaveSpeed2D(vector<vector<EulerMultiphysicsStateVector> > & currentCells, EulerMaterialParameters material1Parameters,
                                            EulerMaterialParameters material2Parameters);

    static double computeStableTimeStep(vector<EulerStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                        EulerMaterialParameters materialParameters);
    static double computeStableTimeStep(vector<EulerMultiphysicsStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                        int currentIteration, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static double computeStableTimeStep2D(vector<vector<EulerStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                          int currentIteration, EulerMaterialParameters materialParameters);
    static double computeStableTimeStep2D(vector<vector<EulerMultiphysicsStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                          int currentIteration, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static double computeStableTimeStep(double timeStep, double currentTime, double finalTime, int currentIteration);

    static EulerStateVector evolveStateByHalfXTimeStep(EulerStateVector leftStateVector, EulerStateVector middleStateVector, EulerStateVector rightStateVector, double cellSpacing, double timeStep,
                                                       double bias, int slopeLimiter, int side, EulerMaterialParameters materialParameters);
    static EulerMultiphysicsStateVector evolveStateByHalfXTimeStep(EulerMultiphysicsStateVector leftStateVector, EulerMultiphysicsStateVector middleStateVector,
                                                                   EulerMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter, int side,
                                                                   EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static EulerStateVector evolveStateByHalfYTimeStep(EulerStateVector topStateVector, EulerStateVector middleStateVector, EulerStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                       double bias, int slopeLimiter, int side, EulerMaterialParameters materialParameters);
    static EulerMultiphysicsStateVector evolveStateByHalfYTimeStep(EulerMultiphysicsStateVector topStateVector, EulerMultiphysicsStateVector middleStateVector,
                                                                   EulerMultiphysicsStateVector bottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                   int side, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static EulerStateVector evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                      EulerMaterialParameters materialParameters);
    static EulerMultiphysicsStateVector evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                                 EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> computeEvolutionVector(vector<double> leftFluxVector, vector<double> rightFluxVector, double cellSpacing, double timeStep);
};

#endif // SOLVERS_H

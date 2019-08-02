#ifndef MHDSOLVERS_H
#define MHDSOLVERS_H

#include "secondordersolver.h"
using namespace std;

class MHDSolvers
{
public:
    MHDSolvers();

    static vector<MHDStateVector> insertBoundaryCells(vector<MHDStateVector> & currentCells, int boundarySize);
    static vector<MHDMultiphysicsStateVector> insertBoundaryCells(vector<MHDMultiphysicsStateVector> & currentCells, int boundarySize);
    static vector<MHDIntermediateStateVector> insertBoundaryCells(vector<MHDIntermediateStateVector> & currentCells, int boundarySize);
    static vector<MHDReducedStateVector> insertBoundaryCells(vector<MHDReducedStateVector> & currentCells, int boundarySize);

    static vector<vector<MHDStateVector> > insertBoundaryCells2D(vector<vector<MHDStateVector> > & currentCells, int boundarySize);
    static vector<vector<MHDIntermediateStateVector> > insertBoundaryCells2D(vector<vector<MHDIntermediateStateVector> > & currentCells, int boundarySize);
    static vector<vector<MHDReducedStateVector> > insertBoundaryCells2D(vector<vector<MHDReducedStateVector> > & currentCells, int boundarySize);

    static double computeMaximumWaveSpeed(vector<MHDStateVector> & currentCells, MHDMaterialParameters materialParameters);
    static double computeMaximumWaveSpeed(vector<MHDMultiphysicsStateVector> & currentCells, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeMaximumWaveSpeed(vector<MHDIntermediateStateVector> & currentCells, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeMaximumWaveSpeed(vector<MHDReducedStateVector> & currentCells, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeMaximumWaveSpeed2D(vector<vector<MHDStateVector> > & currentCells, MHDMaterialParameters materialParameters);
    static double computeMaximumWaveSpeed2D(vector<vector<MHDIntermediateStateVector> > & currentCells, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeMaximumWaveSpeed2D(vector<vector<MHDReducedStateVector> > & currentCells, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeStableTimeStep(vector<MHDStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                        MHDMaterialParameters materialParameters);
    static double computeStableTimeStep(vector<MHDMultiphysicsStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                        int currentIteration, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeStableTimeStep(vector<MHDIntermediateStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                        int currentIteration, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeStableTimeStep(vector<MHDReducedStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                        MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeStableTimeStep2D(vector<vector<MHDStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                          int currentIteration, MHDMaterialParameters materialParameters);
    static double computeStableTimeStep2D(vector<vector<MHDIntermediateStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                          int currentIteration, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeStableTimeStep2D(vector<vector<MHDReducedStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                          int currentIteration, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static MHDStateVector evolveStateByHalfXTimeStep(MHDStateVector leftStateVector, MHDStateVector middleStateVector, MHDStateVector rightStateVector, double cellSpacing, double timeStep,
                                                     double bias, int slopeLimiter, int side, MHDMaterialParameters materialParameters);
    static MHDMultiphysicsStateVector evolveStateByHalfXTimeStep(MHDMultiphysicsStateVector leftStateVector, MHDMultiphysicsStateVector middleStateVector,
                                                                 MHDMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter, int side,
                                                                 MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static MHDIntermediateStateVector evolveStateByHalfXTimeStep(MHDIntermediateStateVector leftStateVector, MHDIntermediateStateVector middleStateVector,
                                                                 MHDIntermediateStateVector rightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter, int side,
                                                                 MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static MHDReducedStateVector evolveStateByHalfXTimeStep(MHDReducedStateVector leftStateVector, MHDReducedStateVector middleStateVector, MHDReducedStateVector rightStateVector,
                                                            double cellSpacing, double timeStep, double bias, int slopeLimiter, int side, MHDMaterialParameters material1Parameters,
                                                            MHDMaterialParameters material2Parameters);

    static MHDStateVector evolveStateByHalfYTimeStep(MHDStateVector topStateVector, MHDStateVector middleStateVector, MHDStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                     double bias, int slopeLimiter, int side, MHDMaterialParameters materialParameters);
    static MHDIntermediateStateVector evolveStateByHalfYTimeStep(MHDIntermediateStateVector topStateVector, MHDIntermediateStateVector middleStateVector,
                                                                 MHDIntermediateStateVector bottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter, int side,
                                                                 MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static MHDReducedStateVector evolveStateByHalfYTimeStep(MHDReducedStateVector topStateVector, MHDReducedStateVector middleStateVector, MHDReducedStateVector bottomStateVector,
                                                            double cellSpacing, double timeStep, double bias, int slopeLimiter, int side, MHDMaterialParameters material1Parameters,
                                                            MHDMaterialParameters material2Parameters);

    static MHDStateVector evolveStateByFractionalXTimeStep(double stepFraction, MHDStateVector leftStateVector, MHDStateVector middleStateVector, MHDStateVector rightStateVector,
                                                           double cellSpacing, double timeStep, double bias, int slopeLimiter, MHDMaterialParameters materialParameters);
    static MHDIntermediateStateVector evolveStateByFractionalXTimeStep(double stepFraction, MHDIntermediateStateVector leftStateVector, MHDIntermediateStateVector middleStateVector,
                                                                       MHDIntermediateStateVector rightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                       MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static MHDReducedStateVector evolveStateByFractionalXTimeStep(double stepFraction, MHDReducedStateVector leftStateVector, MHDReducedStateVector middleStateVector,
                                                                  MHDReducedStateVector rightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                  MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static MHDStateVector evolveStateByFractionalYTimeStep(double stepFraction, MHDStateVector topStateVector, MHDStateVector middleStateVector, MHDStateVector bottomStateVector,
                                                           double cellSpacing, double timeStep, double bias, int slopeLimiter, MHDMaterialParameters materialParameters);
    static MHDIntermediateStateVector evolveStateByFractionalYTimeStep(double stepFraction, MHDIntermediateStateVector topStateVector, MHDIntermediateStateVector middleStateVector,
                                                                       MHDIntermediateStateVector bottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                       MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static MHDReducedStateVector evolveStateByFractionalYTimeStep(double stepFraction, MHDReducedStateVector topStateVector, MHDReducedStateVector middleStateVector,
                                                                  MHDReducedStateVector bottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                  MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static MHDStateVector evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                    MHDMaterialParameters materialParameters);
    static MHDMultiphysicsStateVector evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                                MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static MHDIntermediateStateVector evolveStateByHalfTimeStepIntermediate(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector,
                                                                            int side, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static MHDReducedStateVector evolveStateByHalfTimeStepReduced(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                                  MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static MHDStateVector evolveStateByFractionalTimeStep(vector<double> middleConservedVariableVector, vector<double> conservedVariableVectorEvolution, MHDMaterialParameters materialParameters);
    static MHDIntermediateStateVector evolveStateByFractionalTimeStepIntermediate(vector<double> middleConservedVariableVector, vector<double> conservedVariableVectorEvolution,
                                                                                  MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static MHDReducedStateVector evolveStateByFractionalTimeStepReduced(vector<double> middleConservedVariableVector, vector<double> conservedVariableVectorEvolution,
                                                                        MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
};

#endif // MHDSOLVERS_H

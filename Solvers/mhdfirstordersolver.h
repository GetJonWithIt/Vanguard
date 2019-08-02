#ifndef MHDFIRSTORDERSOLVER_H
#define MHDFIRSTORDERSOLVER_H

#include "mhdforcingsolver.h"
using namespace std;

class MHDFirstOrderSolver
{
public:
    MHDFirstOrderSolver();

    static vector<double> computeXLaxFriedrichsFlux(MHDStateVector leftStateVector, MHDStateVector rightStateVector, double cellSpacing, double timeStep, MHDMaterialParameters materialParameters);
    static vector<double> computeXLaxFriedrichsFlux(MHDMultiphysicsStateVector leftStateVector, MHDMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                                    MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static vector<double> computeXLaxFriedrichsFlux(MHDIntermediateStateVector leftStateVector, MHDIntermediateStateVector rightStateVector, double cellSpacing, double timeStep,
                                                    MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static vector<double> computeXLaxFriedrichsFlux(MHDReducedStateVector leftStateVector, MHDReducedStateVector rightStateVector, double cellSpacing, double timeStep,
                                                    MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static vector<double> computeYLaxFriedrichsFlux(MHDStateVector topStateVector, MHDStateVector bottomStateVector, double cellSpacing, double timeStep, MHDMaterialParameters materialParameters);
    static vector<double> computeYLaxFriedrichsFlux(MHDIntermediateStateVector topStateVector, MHDIntermediateStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                    MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static vector<double> computeYLaxFriedrichsFlux(MHDReducedStateVector topStateVector, MHDReducedStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                    MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static vector<double> computeXRichtmyerFlux(MHDStateVector leftStateVector, MHDStateVector rightStateVector, double cellSpacing, double timeStep, MHDMaterialParameters materialParameters);
    static vector<double> computeXRichtmyerFlux(MHDMultiphysicsStateVector leftStateVector, MHDMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                                MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static vector<double> computeXRichtmyerFlux(MHDIntermediateStateVector leftStateVector, MHDIntermediateStateVector rightStateVector, double cellSpacing, double timeStep,
                                                MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static vector<double> computeXRichtmyerFlux(MHDReducedStateVector leftStateVector, MHDReducedStateVector rightStateVector, double cellSpacing, double timeStep,
                                                MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static vector<double> computeYRichtmyerFlux(MHDStateVector topStateVector, MHDStateVector bottomStateVector, double cellSpacing, double timeStep, MHDMaterialParameters materialParameters);
    static vector<double> computeYRichtmyerFlux(MHDIntermediateStateVector topStateVector, MHDIntermediateStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static vector<double> computeYRichtmyerFlux(MHDReducedStateVector topStateVector, MHDReducedStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static vector<double> computeXFORCEFlux(MHDStateVector leftStateVector, MHDStateVector rightStateVector, double cellSpacing, double timeStep, MHDMaterialParameters materialParameters);
    static vector<double> computeXFORCEFlux(MHDMultiphysicsStateVector leftStateVector, MHDMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                            MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static vector<double> computeXFORCEFlux(MHDIntermediateStateVector leftStateVector, MHDIntermediateStateVector rightStateVector, double cellSpacing, double timeStep,
                                            MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static vector<double> computeXFORCEFlux(MHDReducedStateVector leftStateVector, MHDReducedStateVector rightStateVector, double cellSpacing, double timeStep,
                                            MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static vector<double> computeYFORCEFlux(MHDStateVector topStateVector, MHDStateVector bottomStateVector, double cellSpacing, double timeStep, MHDMaterialParameters materialParameters);
    static vector<double> computeYFORCEFlux(MHDIntermediateStateVector topStateVector, MHDIntermediateStateVector bottomStateVector, double cellSpacing, double timeStep,
                                            MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static vector<double> computeYFORCEFlux(MHDReducedStateVector topStateVector, MHDReducedStateVector bottomStateVector, double cellSpacing, double timeStep,
                                            MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static void computeFORCETimeStep(vector<MHDStateVector> & currentCells, vector<MHDStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                     MHDMaterialParameters materialParameters);
    static void computeFORCETimeStep(vector<MHDMultiphysicsStateVector> & currentCells, vector<MHDMultiphysicsStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                     MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static void computeFORCETimeStep(vector<MHDIntermediateStateVector> & currentCells, vector<MHDIntermediateStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                     MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static void computeFORCETimeStep(vector<MHDReducedStateVector> & currentCells, vector<MHDReducedStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                     MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static void computeXFORCETimeStep2D(vector<vector<MHDStateVector> > & currentCells, vector<vector<MHDStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                        MHDMaterialParameters materialParameters);
    static void computeXFORCETimeStep2D(vector<vector<MHDIntermediateStateVector> > & currentCells, vector<vector<MHDIntermediateStateVector> > & currentCellsWithBoundary,
                                        double cellSpacing, double timeStep, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static void computeXFORCETimeStep2D(vector<vector<MHDReducedStateVector> > & currentCells, vector<vector<MHDReducedStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                        double timeStep, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static void computeYFORCETimeStep2D(vector<vector<MHDStateVector> > & currentCells, vector<vector<MHDStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                        MHDMaterialParameters materialParameters);
    static void computeYFORCETimeStep2D(vector<vector<MHDIntermediateStateVector> > & currentCells, vector<vector<MHDIntermediateStateVector> > & currentCellsWithBoundary,
                                        double cellSpacing, double timeStep, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static void computeYFORCETimeStep2D(vector<vector<MHDReducedStateVector> > & currentCells, vector<vector<MHDReducedStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                        double timeStep, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static vector<MHDStateVector> solve(vector<MHDStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                        MHDMaterialParameters materialParameters);
    static vector<MHDMultiphysicsStateVector> solve(vector<MHDMultiphysicsStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                    int subcyclingIterations, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static vector<MHDIntermediateStateVector> solve(vector<MHDIntermediateStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                    int subcyclingIterations, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static vector<MHDReducedStateVector> solve(vector<MHDReducedStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                               MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static vector<vector<MHDStateVector> > solve2D(vector<vector<MHDStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                                   MHDMaterialParameters materialParameters);
    static vector<vector<MHDIntermediateStateVector> > solve2D(vector<vector<MHDIntermediateStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                               int subcyclingIterations, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static vector<vector<MHDReducedStateVector> > solve2D(vector<vector<MHDReducedStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                          int subcyclingIterations, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
};

#endif // MHDFIRSTORDERSOLVER_H

#ifndef HPRFIRSTORDERSOLVER_H
#define HPRFIRSTORDERSOLVER_H

#include "hprforcingsolver.h"
using namespace std;

class HPRFirstOrderSolver
{
public:
    HPRFirstOrderSolver();

    static vector<double> computeXLaxFriedrichsFlux(HPRStateVector leftStateVector, HPRStateVector rightStateVector, double cellSpacing, double timeStep, HPRMaterialParameters materialParameters);
    static vector<double> computeXLaxFriedrichsFlux(HPRMultiphysicsStateVector leftStateVector, HPRMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                                    HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static vector<double> computeXLaxFriedrichsFlux(HPRIntermediateStateVector leftStateVector, HPRIntermediateStateVector rightStateVector, double cellSpacing, double timeStep,
                                                    HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static vector<double> computeXLaxFriedrichsFlux(HPRReducedStateVector leftStateVector, HPRReducedStateVector rightStateVector, double cellSpacing, double timeStep,
                                                    HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static vector<double> computeYLaxFriedrichsFlux(HPRStateVector topStateVector, HPRStateVector bottomStateVector, double cellSpacing, double timeStep, HPRMaterialParameters materialParameters);
    static vector<double> computeYLaxFriedrichsFlux(HPRIntermediateStateVector topStateVector, HPRIntermediateStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                    HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static vector<double> computeYLaxFriedrichsFlux(HPRReducedStateVector topStateVector, HPRReducedStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                    HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static vector<double> computeXRichtmyerFlux(HPRStateVector leftStateVector, HPRStateVector rightStateVector, double cellSpacing, double timeStep, HPRMaterialParameters materialParameters);
    static vector<double> computeXRichtmyerFlux(HPRMultiphysicsStateVector leftStateVector, HPRMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                                HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static vector<double> computeXRichtmyerFlux(HPRIntermediateStateVector leftStateVector, HPRIntermediateStateVector rightStateVector, double cellSpacing, double timeStep,
                                                HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static vector<double> computeXRichtmyerFlux(HPRReducedStateVector leftStateVector, HPRReducedStateVector rightStateVector, double cellSpacing, double timeStep,
                                                HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static vector<double> computeYRichtmyerFlux(HPRStateVector topStateVector, HPRStateVector bottomStateVector, double cellSpacing, double timeStep, HPRMaterialParameters materialParameters);
    static vector<double> computeYRichtmyerFlux(HPRIntermediateStateVector topStateVector, HPRIntermediateStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static vector<double> computeYRichtmyerFlux(HPRReducedStateVector topStateVector, HPRReducedStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static vector<double> computeXFORCEFlux(HPRStateVector leftStateVector, HPRStateVector rightStateVector, double cellSpacing, double timeStep, HPRMaterialParameters materialParameters);
    static vector<double> computeXFORCEFlux(HPRMultiphysicsStateVector leftStateVector, HPRMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                            HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static vector<double> computeXFORCEFlux(HPRIntermediateStateVector leftStateVector, HPRIntermediateStateVector rightStateVector, double cellSpacing, double timeStep,
                                            HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static vector<double> computeXFORCEFlux(HPRReducedStateVector leftStateVector, HPRReducedStateVector rightStateVector, double cellSpacing, double timeStep,
                                            HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static vector<double> computeYFORCEFlux(HPRStateVector topStateVector, HPRStateVector bottomStateVector, double cellSpacing, double timeStep, HPRMaterialParameters materialParameters);
    static vector<double> computeYFORCEFlux(HPRIntermediateStateVector topStateVector, HPRIntermediateStateVector bottomStateVector, double cellSpacing, double timeStep,
                                            HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static vector<double> computeYFORCEFlux(HPRReducedStateVector topStateVector, HPRReducedStateVector bottomStateVector, double cellSpacing, double timeStep,
                                            HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static void computeFORCETimeStep(vector<HPRStateVector> & currentCells, vector<HPRStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                     HPRMaterialParameters materialParameters);
    static void computeFORCETimeStep(vector<HPRMultiphysicsStateVector> & currentCells, vector<HPRMultiphysicsStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                     HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static void computeFORCETimeStep(vector<HPRIntermediateStateVector> & currentCells, vector<HPRIntermediateStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                     HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static void computeFORCETimeStep(vector<HPRReducedStateVector> & currentCells, vector<HPRReducedStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                     HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static void computeXFORCETimeStep2D(vector<vector<HPRStateVector> > & currentCells, vector<vector<HPRStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                        HPRMaterialParameters materialParameters);
    static void computeXFORCETimeStep2D(vector<vector<HPRIntermediateStateVector> > & currentCells, vector<vector<HPRIntermediateStateVector> > & currentCellsWithBoundary,
                                        double cellSpacing, double timeStep, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static void computeXFORCETimeStep2D(vector<vector<HPRReducedStateVector> > & currentCells, vector<vector<HPRReducedStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                        double timeStep, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static void computeYFORCETimeStep2D(vector<vector<HPRStateVector> > & currentCells, vector<vector<HPRStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                        HPRMaterialParameters materialParameters);
    static void computeYFORCETimeStep2D(vector<vector<HPRIntermediateStateVector> > & currentCells, vector<vector<HPRIntermediateStateVector> > & currentCellsWithBoundary,
                                        double cellSpacing, double timeStep, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static void computeYFORCETimeStep2D(vector<vector<HPRReducedStateVector> > & currentCells, vector<vector<HPRReducedStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                        double timeStep, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static vector<HPRStateVector> solve(vector<HPRStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                        HPRMaterialParameters materialParameters);
    static vector<HPRMultiphysicsStateVector> solve(vector<HPRMultiphysicsStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                                    HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static vector<HPRIntermediateStateVector> solve(vector<HPRIntermediateStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                                    HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static vector<HPRReducedStateVector> solve(vector<HPRReducedStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                               HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static vector<vector<HPRStateVector> > solve2D(vector<vector<HPRStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                                   HPRMaterialParameters materialParameters);
    static vector<vector<HPRIntermediateStateVector> > solve2D(vector<vector<HPRIntermediateStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                               int subcyclingIterations, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
    static vector<vector<HPRReducedStateVector> > solve2D(vector<vector<HPRReducedStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                          int subcyclingIterations, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
};

#endif // HPRFIRSTORDERSOLVER_H

#ifndef FIRSTORDERSOLVER_H
#define FIRSTORDERSOLVER_H

#include "solvers.h"
#include "multiphysicssolvers.h"

class FirstOrderSolver
{
public:
    FirstOrderSolver();

    static vector<double> computeXLaxFriedrichsFlux(EulerStateVector leftStateVector, EulerStateVector rightStateVector, double cellSpacing, double timeStep,
                                                    EulerMaterialParameters materialParameters);
    static vector<double> computeXLaxFriedrichsFlux(EulerMultiphysicsStateVector leftStateVector, EulerMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                                    EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static vector<double> computeXLaxFriedrichsFlux(EulerReducedStateVector leftStateVector, EulerReducedStateVector rightStateVector, double cellSpacing, double timeStep,
                                                    EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> computeYLaxFriedrichsFlux(EulerStateVector topStateVector, EulerStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                    EulerMaterialParameters materialParameters);
    static vector<double> computeYLaxFriedrichsFlux(EulerMultiphysicsStateVector topStateVector, EulerMultiphysicsStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                    EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static vector<double> computeYLaxFriedrichsFlux(EulerReducedStateVector topStateVector, EulerReducedStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                    EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> computeLaxFriedrichsFlux(vector<double> leftConservedVariableVector, vector<double> rightConservedVariableVector, vector<double> leftFluxVector,
                                                   vector<double> rightFluxVector, double cellSpacing, double timeStep);

    static vector<double> computeXRichtmyerFlux(EulerStateVector leftStateVector, EulerStateVector rightStateVector, double cellSpacing, double timeStep,
                                                EulerMaterialParameters materialParameters);
    static vector<double> computeXRichtmyerFlux(EulerMultiphysicsStateVector leftStateVector, EulerMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                                EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static vector<double> computeXRichtmyerFlux(EulerReducedStateVector leftStateVector, EulerReducedStateVector rightStateVector, double cellSpacing, double timeStep,
                                                EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> computeYRichtmyerFlux(EulerStateVector topStateVector, EulerStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                EulerMaterialParameters materialParameters);
    static vector<double> computeYRichtmyerFlux(EulerMultiphysicsStateVector topStateVector, EulerMultiphysicsStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static vector<double> computeYRichtmyerFlux(EulerReducedStateVector topStateaVector, EulerReducedStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> computeRichtmyerFlux(vector<double> leftConservedVariableVector, vector<double> rightConservedVariableVector, vector<double> leftFluxVector,
                                               vector<double> rightFluxVector,
                                               double cellSpacing, double timeStep);

    static vector<double> computeXFORCEFlux(EulerStateVector leftStateVector, EulerStateVector rightStateVector, double cellSpacing, double timeStep,
                                            EulerMaterialParameters materialParameters);
    static vector<double> computeXFORCEFlux(EulerMultiphysicsStateVector leftStateVector, EulerMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                            EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static vector<double> computeXFORCEFlux(EulerReducedStateVector leftStateVector, EulerReducedStateVector rightStateVector, double cellSpacing, double timeStep,
                                            EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> computeYFORCEFlux(EulerStateVector topStateVector, EulerStateVector bottomStateVector, double cellSpacing, double timeStep,
                                            EulerMaterialParameters materialParameters);
    static vector<double> computeYFORCEFlux(EulerMultiphysicsStateVector topStateVector, EulerMultiphysicsStateVector bottomStateVector, double cellspacing, double timeStep,
                                            EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static vector<double> computeYFORCEFlux(EulerReducedStateVector topStateVector, EulerReducedStateVector bottomStateVector, double cellSpacing, double timeStep,
                                            EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> computeFORCEFlux(vector<double> laxFriedrichsFlux, vector<double> richtmyerFlux);

    static void computeFORCETimeStep(vector<EulerStateVector> & currentCells, vector<EulerStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                     EulerMaterialParameters materialParameters);
    static void computeFORCETimeStep(vector<EulerMultiphysicsStateVector> & currentCells, vector<EulerMultiphysicsStateVector> & currentCellsWithBoundary, double cellSpacing,
                                     double timeStep, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static void computeFORCETimeStep(vector<EulerReducedStateVector> & currentCells, vector<EulerReducedStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                     EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static void computeXFORCETimeStep2D(vector<vector<EulerStateVector> > & currentCells, vector<vector<EulerStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                        double timeStep, EulerMaterialParameters materialParameters);
    static void computeXFORCETimeStep2D(vector<vector<EulerMultiphysicsStateVector> > & currentCells, vector<vector<EulerMultiphysicsStateVector> > & currentCellsWithBoundary,
                                        double cellSpacing, double timeStep, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static void computeXFORCETimeStep2D(vector<vector<EulerReducedStateVector> > & currentCells, vector<vector<EulerReducedStateVector> > & currentCellsWithBoundary,
                                        double cellSpacing, double timeStep, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static void computeYFORCETimeStep2D(vector<vector<EulerStateVector> > & currentCells, vector<vector<EulerStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                        double timeStep, EulerMaterialParameters materialParameters);
    static void computeYFORCETimeStep2D(vector<vector<EulerMultiphysicsStateVector> > & currentCells, vector<vector<EulerMultiphysicsStateVector> > & currentCellsWithBoundary,
                                        double cellSpacing, double timeStep, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static void computeYFORCETimeStep2D(vector<vector<EulerReducedStateVector> > & currentCells, vector<vector<EulerReducedStateVector> > & currentCellsWithBoundary,
                                        double cellSpacing, double timeStep, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> computeFORCEUpdate(vector<double> conservedVariableVector, vector<double> leftFluxVector, vector<double> rightFluxVector, double cellSpacing,
                                             double timeStep);

    static vector<EulerStateVector> solve(vector<EulerStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                          EulerMaterialParameters materialParameters);
    static vector<EulerMultiphysicsStateVector> solve(vector<EulerMultiphysicsStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                      int subcyclingIterations, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static vector<EulerReducedStateVector> solve(vector<EulerReducedStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                 int subcyclingIterations, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<vector<EulerStateVector> > solve2D(vector<vector<EulerStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                     int subcyclingIterations, EulerMaterialParameters materialParameters);
    static vector<vector<EulerMultiphysicsStateVector> > solve2D(vector<vector<EulerMultiphysicsStateVector> > & initialCells, double cellSpacing, double CFLCoefficient,
                                                                 double finalTime, int subcyclingIterations, EulerMaterialParameters material1Parameters,
                                                                 EulerMaterialParameters material2Parameters);
    static vector<vector<EulerReducedStateVector> > solve2D(vector<vector<EulerReducedStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                            int subcyclingIterations, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
};

#endif // FIRSTORDERSOLVER_H

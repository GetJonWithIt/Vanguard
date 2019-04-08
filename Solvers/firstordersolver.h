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

    static vector<double> computeXLaxFriedrichsFlux(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, double cellSpacing, double timeStep,
                                                    HyperelasticMaterialParameters materialParameters);
    static vector<double> computeXLaxFriedrichsFlux(ElasticMultiphysicsStateVector leftStateVector, ElasticMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                                    HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static vector<double> computeXLaxFriedrichsFlux(ElasticReducedStateVector leftStateVector, ElasticReducedStateVector rightStateVector, double cellSpacing, double timeStep,
                                                    HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeYLaxFriedrichsFlux(EulerStateVector topStateVector, EulerStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                    EulerMaterialParameters materialParameters);
    static vector<double> computeYLaxFriedrichsFlux(EulerMultiphysicsStateVector topStateVector, EulerMultiphysicsStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                    EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> computeYLaxFriedrichsFlux(ElasticStateVector topStateVector, ElasticStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                    HyperelasticMaterialParameters materialParameters);
    static vector<double> computeYLaxFriedrichsFlux(ElasticReducedStateVector topStateVector, ElasticReducedStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                    HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeLaxFriedrichsFlux(vector<double> leftConservedVariableVector, vector<double> rightConservedVariableVector, vector<double> leftFluxVector,
                                                   vector<double> rightFluxVector, double cellSpacing, double timeStep);

    static vector<double> computeXRichtmyerFlux(EulerStateVector leftStateVector, EulerStateVector rightStateVector, double cellSpacing, double timeStep,
                                                EulerMaterialParameters materialParameters);
    static vector<double> computeXRichtmyerFlux(EulerMultiphysicsStateVector leftStateVector, EulerMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                                EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static vector<double> computeXRichtmyerFlux(EulerReducedStateVector leftStateVector, EulerReducedStateVector rightStateVector, double cellSpacing, double timeStep,
                                                EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> computeXRichtmyerFlux(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, double cellSpacing, double timeStep,
                                                HyperelasticMaterialParameters materialParameters);
    static vector<double> computeXRichtmyerFlux(ElasticMultiphysicsStateVector leftStateVector, ElasticMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                                HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static vector<double> computeXRichtmyerFlux(ElasticReducedStateVector leftStateVector, ElasticReducedStateVector rightStateVector, double cellSpacing, double timeStep,
                                                HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeYRichtmyerFlux(EulerStateVector topStateVector, EulerStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                EulerMaterialParameters materialParameters);
    static vector<double> computeYRichtmyerFlux(EulerMultiphysicsStateVector topStateVector, EulerMultiphysicsStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> computeYRichtmyerFlux(ElasticStateVector topStateVector, ElasticStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                HyperelasticMaterialParameters materialParameters);
    static vector<double> computeYRichtmyerFlux(ElasticReducedStateVector topStateVector, ElasticReducedStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeRichtmyerFlux(vector<double> leftConservedVariableVector, vector<double> rightConservedVariableVector, vector<double> leftFluxVector,
                                               vector<double> rightFluxVector,
                                               double cellSpacing, double timeStep);

    static vector<double> computeXFORCEFlux(EulerStateVector leftStateVector, EulerStateVector rightStateVector, double cellSpacing, double timeStep, EulerMaterialParameters materialParameters);
    static vector<double> computeXFORCEFlux(EulerMultiphysicsStateVector leftStateVector, EulerMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                            EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static vector<double> computeXFORCEFlux(EulerReducedStateVector leftStateVector, EulerReducedStateVector rightStateVector, double cellSpacing, double timeStep,
                                            EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> computeXFORCEFlux(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, double cellSpacing, double timeStep,
                                            HyperelasticMaterialParameters materialParameters);
    static vector<double> computeXFORCEFlux(ElasticMultiphysicsStateVector leftStateVector, ElasticMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                            HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static vector<double> computeXFORCEFlux(ElasticReducedStateVector leftStateVector, ElasticReducedStateVector rightStateVector, double cellSpacing, double timeStep,
                                            HyperelasticMaterialParameters material1Parametres, HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeYFORCEFlux(EulerStateVector topStateVector, EulerStateVector bottomStateVector, double cellSpacing, double timeStep, EulerMaterialParameters materialParameters);
    static vector<double> computeYFORCEFlux(EulerMultiphysicsStateVector topStateVector, EulerMultiphysicsStateVector bottomStateVector, double cellspacing, double timeStep,
                                            EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> computeYFORCEFlux(ElasticStateVector topStateVector, ElasticStateVector bottomStateVector, double cellSpacing, double timeStep,
                                            HyperelasticMaterialParameters materialParameters);
    static vector<double> computeYFORCEFlux(ElasticReducedStateVector topStateVector, ElasticReducedStateVector bottomStateVector, double cellSpacing, double timeStep,
                                            HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeFORCEFlux(vector<double> laxFriedrichsFlux, vector<double> richtmyerFlux);

    static void computeFORCETimeStep(vector<EulerStateVector> & currentCells, vector<EulerStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                     EulerMaterialParameters materialParameters);
    static void computeFORCETimeStep(vector<EulerMultiphysicsStateVector> & currentCells, vector<EulerMultiphysicsStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                     EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static void computeFORCETimeStep(vector<EulerReducedStateVector> & currentCells, vector<EulerReducedStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                     EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static void computeFORCETimeStep(vector<ElasticStateVector> & currentCells, vector<ElasticStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                     HyperelasticMaterialParameters materialParameters);
    static void computeFORCETimeStep(vector<ElasticMultiphysicsStateVector> & currentCells, vector<ElasticMultiphysicsStateVector> & currentCellsWithBoundary, double cellSpacing,
                                     double timeStep, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static void computeFORCETimeStep(vector<ElasticReducedStateVector> & currentCells, vector<ElasticReducedStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                     HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static void computeXFORCETimeStep2D(vector<vector<EulerStateVector> > & currentCells, vector<vector<EulerStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                        EulerMaterialParameters materialParameters);
    static void computeXFORCETimeStep2D(vector<vector<EulerMultiphysicsStateVector> > & currentCells, vector<vector<EulerMultiphysicsStateVector> > & currentCellsWithBoundary,
                                        double cellSpacing, double timeStep, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static void computeXFORCETimeStep2D(vector<vector<ElasticStateVector> > & currentCells, vector<vector<ElasticStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                        double timeStep, HyperelasticMaterialParameters materialParameters);
    static void computeXFORCETimeStep2D(vector<vector<ElasticReducedStateVector> > & currentCells, vector<vector<ElasticReducedStateVector> > & currentCellsWithBoundary,
                                        double cellSpacing, double timeStep, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static void computeYFORCETimeStep2D(vector<vector<EulerStateVector> > & currentCells, vector<vector<EulerStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                        EulerMaterialParameters materialParameters);
    static void computeYFORCETimeStep2D(vector<vector<EulerMultiphysicsStateVector> > & currentCells, vector<vector<EulerMultiphysicsStateVector> > & currentCellsWithBoundary,
                                        double cellSpacing, double timeStep, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static void computeYFORCETimeStep2D(vector<vector<ElasticStateVector> > & currentCells, vector<vector<ElasticStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                        double timeStep, HyperelasticMaterialParameters materialParameters);
    static void computeYFORCETimeStep2D(vector<vector<ElasticReducedStateVector> > & currentCells, vector<vector<ElasticReducedStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                        double timeStep, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeFORCEUpdate(vector<double> conservedVariableVector, vector<double> leftFluxVector, vector<double> rightFluxVector, double cellSpacing, double timeStep);

    static vector<EulerStateVector> solve(vector<EulerStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                          EulerMaterialParameters materialParameters);
    static vector<EulerMultiphysicsStateVector> solve(vector<EulerMultiphysicsStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                      int subcyclingIterations, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static vector<EulerReducedStateVector> solve(vector<EulerReducedStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                                 EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<ElasticStateVector> solve(vector<ElasticStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                            HyperelasticMaterialParameters materialParameters);
    static vector<ElasticMultiphysicsStateVector> solve(vector<ElasticMultiphysicsStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                        int subcyclingIterations, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static vector<ElasticReducedStateVector> solve(vector<ElasticReducedStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                                   HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<vector<EulerStateVector> > solve2D(vector<vector<EulerStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                                     EulerMaterialParameters materialParameters);
    static vector<vector<EulerMultiphysicsStateVector> > solve2D(vector<vector<EulerMultiphysicsStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                                 int subcyclingIterations, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<vector<ElasticStateVector> > solve2D(vector<vector<ElasticStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                                     HyperelasticMaterialParameters materialParameters);
    static vector<vector<ElasticReducedStateVector> > solve2D(vector<vector<ElasticReducedStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                              int subcyclingIterations, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
};

#endif // FIRSTORDERSOLVER_H

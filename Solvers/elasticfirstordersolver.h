#ifndef ELASTICFIRSTORDERSOLVER_H
#define ELASTICFIRSTORDERSOLVER_H

#include "elasticsolvers.h"
using namespace std;

class ElasticFirstOrderSolver
{
public:
    ElasticFirstOrderSolver();

    static vector<double> computeXLaxFriedrichsFlux(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, double cellSpacing, double timeStep,
                                                    HyperelasticMaterialParameters materialParameters);
    static vector<double> computeXLaxFriedrichsFlux(ElasticMultiphysicsStateVector leftStateVector, ElasticMultiphysicsStateVector rightStateVector, double cellSpacing,
                                                    double timeStep, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static vector<double> computeXLaxFriedrichsFlux(ElasticReducedStateVector leftStateVector, ElasticReducedStateVector rightStateVector, double cellSpacing, double timeStep,
                                                    HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeYLaxFriedrichsFlux(ElasticStateVector topStateVector, ElasticStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                    HyperelasticMaterialParameters materialParameters);
    static vector<double> computeYLaxFriedrichsFlux(ElasticReducedStateVector topStateVector, ElasticReducedStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                    HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeXRichtmyerFlux(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, double cellSpacing, double timeStep,
                                                HyperelasticMaterialParameters materialParameters);
    static vector<double> computeXRichtmyerFlux(ElasticMultiphysicsStateVector leftStateVector, ElasticMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                                HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static vector<double> computeXRichtmyerFlux(ElasticReducedStateVector leftStateVector, ElasticReducedStateVector rightStateVector, double cellSpacing, double timeStep,
                                                HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeYRichtmyerFlux(ElasticStateVector topStateVector, ElasticStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                HyperelasticMaterialParameters materialParameters);
    static vector<double> computeYRichtmyerFlux(ElasticReducedStateVector topStateVector, ElasticReducedStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeXFORCEFlux(ElasticStateVector leftStateVector, ElasticStateVector rightStateVector, double cellSpacing, double timeStep,
                                            HyperelasticMaterialParameters materialParameters);
    static vector<double> computeXFORCEFlux(ElasticMultiphysicsStateVector leftStateVector, ElasticMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep,
                                            HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static vector<double> computeXFORCEFlux(ElasticReducedStateVector leftStateVector, ElasticReducedStateVector rightStateVector, double cellSpacing, double timeStep,
                                            HyperelasticMaterialParameters material1Parametres, HyperelasticMaterialParameters material2Parameters);

    static vector<double> computeYFORCEFlux(ElasticStateVector topStateVector, ElasticStateVector bottomStateVector, double cellSpacing, double timeStep,
                                            HyperelasticMaterialParameters materialParameters);
    static vector<double> computeYFORCEFlux(ElasticReducedStateVector topStateVector, ElasticReducedStateVector bottomStateVector, double cellSpacing, double timeStep,
                                            HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static void computeFORCETimeStep(vector<ElasticStateVector> & currentCells, vector<ElasticStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                     HyperelasticMaterialParameters materialParameters);
    static void computeFORCETimeStep(vector<ElasticMultiphysicsStateVector> & currentCells, vector<ElasticMultiphysicsStateVector> & currentCellsWithBoundary, double cellSpacing,
                                     double timeStep, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);
    static void computeFORCETimeStep(vector<ElasticReducedStateVector> & currentCells, vector<ElasticReducedStateVector> & currentCellsWithBoundary, double cellSpacing,
                                     double timeStep, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static void computeXFORCETimeStep2D(vector<vector<ElasticStateVector> > & currentCells, vector<vector<ElasticStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                        double timeStep, HyperelasticMaterialParameters materialParameters);
    static void computeXFORCETimeStep2D(vector<vector<ElasticReducedStateVector> > & currentCells, vector<vector<ElasticReducedStateVector> > & currentCellsWithBoundary,
                                        double cellSpacing, double timeStep, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static void computeYFORCETimeStep2D(vector<vector<ElasticStateVector> > & currentCells, vector<vector<ElasticStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                        double timeStep, HyperelasticMaterialParameters materialParameters);
    static void computeYFORCETimeStep2D(vector<vector<ElasticReducedStateVector> > & currentCells, vector<vector<ElasticReducedStateVector> > & currentCellsWithBoundary,
                                        double cellSpacing, double timeStep, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<ElasticStateVector> solve(vector<ElasticStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                            HyperelasticMaterialParameters materialParameters);
    static vector<ElasticMultiphysicsStateVector> solve(vector<ElasticMultiphysicsStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                        int subcyclingIterations, HyperelasticMaterialParameters material1Parameters,
                                                        HyperelasticMaterialParameters material2Parameters);
    static vector<ElasticReducedStateVector> solve(vector<ElasticReducedStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                   int subcyclingIterations, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters);

    static vector<vector<ElasticStateVector> > solve2D(vector<vector<ElasticStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                       int subcyclingIterations, HyperelasticMaterialParameters materialParameters);
    static vector<vector<ElasticReducedStateVector> > solve2D(vector<vector<ElasticReducedStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime,
                                                              int subcyclingIterations, HyperelasticMaterialParameters material1Parameters,
                                                              HyperelasticMaterialParameters material2Parameters);
};

#endif // ELASTICFIRSTORDERSOLVER_H

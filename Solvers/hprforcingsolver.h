#ifndef HPRFORCINGSOLVER_H
#define HPRFORCINGSOLVER_H

#include "hprsolvers.h"
using namespace std;

class HPRForcingSolver
{
public:
    HPRForcingSolver();

    static vector<double> evolveConservedVariableVector(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector, vector<double> rightConservedVariableVector,
                                                        double cellSpacing, double timeStep, double bias, int slopeLimiter, HPRMaterialParameters materialParameters);
    static vector<double> evolveReducedConservedVariableVector(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector,
                                                               vector<double> rightConservedVariableVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                               HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static vector<double> evolveConservedVariableVector2D(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector, vector<double> rightConservedVariableVector,
                                                          vector<double> topConservedVariableVector, vector<double> bottomConservedVariableVector, double cellSpacing, double timeStep,
                                                          double bias, int slopeLimiter, HPRMaterialParameters materialParameters);
    static vector<double> evolveReducedConservedVariableVector2D(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector,
                                                                 vector<double> rightConservedVariableVector, vector<double> topConservedVariableVector,
                                                                 vector<double> bottomConservedVariableVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                 HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static vector<double> computeEvolvedConservedVariableVector(vector<double> middleConservedVariableVector, vector<double> firstStep, vector<double> secondStep, vector<double> thirdStep,
                                                                vector<double> fourthStep);

    static void computeRungeKuttaTimeStep(vector<HPRStateVector> & currentCells, vector<HPRStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep, double bias,
                                          int slopeLimiter, HPRMaterialParameters materialParameters);
    static void computeRungeKuttaTimeStep(vector<HPRReducedStateVector> & currentCells, vector<HPRReducedStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                          double bias, int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static void computeRungeKuttaTimeStep2D(vector<vector<HPRStateVector> > & currentCells, vector<vector<HPRStateVector> > & currentCellsWithBoundary, double cellSpacing, double timStep,
                                            double bias, int slopeLimiter, HPRMaterialParameters materialParameters);
    static void computeRungeKuttaTimeStep2D(vector<vector<HPRReducedStateVector> > & currentCells, vector<vector<HPRReducedStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                            double timeStep, double bias, int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);
};

#endif // HPRFORCINGSOLVER_H

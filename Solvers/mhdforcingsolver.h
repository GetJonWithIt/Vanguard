#ifndef MHDFORCINGSOLVER_H
#define MHDFORCINGSOLVER_H

#include "mhdsolvers.h"
#include "hprforcingsolver.h"
using namespace std;

class MHDForcingSolver
{
public:
    MHDForcingSolver();

    static vector<double> evolveConservedVariableVector(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector, vector<double> rightConservedVariableVector,
                                                        double cellSpacing, double timeStep, double bias, int slopeLimiter, MHDMaterialParameters materialParameters);
    static vector<double> evolveReducedConservedVariableVector(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector,
                                                               vector<double> rightConservedVariableVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                               MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static vector<double> evolveConservedVariableVector2D(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector, vector<double> rightConservedVariableVector,
                                                          vector<double> topConservedVariableVector, vector<double> bottomConservedVariableVector, double cellSpacing, double timeStep,
                                                          double bias, int slopeLimiter, MHDMaterialParameters materialParameters);
    static vector<double> evolveReducedConservedVariableVector2D(vector<double> leftConservedVariableVector, vector<double> middleConservedVariableVector,
                                                                 vector<double> rightConservedVariableVector, vector<double> topConservedVariableVector,
                                                                 vector<double> bottomConservedVariableVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                 MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static void computeRungeKuttaTimeStep(vector<MHDStateVector> & currentCells, vector<MHDStateVector> & currentCellsWithBoundary, double cellSpacing, double tiemStep, double bias,
                                          int slopeLimiter, MHDMaterialParameters materialParameters);
    static void computeRungeKuttaTimeStep(vector<MHDReducedStateVector> & currentCells, vector<MHDReducedStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                          double bias, int slopeLimiter, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static void computeRungeKuttaTimeStep2D(vector<vector<MHDStateVector> > & currentCells, vector<vector<MHDStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                            double bias, int slopeLimiter, MHDMaterialParameters materialParameters);
    static void computeRungeKuttaTimeStep2D(vector<vector<MHDReducedStateVector> > & currentCells, vector<vector<MHDReducedStateVector> > & currentCellsWithBoundary, double cellSpacing,
                                            double timeStep, double bias, int slopeLimiter, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
};

#endif // MHDFORCINGSOLVER_H

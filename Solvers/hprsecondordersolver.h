#ifndef HPRSECONDORDERSOLVER_H
#define HPRSECONDORDERSOLVER_H

#include "hprfirstordersolver.h"
using namespace std;

class HPRSecondOrderSolver
{
public:
    HPRSecondOrderSolver();

    static vector<double> computeXSLICFlux(HPRStateVector leftLeftStateVector, HPRStateVector leftStateVector, HPRStateVector rightStateVector, HPRStateVector rightRightStateVector,
                                           double cellSpacing, double timeStep, double bias, int slopeLimiter, HPRMaterialParameters materialParameters);
    static vector<double> computeXSLICFlux(HPRReducedStateVector leftLeftStateVector, HPRReducedStateVector leftStateVector, HPRReducedStateVector rightStateVector,
                                           HPRReducedStateVector rightRightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                           HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static vector<double> computeYSLICFlux(HPRStateVector topTopStateVector, HPRStateVector topStateVector, HPRStateVector bottomStateVector, HPRStateVector bottomBottomStateVector,
                                           double cellSpacing, double timeStep, double bias, int slopeLimiter, HPRMaterialParameters materialParameters);

    static void computeSLICTimeStep(vector<HPRStateVector> & currentCells, vector<HPRStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep, double bias,int slopeLimiter,
                                    HPRMaterialParameters materialParameters);
    static void computeSLICTimeStep(vector<HPRReducedStateVector> & currentCells, vector<HPRReducedStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep, double bias,
                                    int slopeLimiter, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static void computeXSLICTimeStep2D(vector<vector<HPRStateVector> > & currentCells, vector<vector<HPRStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep, double bias,
                                       int slopeLimiter, HPRMaterialParameters materialParameters);

    static void computeYSLICTimeStep2D(vector<vector<HPRStateVector> > & currentCells, vector<vector<HPRStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep, double bias,
                                       int slopeLimiter, HPRMaterialParameters materialParameters);

    static vector<HPRStateVector> solve(vector<HPRStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                        int subcyclingIterations, HPRMaterialParameters materialParameters);
    static vector<HPRReducedStateVector> solve(vector<HPRReducedStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                               int subcyclingIterations, int reinitialisationFrequency, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters);

    static vector<vector<HPRStateVector> > solve2D(vector<vector<HPRStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                                   int subcyclingIterations, HPRMaterialParameters materialParameters);
};

#endif // HPRSECONDORDERSOLVER_H

#ifndef HPRFIRSTORDERSOLVER_H
#define HPRFIRSTORDERSOLVER_H

#include "hprforcingsolver.h"
using namespace std;

class HPRFirstOrderSolver
{
public:
    HPRFirstOrderSolver();

    static vector<double> computeXLaxFriedrichsFlux(HPRStateVector leftStateVector, HPRStateVector rightStateVector, double cellSpacing, double timeStep, HPRMaterialParameters materialParameters);

    static vector<double> computeYLaxFriedrichsFlux(HPRStateVector topStateVector, HPRStateVector bottomStateVector, double cellSpacing, double timeStep, HPRMaterialParameters materialParameters);

    static vector<double> computeXRichtmyerFlux(HPRStateVector leftStateVector, HPRStateVector rightStateVector, double cellSpacing, double timeStep, HPRMaterialParameters materialParameters);

    static vector<double> computeYRichtmyerFlux(HPRStateVector topStateVector, HPRStateVector bottomStateVector, double cellSpacing, double timeStep, HPRMaterialParameters materialParameters);

    static vector<double> computeXFORCEFlux(HPRStateVector leftStateVector, HPRStateVector rightStateVector, double cellSpacing, double timeStep, HPRMaterialParameters materialParameters);

    static vector<double> computeYFORCEFlux(HPRStateVector topStateVector, HPRStateVector bottomStateVector, double cellSpacing, double timeStep, HPRMaterialParameters materialParameters);

    static void computeFORCETimeStep(vector<HPRStateVector> & currentCells, vector<HPRStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                     HPRMaterialParameters materialParameters);

    static void computeXFORCETimeStep2D(vector<vector<HPRStateVector> > & currentCells, vector<vector<HPRStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                        HPRMaterialParameters materialParameters);

    static void computeYFORCETimeStep2D(vector<vector<HPRStateVector> > & currentCells, vector<vector<HPRStateVector> > & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                        HPRMaterialParameters materialParameters);

    static vector<HPRStateVector> solve(vector<HPRStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                        HPRMaterialParameters materialParameters);

    static vector<vector<HPRStateVector> > solve2D(vector<vector<HPRStateVector> > & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                                   HPRMaterialParameters materialParameters);
};

#endif // HPRFIRSTORDERSOLVER_H

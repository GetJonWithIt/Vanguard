#ifndef MHDFIRSTORDERSOLVER_H
#define MHDFIRSTORDERSOLVER_H

#include "mhdsolvers.h"
using namespace std;

class MHDFirstOrderSolver
{
public:
    MHDFirstOrderSolver();

    static vector<double> computeXLaxFriedrichsFlux(MHDStateVector leftStateVector, MHDStateVector rightStateVector, double cellSpacing, double timeStep, MHDMaterialParameters materialParameters);

    static vector<double> computeXRichtmyerFlux(MHDStateVector leftStateVector, MHDStateVector rightStateVector, double cellSpacing, double timeStep, MHDMaterialParameters materialParameters);

    static vector<double> computeXFORCEFlux(MHDStateVector leftStateVector, MHDStateVector rightStateVector, double cellSpacing, double timeStep, MHDMaterialParameters materialParameters);

    static void computeFORCETimeStep(vector<MHDStateVector> & currentCells, vector<MHDStateVector> & currentCellsWithBoundary, double cellSpacing, double timeStep,
                                     MHDMaterialParameters materialParameters);

    static vector<MHDStateVector> solve(vector<MHDStateVector> & initialCells, double cellSpacing, double CFLCoefficient, double finalTime, int subcyclingIterations,
                                        MHDMaterialParameters materialParameters);
};

#endif // MHDFIRSTORDERSOLVER_H

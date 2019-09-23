#ifndef MHDRGFMSOLVER_H
#define MHDRGFMSOLVER_H

#include "mhdmultimaterialsystem.h"
#include "rgfmsolver.h"
using namespace std;

class MHDRGFMSolver
{
public:
    MHDRGFMSolver();

    static MHDMultimaterialSystem applyRGFMBoundaryConditions(MHDMultimaterialSystem multimaterialSystem, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static vector<double> updateLevelSetFunction(vector<double> levelSetFunction, double cellSpacing, double timeStep, vector<MHDStateVector> material1Cells,
                                                 vector<MHDStateVector> material2Cells);

    static MHDMultimaterialSystem solve(MHDMultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                        int subcyclingIterations, int reinitialisationFrequency, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
};

#endif // MHDRGFMSOLVER_H

#ifndef RGFMSOLVER_H
#define RGFMSOLVER_H

#include "multimaterialsystem.h"
using namespace std;

class RGFMSolver
{
public:
    RGFMSolver();

    static vector<double> insertBoundaryCells(vector<double> & levelSetFunction, int boundarySize);
    static vector<vector<double> > insertBoundaryCells2D(vector<vector<double> > & levelSetFunction, int boundarySize);

    static MultimaterialSystem applyRGFMBoundaryConditionsExact(MultimaterialSystem multimaterialSystem, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static MultimaterialSystem applyRGFMBoundaryConditionsExact2D(MultimaterialSystem multimaterialSystem, EulerMaterialParameters material1Parameters,
                                                                  EulerMaterialParameters material2Parameters);

    static MultimaterialSystem applyRGFMBoundaryConditionsHLLC(MultimaterialSystem multimaterialSystem, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static MultimaterialSystem applyRGFMBoundaryConditionsHLLC2D(MultimaterialSystem multimaterialSystem, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static vector<double> updateLevelSetFunction(vector<double> levelSetFunction, double cellSpacing, double timeStep, vector<EulerStateVector> material1Cells,
                                                 vector<EulerStateVector> material2Cells);

    static vector<vector<double> > update2DLevelSetFunctionX(vector<vector<double> > levelSetFunction, double cellSpacing, double timeStep, vector<vector<EulerStateVector> > material1Cells,
                                                             vector<vector<EulerStateVector> > material2Cells);
    static vector<vector<double> > update2DLevelSetFunctionY(vector<vector<double> > levelSetFunction, double cellSpacing, double timeStep, vector<vector<EulerStateVector> > material1Cells,
                                                             vector<vector<EulerStateVector> > material2Cells);

    static MultimaterialSystem solveExact(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                          int subcyclingIterations, int reinitialisationFrequency, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static MultimaterialSystem solveExact2D(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                            int subcyclingIterations, int reinitialisationFrequency, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);

    static MultimaterialSystem solveHLLC(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                         int subcyclingIterations, int reinitialisationFrequency, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static MultimaterialSystem solveHLLC2D(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                           int subcyclingIterations, int reinitialisationFrequency, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
};

#endif // RGFMSOLVER_H

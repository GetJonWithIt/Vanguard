#ifndef EXACTSOLVER_H
#define EXACTSOLVER_H

#include "Solvers/secondordersolver.h"
using namespace std;

class ExactSolver
{
public:
    ExactSolver();

    static double computeACoefficient(double density, double adiabaticIndex);
    static double computeBCoefficient(double pressure, double adiabaticIndex, double stiffeningParameter);

    static double computeStarRegionSoundSpeed(double starRegionPressure, EulerStateVector stateVector, EulerMaterialParameters materialParameters);

    static vector<double> computeRarefactionXPrimitiveVariableVector(double waveSpeed, double soundSpeed, EulerStateVector stateVector, EulerMaterialParameters materialParameters);
    static vector<double> computeRarefactionYPrimitiveVariableVector(double waveSpeed, double soundSpeed, EulerStateVector stateVector, EulerMaterialParameters materialParameters);

    static double computeStarRegionShockDensity(double starRegionPressure, EulerStateVector stateVector, EulerMaterialParameters materialParameters);
    static double computeStarRegionRarefactionDensity(double starRegionPressure, EulerStateVector stateVector, EulerMaterialParameters materialParameters);

    static double computeStarRegionXVelocity(double starRegionPressure, EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                             EulerMaterialParameters material2Parameters);
    static double computeStarRegionYVelocity(double starRegionPressure, EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                             EulerMaterialParameters material2Parameters);

    static double computeStarRegionXPressure(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                             EulerMaterialParameters material2Parameters);
    static double computeStarRegionYPressure(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                             EulerMaterialParameters material2Parameters);

    static double computeWaveJumpFunctionComponent(double newPressure, EulerStateVector stateVector, EulerMaterialParameters materialParameters);

    static double computeXWaveJumpFunction(double newPressure, EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                           EulerMaterialParameters material2Parameters);
    static double computeYWaveJumpFunction(double newPressure, EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                           EulerMaterialParameters material2Parameters);

    static double computeWaveJumpFunctionComponentDerivative(double newPressure, EulerStateVector stateVector, EulerMaterialParameters materialParameters);
    static double computeWaveJumpFunctionDerivative(double newPressure, EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                                    EulerMaterialParameters material2Parameters);

    static double computePressureChange(double oldPressure, double newPressure);

    static double computeShockSpeed(double starRegionPressure, EulerStateVector stateVector, EulerMaterialParameters materialParameters);

    static double computeStarRegionDensity(double starRegionPressure, EulerStateVector stateVector, EulerMaterialParameters materialParameters);
    static EulerStateVector solveX(double position, double time, double interfaceLocation, EulerStateVector leftStateVector, EulerStateVector rightStateVector,
                                   EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
    static EulerStateVector solveY(double position, double time, double interfaceLocation, EulerStateVector topStateVector, EulerStateVector bottomStateVector,
                                   EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters);
};

#endif // EXACTSOLVER_H

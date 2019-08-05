#ifndef HLLCSOLVER_H
#define HLLCSOLVER_H

#include "Solvers/secondordersolver.h"
using namespace std;

class HLLCSolver
{
public:
    HLLCSolver();

    static double computeStarRegionDensity(double density, double waveSpeed, double velocity, double starRegionWaveSpeed);

    static double computeLeftStarRegionDensity(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                               EulerMaterialParameters material2Parameters);
    static double computeRightStarRegionDensity(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                                EulerMaterialParameters material2Parameters);

    static double computeTopStarRegionDensity(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                              EulerMaterialParameters material2Parameters);
    static double computeBottomStarRegionDensity(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                                 EulerMaterialParameters material2Parameters);

    static vector<double> computeStarRegionXConservedVariableVector(EulerStateVector stateVector, double waveSpeed, double starRegionWaveSpeed, EulerMaterialParameters materialParameters);
    static vector<double> computeStarRegionYConservedVariableVector(EulerStateVector stateVector, double waveSpeed, double starRegionWaveSpeed, EulerMaterialParameters materialParameters);

    static EulerStateVector solveX(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                   EulerMaterialParameters material2Parameters);
    static EulerStateVector solveY(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                   EulerMaterialParameters material2Parameters);

    static vector<double> computeStarRegionFlux(EulerStateVector stateVector, double waveSpeed, double starRegionWaveSpeed, EulerMaterialParameters materialParameters);
    static vector<double> computeHLLCFlux(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                          EulerMaterialParameters material2Parameters);

    static double computePVRSStarRegionXVelocity(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                                 EulerMaterialParameters material2Parameters);
    static double computePVRSStarRegionXPressure(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                                 EulerMaterialParameters material2Parameters);

    static double computePVRSStarRegionYVelocity(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                                 EulerMaterialParameters material2Parameters);
    static double computePVRSStarRegionYPressure(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                                 EulerMaterialParameters material2Parameters);

    static double computeTwoShockStarRegionXVelocity(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                                     EulerMaterialParameters material2Parameters);
    static double computeTwoShockStarRegionXPressure(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                                     EulerMaterialParameters material2Parameters);

    static double computeTwoShockStarRegionYVelocity(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                                     EulerMaterialParameters material2Parameters);
    static double computeTwoShockStarRegionYPressure(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                                     EulerMaterialParameters material2Parameters);

    static double computeTwoRarefactionStarRegionXVelocity(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                                           EulerMaterialParameters material2Parameters);
    static double computeTwoRarefactionStarRegionXPressure(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                                           EulerMaterialParameters material2Parameters);

    static double computeTwoRarefactionStarRegionYVelocity(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                                           EulerMaterialParameters material2Parameters);
    static double computeTwoRarefactionStarRegionYPressure(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                                           EulerMaterialParameters material2Parameters);

    static double computeStarRegionXVelocity(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                             EulerMaterialParameters material2Parameters);
    static double computeStarRegionXPressure(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                             EulerMaterialParameters material2Parameters);

    static double computeStarRegionYVelocity(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                             EulerMaterialParameters material2Parameters);
    static double computeStarRegionYPressure(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                             EulerMaterialParameters material2Parameters);

    static double computeWaveSpeedWeighting(double starRegionPressure, EulerStateVector stateVector, EulerMaterialParameters materialParameters);

    static double computeLeftWaveSpeed(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                       EulerMaterialParameters material2Parameters);
    static double computeRightWaveSpeed(EulerStateVector leftStateVector, EulerStateVector rightStateVector, EulerMaterialParameters material1Parameters,
                                        EulerMaterialParameters material2Parameters);

    static double computeTopWaveSpeed(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                      EulerMaterialParameters material2Parameters);
    static double computeBottomWaveSpeed(EulerStateVector topStateVector, EulerStateVector bottomStateVector, EulerMaterialParameters material1Parameters,
                                         EulerMaterialParameters material2Parameters);
};

#endif // HLLCSOLVER_H

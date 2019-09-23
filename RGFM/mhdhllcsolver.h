#ifndef MHDHLLCSOLVER_H
#define MHDHLLCSOLVER_H

#include "Solvers/mhdsecondordersolver.h"
using namespace std;

class MHDHLLCSolver
{
public:
    MHDHLLCSolver();

    static double computeStarRegionDensity(double density, double waveSpeed, double velocity, double starRegionWaveSpeed);

    static double computeLeftStarRegionDensity(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                               MHDMaterialParameters material2Parameters);
    static double computeRightStarRegionDensity(MHDStateVector leftStatevector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                MHDMaterialParameters material2Parameters);

    static vector<double> computeLeftStarRegionConservedVariableVector(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                                       MHDMaterialParameters material2Parameters);
    static vector<double> computeRightStarRegionConservedVariableVector(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                                        MHDMaterialParameters material2Parameters);

    static MHDStateVector solveX(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeTotalPressure(MHDStateVector stateVector);

    static double computeStarRegionXVelocity(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                             MHDMaterialParameters material2Parameters);

    static double computeLeftStarRegionTotalPressure(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                     MHDMaterialParameters material2Parameters);
    static double computeRightStarRegionTotalPressure(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                      MHDMaterialParameters material2Parameters);

    static double computeLeftStarRegionYMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                 MHDMaterialParameters material2Parameters);
    static double computeRightStarRegionYMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                  MHDMaterialParameters material2Parameters);

    static double computeLeftStarRegionZMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                 MHDMaterialParameters material2Parameters);
    static double computeRightStarRegionZMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                  MHDMaterialParameters material2Parameters);

    static double computeXStarRegionXMagneticField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                   MHDMaterialParameters material2Parameters);
    static double computeXStarRegionYMagneticField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                   MHDMaterialParameters material2Parameters);
    static double computeXStarRegionZMagneticField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                   MHDMaterialParameters material2Parameters);

    static double computeStarRegionXMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeStarRegionYMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeStarRegionZMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeXStarRegionDensity(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeXHLLDensity(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeXHLLXMagneticField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeXHLLYMagneticField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeXHLLZMagneticField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeHLLXMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeHLLYMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeHLLZMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeLeftStarRegionTotalEnergy(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                   MHDMaterialParameters material2Parameters);
    static double computeRightStarRegionTotalEnergy(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                    MHDMaterialParameters material2Parameters);

    static double computeLeftWaveSpeed(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeRightWaveSpeed(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
};

#endif // MHDHLLCSOLVER_H

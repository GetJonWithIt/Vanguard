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

    static double computeTopStarRegionDensity(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                              MHDMaterialParameters material2Parameters);
    static double computeBottomStarRegionDensity(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                 MHDMaterialParameters material2Parameters);

    static vector<double> computeLeftStarRegionConservedVariableVector(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                                       MHDMaterialParameters material2Parameters);
    static vector<double> computeRightStarRegionConservedVariableVector(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                                        MHDMaterialParameters material2Parameters);

    static vector<double> computeTopStarRegionConservedVariableVector(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                                      MHDMaterialParameters material2Parameters);
    static vector<double> computeBottomStarRegionConservedVariableVector(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                                         MHDMaterialParameters material2Parameters);

    static MHDStateVector solveX(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static MHDStateVector solveY(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeTotalPressure(MHDStateVector stateVector);

    static double computeStarRegionXVelocity(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                             MHDMaterialParameters material2Parameters);
    static double computeStarRegionYVelocity(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                             MHDMaterialParameters material2Parameters);

    static double computeLeftStarRegionTotalPressure(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                     MHDMaterialParameters material2Parameters);
    static double computeRightStarRegionTotalPressure(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                      MHDMaterialParameters material2Parameters);

    static double computeTopStarRegionTotalPressure(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                    MHDMaterialParameters material2Parameters);
    static double computeBottomStarRegionTotalPressure(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                       MHDMaterialParameters material2Parameters);

    static double computeLeftStarRegionYMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                 MHDMaterialParameters material2Parameters);
    static double computeRightStarRegionYMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                  MHDMaterialParameters material2Parameters);

    static double computeLeftStarRegionZMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                 MHDMaterialParameters material2Parameters);
    static double computeRightStarRegionZMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                  MHDMaterialParameters material2Parameters);

    static double computeTopStarRegionXMomentum(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                MHDMaterialParameters material2Parameters);
    static double computeBottomStarRegionXMomentum(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                   MHDMaterialParameters material2Parameters);

    static double computeTopStarRegionZMomentum(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                MHDMaterialParameters material2Parameters);
    static double computeBottomStarRegionZMomentum(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                   MHDMaterialParameters material2Parameters);

    static double computeXStarRegionXMagneticField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                   MHDMaterialParameters material2Parameters);
    static double computeXStarRegionYMagneticField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                   MHDMaterialParameters material2Parameters);
    static double computeXStarRegionZMagneticField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                   MHDMaterialParameters material2Parameters);

    static double computeYStarRegionXMagneticField(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                   MHDMaterialParameters material2Parameters);
    static double computeYStarRegionYMagneticField(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                   MHDMaterialParameters material2Parameters);
    static double computeYStarRegionZMagneticField(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                   MHDMaterialParameters material2Parameters);

    static double computeXStarRegionAuxiliaryField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                   MHDMaterialParameters material2Parameters);
    static double computeYStarRegionAuxiliaryField(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                   MHDMaterialParameters material2Parameters);

    static double computeXStarRegionXMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeXStarRegionYMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeXStarRegionZMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeYStarRegionXMomentum(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeYStarRegionYMomentum(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeYStarRegionZMomentum(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeXStarRegionDensity(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeYStarRegionDensity(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeXHLLDensity(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeYHLLDensity(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeXHLLXMagneticField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeXHLLYMagneticField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeXHLLZMagneticField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeYHLLXMagneticField(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeYHLLYMagneticField(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeYHLLZMagneticField(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeXHLLAuxiliaryField(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeYHLLAuxiliaryField(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeXHLLXMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeXHLLYMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeXHLLZMomentum(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeYHLLXMomentum(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeYHLLYMomentum(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeYHLLZMomentum(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeLeftStarRegionTotalEnergy(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                   MHDMaterialParameters material2Parameters);
    static double computeRightStarRegionTotalEnergy(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters,
                                                    MHDMaterialParameters material2Parameters);

    static double computeTopStarRegionTotalEnergy(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                  MHDMaterialParameters material2Parameters);
    static double computeBottomStarRegionTotalEnergy(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters,
                                                     MHDMaterialParameters material2Parameters);

    static double computeLeftWaveSpeed(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeRightWaveSpeed(MHDStateVector leftStateVector, MHDStateVector rightStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);

    static double computeTopWaveSpeed(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
    static double computeBottomWaveSpeed(MHDStateVector topStateVector, MHDStateVector bottomStateVector, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters);
};

#endif // MHDHLLCSOLVER_H

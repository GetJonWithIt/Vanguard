#include "mhdsolvers.h"

MHDSolvers::MHDSolvers()
{
}

vector<MHDStateVector> MHDSolvers::insertBoundaryCells(vector<MHDStateVector> & currentCells, int boundarySize)
{
    int cellCount = currentCells.size();
    vector<MHDStateVector> currentCellsWithBoundary(cellCount + (2 * boundarySize));

    if (boundarySize == 1)
    {
        currentCellsWithBoundary[0] = currentCells[0];
        currentCellsWithBoundary[cellCount + 1] = currentCells[cellCount - 1];
    }
    else if (boundarySize == 2)
    {
        currentCellsWithBoundary[0] = currentCells[1];
        currentCellsWithBoundary[1] = currentCells[0];
        currentCellsWithBoundary[cellCount + 2] = currentCells[cellCount - 1];
        currentCellsWithBoundary[cellCount + 3] = currentCells[cellCount - 2];
    }

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        currentCellsWithBoundary[i + boundarySize] = currentCells[i];
    }

    return currentCellsWithBoundary;
}

double MHDSolvers::computeMaximumWaveSpeed(vector<MHDStateVector> & currentCells, MHDMaterialParameters materialParameters)
{
    double maximumWaveSpeed = 0.0;
    int cellCount = currentCells.size();

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        double waveSpeed = abs(currentCells[i].getXVelocity()) + max(max(abs(currentCells[i].computeFastMagnetoAcousticSpeed(materialParameters)),
                                                                         abs(currentCells[i].computeAlfvenWaveSpeed())),
                                                                     abs(currentCells[i].computeSlowMagnetoAcousticSpeed(materialParameters)));

        if (waveSpeed > maximumWaveSpeed)
        {
            maximumWaveSpeed = waveSpeed;
        }
    }

    return maximumWaveSpeed;
}

double MHDSolvers::computeStableTimeStep(vector<MHDStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                         MHDMaterialParameters materialParameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed(currentCells, materialParameters));

    return Solvers::computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

MHDStateVector MHDSolvers::evolveStateByHalfXTimeStep(MHDStateVector leftStateVector, MHDStateVector middleStateVector, MHDStateVector rightStateVector, double cellSpacing, double timeStep,
                                                      double bias, int slopeLimiter, int side, MHDMaterialParameters materialParameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, slopeLimiter, materialParameters);
    vector<double> leftExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(materialParameters), VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(materialParameters), VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> leftFluxVector = MHDStateVector::computeXFluxVector(leftExtrapolatedValue, materialParameters);
    vector<double> rightFluxVector = MHDStateVector::computeXFluxVector(rightExtrapolatedValue, materialParameters);
    vector<double> evolutionVector = Solvers::computeEvolutionVector(leftFluxVector, rightFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStep(leftExtrapolatedValue, rightExtrapolatedValue, evolutionVector, side, materialParameters);
}

MHDStateVector MHDSolvers::evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                     MHDMaterialParameters materialParameters)
{
    MHDStateVector evolvedStateVector;

    if (side == 0)
    {
        evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(leftExtrapolatedValue, evolutionVector), materialParameters);
    }
    else
    {
        evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(rightExtrapolatedValue, evolutionVector), materialParameters);
    }

    return evolvedStateVector;
}

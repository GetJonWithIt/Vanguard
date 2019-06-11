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

vector<MHDReducedStateVector> MHDSolvers::insertBoundaryCells(vector<MHDReducedStateVector> & currentCells, int boundarySize)
{
    int cellCount = currentCells.size();
    vector<MHDReducedStateVector> currentCellsWithBoundary(cellCount + (2 * boundarySize));

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

vector<vector<MHDStateVector> > MHDSolvers::insertBoundaryCells2D(vector<vector<MHDStateVector> > & currentCells, int boundarySize)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();
    vector<vector<MHDStateVector> > currentCellsWithBoundary(rowCount + (2 * boundarySize), vector<MHDStateVector>(columnCount + (2 * boundarySize)));

    if (boundarySize == 1)
    {
#pragma omp parallel for
        for (int i = 0; i < rowCount; i++)
        {
            currentCellsWithBoundary[i + 1][0] = currentCells[i][0];
            currentCellsWithBoundary[i + 1][columnCount + 1] = currentCells[i][columnCount - 1];
        }

#pragma omp parallel for
        for (int i = 0; i < columnCount; i++)
        {
            currentCellsWithBoundary[0][i + 1] = currentCells[0][i];
            currentCellsWithBoundary[rowCount + 1][i + 1] = currentCells[rowCount - 1][i];
        }
    }
    else if (boundarySize == 2)
    {
#pragma omp parallel for
        for (int i = 0; i < rowCount; i++)
        {
            currentCellsWithBoundary[i + 2][0] = currentCells[i][1];
            currentCellsWithBoundary[i + 2][1] = currentCells[i][0];

            currentCellsWithBoundary[i + 2][columnCount + 2] = currentCells[i][columnCount - 1];
            currentCellsWithBoundary[i + 2][columnCount + 3] = currentCells[i][columnCount - 2];
        }

#pragma omp parallel for
        for (int i = 0; i < columnCount; i++)
        {
            currentCellsWithBoundary[0][i + 2] = currentCells[1][i];
            currentCellsWithBoundary[1][i + 2] = currentCells[0][i];

            currentCellsWithBoundary[rowCount + 2][i + 2] = currentCells[rowCount - 1][i];
            currentCellsWithBoundary[rowCount + 3][i + 2] = currentCells[rowCount - 2][i];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            currentCellsWithBoundary[i + boundarySize][j + boundarySize] = currentCells[i][j];
        }
    }

    return currentCellsWithBoundary;
}

vector<vector<MHDReducedStateVector> > MHDSolvers::insertBoundaryCells2D(vector<vector<MHDReducedStateVector> > & currentCells, int boundarySize)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();
    vector<vector<MHDReducedStateVector> > currentCellsWithBoundary(rowCount + (2 * boundarySize), vector<MHDReducedStateVector>(columnCount + (2 * boundarySize)));

    if (boundarySize == 1)
    {
#pragma omp parallel for
        for (int i = 0; i < rowCount; i++)
        {
            currentCellsWithBoundary[i + 1][0] = currentCells[i][0];
            currentCellsWithBoundary[i + 1][columnCount + 1] = currentCells[i][columnCount - 1];
        }

#pragma omp parallel for
        for (int i = 0; i < columnCount; i++)
        {
            currentCellsWithBoundary[0][i + 1] = currentCells[0][i];
            currentCellsWithBoundary[rowCount + 1][i + 1] = currentCells[rowCount -1][i];
        }
    }
    else if (boundarySize == 2)
    {
#pragma omp parallel for
        for (int i = 0; i < rowCount; i++)
        {
            currentCellsWithBoundary[i + 2][0] = currentCells[i][1];
            currentCellsWithBoundary[i + 2][1] = currentCells[i][0];

            currentCellsWithBoundary[i + 2][columnCount + 2] = currentCells[i][columnCount - 1];
            currentCellsWithBoundary[i + 2][columnCount + 3] = currentCells[i][columnCount - 2];
        }

#pragma omp parallel for
        for (int i = 0; i < columnCount; i++)
        {
            currentCellsWithBoundary[0][i + 2] = currentCells[1][i];
            currentCellsWithBoundary[1][i + 2] = currentCells[0][i];

            currentCellsWithBoundary[rowCount + 2][i + 2] = currentCells[rowCount - 1][i];
            currentCellsWithBoundary[rowCount + 3][i + 2] = currentCells[rowCount - 2][i];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            currentCellsWithBoundary[i + boundarySize][j + boundarySize] = currentCells[i][j];
        }
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
        double waveSpeed = abs(currentCells[i].getXVelocity()) + abs(currentCells[i].computeXFastMagnetoAcousticSpeed(materialParameters));

        if (waveSpeed > maximumWaveSpeed)
        {
            maximumWaveSpeed = waveSpeed;
        }
    }

    return maximumWaveSpeed;
}

double MHDSolvers::computeMaximumWaveSpeed(vector<MHDReducedStateVector> & currentCells, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double maximumWaveSpeed = 0.0;
    int cellCount = currentCells.size();

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        double waveSpeed = abs(currentCells[i].getInterfaceXVelocity()) + max(abs(currentCells[i].computeMaterial1XFastMagnetoAcousticSpeed(material1Parameters)),
                                                                            abs(currentCells[i].computeMaterial2XFastMagnetoAcousticSpeed(material2Parameters)));

        if (waveSpeed > maximumWaveSpeed)
        {
            maximumWaveSpeed = waveSpeed;
        }
    }

    return maximumWaveSpeed;
}

double MHDSolvers::computeMaximumWaveSpeed2D(vector<vector<MHDStateVector> > & currentCells, MHDMaterialParameters materialParameters)
{
    double maximumWaveSpeed = 0.0;
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            double waveSpeed = max(abs(currentCells[i][j].getXVelocity()), abs(currentCells[i][j].getYVelocity())) +
                    max(abs(currentCells[i][j].computeXFastMagnetoAcousticSpeed(materialParameters)), abs(currentCells[i][j].computeYFastMagnetoAcousticSpeed(materialParameters)));

            if (waveSpeed > maximumWaveSpeed)
            {
                maximumWaveSpeed = waveSpeed;
            }
        }
    }

    return maximumWaveSpeed;
}

double MHDSolvers::computeMaximumWaveSpeed2D(vector<vector<MHDReducedStateVector> > & currentCells, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double maximumWaveSpeed = 0.0;
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            double waveSpeed = max(abs(currentCells[i][j].getInterfaceXVelocity()), abs(currentCells[i][j].getInterfaceYVelocity())) +
                    max(max(abs(currentCells[i][j].computeMaterial1XFastMagnetoAcousticSpeed(material1Parameters)),
                            abs(currentCells[i][j].computeMaterial1YFastMagnetoAcousticSpeed(material1Parameters))),
                        max(abs(currentCells[i][j].computeMaterial2XFastMagnetoAcousticSpeed(material2Parameters)),
                            abs(currentCells[i][j].computeMaterial2YFastMagnetoAcousticSpeed(material2Parameters))));

            if (waveSpeed > maximumWaveSpeed)
            {
                maximumWaveSpeed = waveSpeed;
            }
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

double MHDSolvers::computeStableTimeStep(vector<MHDReducedStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                         int currentIteration, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed(currentCells, material1Parameters, material2Parameters));

    return Solvers::computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

double MHDSolvers::computeStableTimeStep2D(vector<vector<MHDStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                           int currentIteration, MHDMaterialParameters materialParameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed2D(currentCells, materialParameters));

    return Solvers::computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

double MHDSolvers::computeStableTimeStep2D(vector<vector<MHDReducedStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                           int currentIteration, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed2D(currentCells, material1Parameters, material2Parameters));

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

MHDReducedStateVector MHDSolvers::evolveStateByHalfXTimeStep(MHDReducedStateVector leftStateVector, MHDReducedStateVector middleStateVector, MHDReducedStateVector rightStateVector,
                                                             double cellSpacing, double timeStep, double bias, int slopeLimiter, int side, MHDMaterialParameters material1Parameters,
                                                             MHDMaterialParameters material2Parameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, slopeLimiter, material1Parameters, material2Parameters);
    vector<double> leftExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                          VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                      VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> leftFluxVector = MHDReducedStateVector::computeXFluxVector(leftExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> rightFluxVector = MHDReducedStateVector::computeXFluxVector(rightExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> evolutionVector = Solvers::computeEvolutionVector(leftFluxVector, rightFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStepReduced(leftExtrapolatedValue, rightExtrapolatedValue, evolutionVector, side, material1Parameters, material2Parameters);
}

MHDStateVector MHDSolvers::evolveStateByHalfYTimeStep(MHDStateVector topStateVector, MHDStateVector middleStateVector, MHDStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                      double bias, int slopeLimiter, int side, MHDMaterialParameters materialParameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(topStateVector, middleStateVector, bottomStateVector, bias, slopeLimiter, materialParameters);
    vector<double> topExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(materialParameters), VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> bottomExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(materialParameters), VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> topFluxVector = MHDStateVector::computeYFluxVector(topExtrapolatedValue, materialParameters);
    vector<double> bottomFluxVector = MHDStateVector::computeYFluxVector(bottomExtrapolatedValue, materialParameters);
    vector<double> evolutionVector = Solvers::computeEvolutionVector(topFluxVector, bottomFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStep(topExtrapolatedValue, bottomExtrapolatedValue, evolutionVector, side, materialParameters);
}

MHDReducedStateVector MHDSolvers::evolveStateByHalfYTimeStep(MHDReducedStateVector topStateVector, MHDReducedStateVector middleStateVector, MHDReducedStateVector bottomStateVector,
                                                             double cellSpacing, double timeStep, double bias, int slopeLimiter, int side, MHDMaterialParameters material1Parameters,
                                                             MHDMaterialParameters material2Parameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(topStateVector, middleStateVector, bottomStateVector, bias, slopeLimiter, material1Parameters, material2Parameters);
    vector<double> topExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                         VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> bottomExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                       VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> topFluxVector = MHDReducedStateVector::computeYFluxVector(topExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> bottomFluxVector = MHDReducedStateVector::computeYFluxVector(bottomExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> evolutionVector = Solvers::computeEvolutionVector(topFluxVector, bottomFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStepReduced(topExtrapolatedValue, bottomExtrapolatedValue, evolutionVector, side, material1Parameters, material2Parameters);
}

MHDStateVector MHDSolvers::evolveStateByFractionalXTimeStep(double stepFraction, MHDStateVector leftStateVector, MHDStateVector middleStateVector, MHDStateVector rightStateVector,
                                                            double cellSpacing, double timeStep, double bias, int slopeLimiter, MHDMaterialParameters materialParameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, slopeLimiter, materialParameters);

    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(materialParameters);
    vector<double> leftConservedVariableVector = VectorAlgebra::subtractVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightConservedVariableVector = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> leftFluxVector = MHDStateVector::computeXFluxVector(leftConservedVariableVector, materialParameters);
    vector<double> rightFluxVector = MHDStateVector::computeXFluxVector(rightConservedVariableVector, materialParameters);
    vector<double> evolutionVector = Solvers::computeFractionalEvolutionVector(stepFraction, leftFluxVector, rightFluxVector, cellSpacing, timeStep);

    return evolveStateByFractionalTimeStep(middleConservedVariableVector, evolutionVector, materialParameters);
}

MHDReducedStateVector MHDSolvers::evolveStateByFractionalXTimeStep(double stepFraction, MHDReducedStateVector leftStateVector, MHDReducedStateVector middleStateVector,
                                                                   MHDReducedStateVector rightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                   MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, slopeLimiter, material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> leftConservedVariableVector = VectorAlgebra::subtractVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightConservedVariableVector = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> leftFluxVector = MHDReducedStateVector::computeXFluxVector(leftConservedVariableVector, material1Parameters, material2Parameters);
    vector<double> rightFluxVector = MHDReducedStateVector::computeXFluxVector(rightConservedVariableVector, material1Parameters, material2Parameters);
    vector<double> evolutionVector = Solvers::computeFractionalEvolutionVector(stepFraction, leftFluxVector, rightFluxVector, cellSpacing, timeStep);

    return evolveStateByFractionalTimeStepReduced(middleConservedVariableVector, evolutionVector, material1Parameters, material2Parameters);
}

MHDStateVector MHDSolvers::evolveStateByFractionalYTimeStep(double stepFraction, MHDStateVector topStateVector, MHDStateVector middleStateVector, MHDStateVector bottomStateVector,
                                                            double cellSpacing, double timeStep, double bias, int slopeLimiter, MHDMaterialParameters materialParameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(topStateVector, middleStateVector, bottomStateVector, bias, slopeLimiter, materialParameters);

    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(materialParameters);
    vector<double> topConservedVariableVector = VectorAlgebra::subtractVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> bottomConservedVariableVector = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> topFluxVector = MHDStateVector::computeYFluxVector(topConservedVariableVector, materialParameters);
    vector<double> bottomFluxVector = MHDStateVector::computeYFluxVector(bottomConservedVariableVector, materialParameters);
    vector<double> evolutionVector = Solvers::computeFractionalEvolutionVector(stepFraction, topFluxVector, bottomFluxVector, cellSpacing, timeStep);

    return evolveStateByFractionalTimeStep(middleConservedVariableVector, evolutionVector, materialParameters);
}

MHDReducedStateVector MHDSolvers::evolveStateByFractionalYTimeStep(double stepFraction, MHDReducedStateVector topStateVector, MHDReducedStateVector middleStateVector,
                                                                   MHDReducedStateVector bottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                   MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(topStateVector, middleStateVector, bottomStateVector, bias, slopeLimiter, material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> topConservedVariableVector = VectorAlgebra::subtractVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> bottomConservedVariableVector = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> topFluxVector = MHDReducedStateVector::computeYFluxVector(topConservedVariableVector, material1Parameters, material2Parameters);
    vector<double> bottomFluxVector = MHDReducedStateVector::computeYFluxVector(bottomConservedVariableVector, material1Parameters, material2Parameters);
    vector<double> evolutionVector = Solvers::computeFractionalEvolutionVector(stepFraction, topFluxVector, bottomFluxVector, cellSpacing, timeStep);

    return evolveStateByFractionalTimeStepReduced(middleConservedVariableVector, evolutionVector, material1Parameters, material2Parameters);
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

MHDReducedStateVector MHDSolvers::evolveStateByHalfTimeStepReduced(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                                   MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    MHDReducedStateVector evolvedStateVector;

    if (side == 0)
    {
        evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(leftExtrapolatedValue, evolutionVector), material1Parameters, material2Parameters);
    }
    else
    {
        evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(rightExtrapolatedValue, evolutionVector), material1Parameters, material2Parameters);
    }

    return evolvedStateVector;
}

MHDStateVector MHDSolvers::evolveStateByFractionalTimeStep(vector<double> middleConservedVariableVector, vector<double> conservedVariableVectorEvolution, MHDMaterialParameters materialParameters)
{
    MHDStateVector evolvedStateVector;
    evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(middleConservedVariableVector, conservedVariableVectorEvolution), materialParameters);

    return evolvedStateVector;
}

MHDReducedStateVector MHDSolvers::evolveStateByFractionalTimeStepReduced(vector<double> middleConservedVariableVector, vector<double> conservedVariableVectorEvolution,
                                                                         MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    MHDReducedStateVector evolvedStateVector;
    evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(middleConservedVariableVector, conservedVariableVectorEvolution), material1Parameters, material2Parameters);

    return evolvedStateVector;
}

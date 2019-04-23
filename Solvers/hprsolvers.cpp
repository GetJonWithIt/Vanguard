#include "hprsolvers.h"

HPRSolvers::HPRSolvers()
{
}

vector<HPRStateVector> HPRSolvers::insertBoundaryCells(vector<HPRStateVector> & currentCells, int boundarySize)
{
    int cellCount = currentCells.size();
    vector<HPRStateVector> currentCellsWithBoundary(cellCount + (2 * boundarySize));

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

    for (int i = 0; i < cellCount; i++)
    {
        currentCellsWithBoundary[i + boundarySize] = currentCells[i];
    }

    return currentCellsWithBoundary;
}

vector<vector<HPRStateVector> > HPRSolvers::insertBoundaryCells2D(vector<vector<HPRStateVector> > & currentCells, int boundarySize)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();
    vector<vector<HPRStateVector> > currentCellsWithBoundary(rowCount + (2 * boundarySize), vector<HPRStateVector>(columnCount + (2 * boundarySize)));

    if (boundarySize == 1)
    {
        for (int i = 0; i < rowCount; i++)
        {
            currentCellsWithBoundary[i + 1][0] = currentCells[i][0];
            currentCellsWithBoundary[i + 1][columnCount + 1] = currentCells[i][columnCount - 1];
        }

        for (int i = 0; i < columnCount; i++)
        {
            currentCellsWithBoundary[0][i + 1] = currentCells[0][i];
            currentCellsWithBoundary[rowCount + 1][i + 1] = currentCells[rowCount - 1][i];
        }
    }
    else if (boundarySize == 2)
    {
        for (int i = 0; i < rowCount; i++)
        {
            currentCellsWithBoundary[i + 2][0] = currentCells[i][1];
            currentCellsWithBoundary[i + 2][1] = currentCells[i][0];

            currentCellsWithBoundary[i + 2][columnCount + 2] = currentCells[i][columnCount - 1];
            currentCellsWithBoundary[i + 2][columnCount + 3] = currentCells[i][columnCount - 2];
        }

        for (int i = 0; i < columnCount; i++)
        {
            currentCellsWithBoundary[0][i + 2] = currentCells[1][i];
            currentCellsWithBoundary[1][i + 2] = currentCells[0][i];

            currentCellsWithBoundary[rowCount + 2][i + 2] = currentCells[rowCount - 1][i];
            currentCellsWithBoundary[rowCount + 3][i + 2] = currentCells[rowCount - 2][i];
        }
    }

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            currentCellsWithBoundary[i + boundarySize][j + boundarySize] = currentCells[i][j];
        }
    }

    return currentCellsWithBoundary;
}

double HPRSolvers::computeMaximumWaveSpeed(vector<HPRStateVector> & currentCells, HPRMaterialParameters materialParameters)
{
    double maximumWaveSpeed = 0.0;
    int cellCount = currentCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        double waveSpeed = abs(currentCells[i].getXVelocity()) + HPRAcousticTensor::computeMaximumWaveSpeed(currentCells[i], materialParameters, 0);

        if (waveSpeed > maximumWaveSpeed)
        {
            maximumWaveSpeed = waveSpeed;
        }
    }

    return maximumWaveSpeed;
}

double HPRSolvers::computeMaximumWaveSpeed2D(vector<vector<HPRStateVector> > & currentCells, HPRMaterialParameters materialParameters)
{
    double maximumWaveSpeed = 0.0;
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            double waveSpeed = max(abs(currentCells[i][j].getXVelocity()), abs(currentCells[i][j].getYVelocity())) + max(
                        HPRAcousticTensor::computeMaximumWaveSpeed(currentCells[i][j], materialParameters, 0),
                        HPRAcousticTensor::computeMaximumWaveSpeed(currentCells[i][j], materialParameters, 1));

            if (waveSpeed > maximumWaveSpeed)
            {
                maximumWaveSpeed = waveSpeed;
            }
        }
    }

    return maximumWaveSpeed;
}

double HPRSolvers::computeStableTimeStep(vector<HPRStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                         HPRMaterialParameters materialParameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed(currentCells, materialParameters));

    return Solvers::computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

double HPRSolvers::computeStableTimeStep2D(vector<vector<HPRStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                           HPRMaterialParameters materialParameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed2D(currentCells, materialParameters));

    return Solvers::computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

HPRStateVector HPRSolvers::evolveStateByHalfXTimeStep(HPRStateVector leftStateVector, HPRStateVector middleStateVector, HPRStateVector rightStateVector, double cellSpacing, double timeStep,
                                                      double bias, int slopeLimiter, int side, HPRMaterialParameters materialParameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, slopeLimiter, materialParameters);
    vector<double> leftExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(materialParameters), VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(materialParameters), VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> leftFluxVector = HPRStateVector::computeXFluxVector(leftExtrapolatedValue, materialParameters);
    vector<double> rightFluxVector = HPRStateVector::computeXFluxVector(rightExtrapolatedValue, materialParameters);
    vector<double> evolutionVector = Solvers::computeEvolutionVector(leftFluxVector, rightFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStep(leftExtrapolatedValue, rightExtrapolatedValue, evolutionVector, side, materialParameters);
}

HPRStateVector HPRSolvers::evolveStateByHalfYTimeStep(HPRStateVector topStateVector, HPRStateVector middleStateVector, HPRStateVector bottomStateVector, double cellSpacing, double timeStep,
                                                      double bias, int slopeLimiter, int side, HPRMaterialParameters materialParameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(topStateVector, middleStateVector, bottomStateVector, bias, slopeLimiter, materialParameters);
    vector<double> topExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(materialParameters), VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> bottomExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(materialParameters), VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> topFluxVector = HPRStateVector::computeYFluxVector(topExtrapolatedValue, materialParameters);
    vector<double> bottomFluxVector = HPRStateVector::computeYFluxVector(bottomExtrapolatedValue, materialParameters);
    vector<double> evolutionVector = Solvers::computeEvolutionVector(topFluxVector, bottomFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStep(topExtrapolatedValue, bottomExtrapolatedValue, evolutionVector, side, materialParameters);
}

HPRStateVector HPRSolvers::evolveStateByFractionalXTimeStep(double stepFraction, HPRStateVector leftStateVector, HPRStateVector middleStateVector, HPRStateVector rightStateVector,
                                                            double cellSpacing, double timeStep, double bias, int slopeLimiter, HPRMaterialParameters materialParameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, slopeLimiter, materialParameters);

    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(materialParameters);
    vector<double> leftConservedVariableVector = VectorAlgebra::subtractVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightConservedVariableVector = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> leftFluxVector = HPRStateVector::computeXFluxVector(leftConservedVariableVector, materialParameters);
    vector<double> rightFluxVector = HPRStateVector::computeXFluxVector(rightConservedVariableVector, materialParameters);
    vector<double> evolutionVector = Solvers::computeFractionalEvolutionVector(stepFraction, leftFluxVector, rightFluxVector, cellSpacing, timeStep);

    return evolveStateByFractionalTimeStep(middleConservedVariableVector, evolutionVector, materialParameters);
}

HPRStateVector HPRSolvers::evolveStateByFractionalYTimeStep(double stepFraction, HPRStateVector topStateVector, HPRStateVector middleStateVector, HPRStateVector bottomStateVector,
                                                            double cellSpacing, double timeStep, double bias, int slopeLimiter, HPRMaterialParameters materialParameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(topStateVector, middleStateVector, bottomStateVector, bias, slopeLimiter, materialParameters);

    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(materialParameters);
    vector<double> topConservedVariableVector = VectorAlgebra::subtractVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> bottomConservedVariableVector = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> topFluxVector = HPRStateVector::computeYFluxVector(topConservedVariableVector, materialParameters);
    vector<double> bottomFluxVector = HPRStateVector::computeYFluxVector(bottomConservedVariableVector, materialParameters);
    vector<double> evolutionVector = Solvers::computeFractionalEvolutionVector(stepFraction, topFluxVector, bottomFluxVector, cellSpacing, timeStep);

    return evolveStateByFractionalTimeStep(middleConservedVariableVector, evolutionVector, materialParameters);
}

HPRStateVector HPRSolvers::evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                     HPRMaterialParameters materialParameters)
{
    HPRStateVector evolvedStateVector;

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

HPRStateVector HPRSolvers::evolveStateByFractionalTimeStep(vector<double> middleConservedVariableVector, vector<double> conservedVariableVectorEvolution, HPRMaterialParameters materialParameters)
{
    HPRStateVector evolvedStateVector;
    evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(middleConservedVariableVector, conservedVariableVectorEvolution), materialParameters);

    return evolvedStateVector;
}

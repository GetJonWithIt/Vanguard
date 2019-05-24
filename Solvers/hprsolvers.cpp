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

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        currentCellsWithBoundary[i + boundarySize] = currentCells[i];
    }

    return currentCellsWithBoundary;
}

vector<HPRIntermediateStateVector> HPRSolvers::insertBoundaryCells(vector<HPRIntermediateStateVector> & currentCells, int boundarySize)
{
    int cellCount = currentCells.size();
    vector<HPRIntermediateStateVector> currentCellsWithBoundary(cellCount + (2 * boundarySize));

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

vector<HPRReducedStateVector> HPRSolvers::insertBoundaryCells(vector<HPRReducedStateVector> & currentCells, int boundarySize)
{
    int cellCount = currentCells.size();
    vector<HPRReducedStateVector> currentCellsWithBoundary(cellCount + (2 * boundarySize));

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

vector<vector<HPRStateVector> > HPRSolvers::insertBoundaryCells2D(vector<vector<HPRStateVector> > & currentCells, int boundarySize)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();
    vector<vector<HPRStateVector> > currentCellsWithBoundary(rowCount + (2 * boundarySize), vector<HPRStateVector>(columnCount + (2 * boundarySize)));

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

vector<vector<HPRIntermediateStateVector> > HPRSolvers::insertBoundaryCells2D(vector<vector<HPRIntermediateStateVector> > & currentCells, int boundarySize)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();
    vector<vector<HPRIntermediateStateVector> > currentCellsWithBoundary(rowCount + (2 * boundarySize), vector<HPRIntermediateStateVector>(columnCount + (2 * boundarySize)));

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
        for (int i = 0 ; i < columnCount; i++)
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

vector<vector<HPRReducedStateVector> > HPRSolvers::insertBoundaryCells2D(vector<vector<HPRReducedStateVector> > & currentCells, int boundarySize)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();
    vector<vector<HPRReducedStateVector> > currentCellsWithBoundary(rowCount + (2 * boundarySize), vector<HPRReducedStateVector>(columnCount + (2 * boundarySize)));

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

double HPRSolvers::computeMaximumWaveSpeed(vector<HPRStateVector> & currentCells, HPRMaterialParameters materialParameters)
{
    double maximumWaveSpeed = 0.0;
    int cellCount = currentCells.size();

#pragma omp parallel for
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

double HPRSolvers::computeMaximumWaveSpeed(vector<HPRIntermediateStateVector> & currentCells, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    double maximumWaveSpeed = 0.0;
    int cellCount = currentCells.size();

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        double waveSpeed = abs(currentCells[i].getInterfaceXVelocity()) + HPRIntermediateAcousticTensor::computeMaximumWaveSpeed(currentCells[i], material1Parameters, material2Parameters, 0);

        if (waveSpeed > maximumWaveSpeed)
        {
            maximumWaveSpeed = waveSpeed;
        }
    }

    return maximumWaveSpeed;
}

double HPRSolvers::computeMaximumWaveSpeed(vector<HPRReducedStateVector> & currentCells, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    double maximumWaveSpeed = 0.0;
    int cellCount = currentCells.size();

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        double waveSpeed = abs(currentCells[i].getInterfaceXVelocity()) + HPRReducedAcousticTensor::computeMaximumWaveSpeed(currentCells[i], material1Parameters, material2Parameters, 0);

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

#pragma omp parallel for
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

double HPRSolvers::computeMaximumWaveSpeed2D(vector<vector<HPRIntermediateStateVector> > & currentCells, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    double maximumWaveSpeed = 0.0;
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            double waveSpeed = max(abs(currentCells[i][j].getInterfaceXVelocity()), abs(currentCells[i][j].getInterfaceYVelocity())) + max(
                        HPRIntermediateAcousticTensor::computeMaximumWaveSpeed(currentCells[i][j], material1Parameters, material2Parameters, 0),
                        HPRIntermediateAcousticTensor::computeMaximumWaveSpeed(currentCells[i][j], material1Parameters, material2Parameters, 1));

            if (waveSpeed > maximumWaveSpeed)
            {
                maximumWaveSpeed = waveSpeed;
            }
        }
    }

    return maximumWaveSpeed;
}

double HPRSolvers::computeMaximumWaveSpeed2D(vector<vector<HPRReducedStateVector> > & currentCells, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    double maximumWaveSpeed = 0.0;
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            double waveSpeed = max(abs(currentCells[i][j].getInterfaceXVelocity()), abs(currentCells[i][j].getInterfaceYVelocity())) + max(
                        HPRReducedAcousticTensor::computeMaximumWaveSpeed(currentCells[i][j], material1Parameters, material2Parameters, 0),
                        HPRReducedAcousticTensor::computeMaximumWaveSpeed(currentCells[i][j], material1Parameters, material2Parameters, 1));

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

double HPRSolvers::computeStableTimeStep(vector<HPRIntermediateStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                         HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed(currentCells, material1Parameters, material2Parameters));

    return Solvers::computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

double HPRSolvers::computeStableTimeStep(vector<HPRReducedStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                         HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed(currentCells, material1Parameters, material2Parameters));

    return Solvers::computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

double HPRSolvers::computeStableTimeStep2D(vector<vector<HPRStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                           HPRMaterialParameters materialParameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed2D(currentCells, materialParameters));

    return Solvers::computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

double HPRSolvers::computeStableTimeStep2D(vector<vector<HPRIntermediateStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                           int currentIteration, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed2D(currentCells, material1Parameters, material2Parameters));

    return Solvers::computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

double HPRSolvers::computeStableTimeStep2D(vector<vector<HPRReducedStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                           int currentIteration, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed2D(currentCells, material1Parameters, material2Parameters));

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

HPRIntermediateStateVector HPRSolvers::evolveStateByHalfXTimeStep(HPRIntermediateStateVector leftStateVector, HPRIntermediateStateVector middleStateVector,
                                                                  HPRIntermediateStateVector rightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter, int side,
                                                                  HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, slopeLimiter, material1Parameters, material2Parameters);
    vector<double> leftExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                          VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                      VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> leftFluxVector = HPRIntermediateStateVector::computeXFluxVector(leftExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> rightFluxVector = HPRIntermediateStateVector::computeXFluxVector(rightExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> evolutionVector = Solvers::computeEvolutionVector(leftFluxVector, rightFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStepIntermediate(leftExtrapolatedValue, rightExtrapolatedValue, evolutionVector, side, material1Parameters, material2Parameters);
}

HPRReducedStateVector HPRSolvers::evolveStateByHalfXTimeStep(HPRReducedStateVector leftStateVector, HPRReducedStateVector middleStateVector, HPRReducedStateVector rightStateVector,
                                                             double cellSpacing, double timeStep, double bias, int slopeLimiter, int side, HPRMaterialParameters material1Parameters,
                                                             HPRMaterialParameters material2Parameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, slopeLimiter, material1Parameters, material2Parameters);
    vector<double> leftExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                          VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                      VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> leftFluxVector = HPRReducedStateVector::computeXFluxVector(leftExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> rightFluxVector = HPRReducedStateVector::computeXFluxVector(rightExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> evolutionVector = Solvers::computeEvolutionVector(leftFluxVector, rightFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStepReduced(leftExtrapolatedValue, rightExtrapolatedValue, evolutionVector, side, material1Parameters, material2Parameters);
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

HPRIntermediateStateVector HPRSolvers::evolveStateByHalfYTimeStep(HPRIntermediateStateVector topStateVector, HPRIntermediateStateVector middleStateVector,
                                                                  HPRIntermediateStateVector bottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter, int side,
                                                                  HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(topStateVector, middleStateVector, bottomStateVector, bias, slopeLimiter, material1Parameters, material2Parameters);
    vector<double> topExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                         VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> bottomExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                       VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> topFluxVector = HPRIntermediateStateVector::computeYFluxVector(topExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> bottomFluxVector = HPRIntermediateStateVector::computeYFluxVector(bottomExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> evolutionVector = Solvers::computeEvolutionVector(topFluxVector, bottomFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStepIntermediate(topExtrapolatedValue, bottomExtrapolatedValue, evolutionVector, side, material1Parameters, material2Parameters);
}

HPRReducedStateVector HPRSolvers::evolveStateByHalfYTimeStep(HPRReducedStateVector topStateVector, HPRReducedStateVector middleStateVector, HPRReducedStateVector bottomStateVector,
                                                             double cellSpacing, double timeStep, double bias, int slopeLimiter, int side, HPRMaterialParameters material1Parameters,
                                                             HPRMaterialParameters material2Parameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(topStateVector, middleStateVector, bottomStateVector, bias, slopeLimiter, material1Parameters, material2Parameters);
    vector<double> topExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                         VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> bottomExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                       VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> topFluxVector = HPRReducedStateVector::computeYFluxVector(topExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> bottomFluxVector = HPRReducedStateVector::computeYFluxVector(bottomExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> evolutionVector = Solvers::computeEvolutionVector(topFluxVector, bottomFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStepReduced(topExtrapolatedValue, bottomExtrapolatedValue, evolutionVector, side, material1Parameters, material2Parameters);
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

HPRIntermediateStateVector HPRSolvers::evolveStateByFractionalXTimeStep(double stepFraction, HPRIntermediateStateVector leftStateVector, HPRIntermediateStateVector middleStateVector,
                                                                        HPRIntermediateStateVector rightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                        HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, slopeLimiter, material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> leftConservedVariableVector = VectorAlgebra::subtractVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightConservedVariableVector = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> leftFluxVector = HPRIntermediateStateVector::computeXFluxVector(leftConservedVariableVector, material1Parameters, material2Parameters);
    vector<double> rightFluxVector = HPRIntermediateStateVector::computeXFluxVector(rightConservedVariableVector, material1Parameters, material2Parameters);
    vector<double> evolutionVector = Solvers::computeFractionalEvolutionVector(stepFraction, leftFluxVector, rightFluxVector, cellSpacing, timeStep);

    return evolveStateByFractionalTimeStepIntermediate(middleConservedVariableVector, evolutionVector, material1Parameters, material2Parameters);
}

HPRReducedStateVector HPRSolvers::evolveStateByFractionalXTimeStep(double stepFraction, HPRReducedStateVector leftStateVector, HPRReducedStateVector middleStateVector,
                                                                   HPRReducedStateVector rightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                   HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, slopeLimiter, material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> leftConservedVariableVector = VectorAlgebra::subtractVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightConservedVariableVector = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> leftFluxVector = HPRReducedStateVector::computeXFluxVector(leftConservedVariableVector, material1Parameters, material2Parameters);
    vector<double> rightFluxVector = HPRReducedStateVector::computeXFluxVector(rightConservedVariableVector, material1Parameters, material2Parameters);
    vector<double> evolutionVector = Solvers::computeFractionalEvolutionVector(stepFraction, leftFluxVector, rightFluxVector, cellSpacing, timeStep);

    return evolveStateByFractionalTimeStepReduced(middleConservedVariableVector, evolutionVector, material1Parameters, material2Parameters);
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

HPRIntermediateStateVector HPRSolvers::evolveStateByFractionalYTimeStep(double stepFraction, HPRIntermediateStateVector topStateVector, HPRIntermediateStateVector middleStateVector,
                                                                        HPRIntermediateStateVector bottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                        HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(topStateVector, middleStateVector, bottomStateVector, bias, slopeLimiter, material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> topConservedVariableVector = VectorAlgebra::subtractVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> bottomConservedVariableVector = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> topFluxVector = HPRIntermediateStateVector::computeYFluxVector(topConservedVariableVector, material1Parameters, material2Parameters);
    vector<double> bottomFluxVector = HPRIntermediateStateVector::computeYFluxVector(bottomConservedVariableVector, material1Parameters, material2Parameters);
    vector<double> evolutionVector = Solvers::computeFractionalEvolutionVector(stepFraction, topFluxVector, bottomFluxVector, cellSpacing, timeStep);

    return evolveStateByFractionalTimeStepIntermediate(middleConservedVariableVector, evolutionVector, material1Parameters, material2Parameters);
}

HPRReducedStateVector HPRSolvers::evolveStateByFractionalYTimeStep(double stepFraction, HPRReducedStateVector topStateVector, HPRReducedStateVector middleStateVector,
                                                                   HPRReducedStateVector bottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                   HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(topStateVector, middleStateVector, bottomStateVector, bias, slopeLimiter, material1Parameters, material2Parameters);

    vector<double> middleConservedVariableVector = middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters);
    vector<double> topConservedVariableVector = VectorAlgebra::subtractVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> bottomConservedVariableVector = VectorAlgebra::addVectors(middleConservedVariableVector, VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> topFluxVector = HPRReducedStateVector::computeYFluxVector(topConservedVariableVector, material1Parameters, material2Parameters);
    vector<double> bottomFluxVector = HPRReducedStateVector::computeYFluxVector(bottomConservedVariableVector, material1Parameters, material2Parameters);
    vector<double> evolutionVector = Solvers::computeFractionalEvolutionVector(stepFraction, topFluxVector, bottomFluxVector, cellSpacing, timeStep);

    return evolveStateByFractionalTimeStepReduced(middleConservedVariableVector, evolutionVector, material1Parameters, material2Parameters);
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

HPRIntermediateStateVector HPRSolvers::evolveStateByHalfTimeStepIntermediate(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector,
                                                                             int side, HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    HPRIntermediateStateVector evolvedStateVector;

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

HPRReducedStateVector HPRSolvers::evolveStateByHalfTimeStepReduced(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                                   HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    HPRReducedStateVector evolvedStateVector;

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

HPRStateVector HPRSolvers::evolveStateByFractionalTimeStep(vector<double> middleConservedVariableVector, vector<double> conservedVariableVectorEvolution, HPRMaterialParameters materialParameters)
{
    HPRStateVector evolvedStateVector;
    evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(middleConservedVariableVector, conservedVariableVectorEvolution), materialParameters);

    return evolvedStateVector;
}

HPRIntermediateStateVector HPRSolvers::evolveStateByFractionalTimeStepIntermediate(vector<double> middleConservedVariableVector, vector<double> conservedVariableVectorEvolution,
                                                                                   HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    HPRIntermediateStateVector evolvedStateVector;
    evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(middleConservedVariableVector, conservedVariableVectorEvolution), material1Parameters, material2Parameters);

    return evolvedStateVector;
}

HPRReducedStateVector HPRSolvers::evolveStateByFractionalTimeStepReduced(vector<double> middleConservedVariableVector, vector<double> conservedVariableVectorEvolution,
                                                                         HPRMaterialParameters material1Parameters, HPRMaterialParameters material2Parameters)
{
    HPRReducedStateVector evolvedStateVector;
    evolvedStateVector.setConservedVariableVector(VectorAlgebra::addVectors(middleConservedVariableVector, conservedVariableVectorEvolution), material1Parameters, material2Parameters);

    return evolvedStateVector;
}

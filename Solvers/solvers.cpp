#include "solvers.h"

Solvers::Solvers()
{
}

vector<EulerStateVector> Solvers::insertBoundaryCells(vector<EulerStateVector> & currentCells, int boundarySize)
{
    int cellCount = currentCells.size();
    vector<EulerStateVector> currentCellsWithBoundary(cellCount + (2 * boundarySize));

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

vector<EulerMultiphysicsStateVector> Solvers::insertBoundaryCells(vector<EulerMultiphysicsStateVector> & currentCells, int boundarySize)
{
    int cellCount = currentCells.size();
    vector<EulerMultiphysicsStateVector> currentCellsWithBoundary(cellCount + (2 * boundarySize));

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

vector<EulerReducedStateVector> Solvers::insertBoundaryCells(vector<EulerReducedStateVector> & currentCells, int boundarySize)
{
    int cellCount = currentCells.size();
    vector<EulerReducedStateVector> currentCellsWithBoundary(cellCount + (2 * boundarySize));

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

vector<vector<EulerStateVector> > Solvers::insertBoundaryCells2D(vector<vector<EulerStateVector> > & currentCells, int boundarySize)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();
    vector<vector<EulerStateVector> > currentCellsWithBoundary(rowCount + (2 * boundarySize), vector<EulerStateVector>(columnCount + (2 * boundarySize)));

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

vector<vector<EulerMultiphysicsStateVector> > Solvers::insertBoundaryCells2D(vector<vector<EulerMultiphysicsStateVector> > & currentCells, int boundarySize)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();
    vector<vector<EulerMultiphysicsStateVector> > currentCellsWithBoundary(rowCount + (2 * boundarySize), vector<EulerMultiphysicsStateVector>(columnCount + (2 * boundarySize)));

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

vector<vector<EulerReducedStateVector> > Solvers::insertBoundaryCells2D(vector<vector<EulerReducedStateVector> > & currentCells, int boundarySize)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();
    vector<vector<EulerReducedStateVector> > currentCellsWithBoundary(rowCount + (2 * boundarySize), vector<EulerReducedStateVector>(columnCount + (2 * boundarySize)));

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

double Solvers::computeMaximumWaveSpeed(vector<EulerStateVector> & currentCells, EulerMaterialParameters materialParameters)
{
    double maximumWaveSpeed = 0.0;
    int cellCount = currentCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        double waveSpeed = abs(currentCells[i].getXVelocity()) + currentCells[i].computeSoundSpeed(materialParameters);

        if (waveSpeed > maximumWaveSpeed)
        {
            maximumWaveSpeed = waveSpeed;
        }
    }

    return maximumWaveSpeed;
}

double Solvers::computeMaximumWaveSpeed(vector<EulerMultiphysicsStateVector> & currentCells, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    double maximumWaveSpeed = 0.0;
    int cellCount = currentCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        double waveSpeed = abs(currentCells[i].getInterfaceXVelocity()) + max(currentCells[i].computeMaterial1SoundSpeed(material1Parameters),
                                                                              currentCells[i].computeMaterial2SoundSpeed(material2Parameters));

        if (waveSpeed > maximumWaveSpeed)
        {
            maximumWaveSpeed = waveSpeed;
        }
    }

    return maximumWaveSpeed;
}

double Solvers::computeMaximumWaveSpeed(vector<EulerReducedStateVector> & currentCells, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    double maximumWaveSpeed = 0.0;
    int cellCount = currentCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        double waveSpeed = abs(currentCells[i].getInterfaceXVelocity()) + max(currentCells[i].computeMaterial1SoundSpeed(material1Parameters),
                                                                              currentCells[i].computeMaterial2SoundSpeed(material2Parameters));

        if (waveSpeed > maximumWaveSpeed)
        {
            maximumWaveSpeed = waveSpeed;
        }
    }

    return maximumWaveSpeed;
}

double Solvers::computeMaximumWaveSpeed2D(vector<vector<EulerStateVector> > & currentCells, EulerMaterialParameters materialParameters)
{
    double maximumWaveSpeed = 0.0;
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            double waveSpeed = max(abs(currentCells[i][j].getXVelocity()), abs(currentCells[i][j].getYVelocity())) + currentCells[i][j].computeSoundSpeed(materialParameters);

            if (waveSpeed > maximumWaveSpeed)
            {
                maximumWaveSpeed = waveSpeed;
            }
        }
    }

    return maximumWaveSpeed;
}

double Solvers::computeMaximumWaveSpeed2D(vector<vector<EulerMultiphysicsStateVector> > & currentCells, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    double maximumWaveSpeed = 0.0;
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            double waveSpeed = max(abs(currentCells[i][j].getInterfaceXVelocity()), abs(currentCells[i][j].getInterfaceYVelocity())) +
                    max(currentCells[i][j].computeMaterial1SoundSpeed(material1Parameters), currentCells[i][j].computeMaterial2SoundSpeed(material2Parameters));

            if (waveSpeed > maximumWaveSpeed)
            {
                maximumWaveSpeed = waveSpeed;
            }
        }
    }

    return maximumWaveSpeed;
}

double Solvers::computeMaximumWaveSpeed2D(vector<vector<EulerReducedStateVector> > & currentCells, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    double maximumWaveSpeed = 0.0;
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            double waveSpeed = max(abs(currentCells[i][j].getInterfaceXVelocity()), abs(currentCells[i][j].getInterfaceYVelocity())) +
                    max(currentCells[i][j].computeMaterial1SoundSpeed(material1Parameters), currentCells[i][j].computeMaterial2SoundSpeed(material2Parameters));

            if (waveSpeed > maximumWaveSpeed)
            {
                maximumWaveSpeed = waveSpeed;
            }
        }
    }

    return maximumWaveSpeed;
}

double Solvers::computeStableTimeStep(vector<EulerStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                      EulerMaterialParameters materialParameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed(currentCells, materialParameters));

    return computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

double Solvers::computeStableTimeStep(vector<EulerMultiphysicsStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                      EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed(currentCells, material1Parameters, material2Parameters));

    return computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

double Solvers::computeStableTimeStep(vector<EulerReducedStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                      EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed(currentCells, material1Parameters, material2Parameters));

    return computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

double Solvers::computeStableTimeStep2D(vector<vector<EulerStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration,
                                        EulerMaterialParameters materialParameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed2D(currentCells, materialParameters));

    return computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

double Solvers::computeStableTimeStep2D(vector<vector<EulerMultiphysicsStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                        int currentIteration, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed2D(currentCells, material1Parameters, material2Parameters));

    return computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

double Solvers::computeStableTimeStep2D(vector<vector<EulerReducedStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                        int currentIteration, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed2D(currentCells, material1Parameters, material2Parameters));

    return computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

double Solvers::computeStableTimeStep(double timeStep, double currentTime, double finalTime, int currentIteration)
{
    double newTimeStep = timeStep;

    if (currentIteration <= 5)
    {
        newTimeStep *= 0.2;
    }

    if ((currentTime + newTimeStep) > finalTime)
    {
        newTimeStep = finalTime - currentTime;
    }

    return newTimeStep;
}

EulerStateVector Solvers::evolveStateByHalfXTimeStep(EulerStateVector leftStateVector, EulerStateVector middleStateVector, EulerStateVector rightStateVector, double cellSpacing,
                                                     double timeStep, double bias, int slopeLimiter, int side, EulerMaterialParameters materialParameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, slopeLimiter, materialParameters);
    vector<double> leftExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(materialParameters), VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(materialParameters), VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> leftFluxVector = EulerStateVector::computeXFluxVector(leftExtrapolatedValue, materialParameters);
    vector<double> rightFluxVector = EulerStateVector::computeXFluxVector(rightExtrapolatedValue, materialParameters);
    vector<double> evolutionVector = computeEvolutionVector(leftFluxVector, rightFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStep(leftExtrapolatedValue, rightExtrapolatedValue, evolutionVector, side, materialParameters);
}

EulerMultiphysicsStateVector Solvers::evolveStateByHalfXTimeStep(EulerMultiphysicsStateVector leftStateVector, EulerMultiphysicsStateVector middleStateVector,
                                                                 EulerMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter, int side,
                                                                 EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, slopeLimiter, material1Parameters, material2Parameters);
    vector<double> leftExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                          VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                      VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> leftFluxVector = EulerMultiphysicsStateVector::computeXFluxVector(leftExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> rightFluxVector = EulerMultiphysicsStateVector::computeXFluxVector(rightExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> evolutionVector = computeEvolutionVector(leftFluxVector, rightFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStep(leftExtrapolatedValue, rightExtrapolatedValue, evolutionVector, side, material1Parameters, material2Parameters);
}

EulerReducedStateVector Solvers::evolveStateByHalfXTimeStep(EulerReducedStateVector leftStateVector, EulerReducedStateVector middleStateVector, EulerReducedStateVector rightStateVector,
                                                            double cellSpacing, double timeStep, double bias, int slopeLimiter, int side, EulerMaterialParameters material1Parameters,
                                                            EulerMaterialParameters material2Parameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, slopeLimiter, material1Parameters, material2Parameters);
    vector<double> leftExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                          VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                      VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> leftFluxVector = EulerReducedStateVector::computeXFluxVector(leftExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> rightFluxVector = EulerReducedStateVector::computeXFluxVector(rightExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> evolutionVector = computeEvolutionVector(leftFluxVector, rightFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStepReduced(leftExtrapolatedValue, rightExtrapolatedValue, evolutionVector, side, material1Parameters, material2Parameters);
}

EulerStateVector Solvers::evolveStateByHalfYTimeStep(EulerStateVector topStateVector, EulerStateVector middleStateVector, EulerStateVector bottomStateVector, double cellSpacing,
                                                     double timeStep, double bias, int slopeLimiter, int side, EulerMaterialParameters materialParameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(topStateVector, middleStateVector, bottomStateVector, bias, slopeLimiter, materialParameters);
    vector<double> topExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(materialParameters), VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> bottomExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(materialParameters), VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> topFluxVector = EulerStateVector::computeYFluxVector(topExtrapolatedValue, materialParameters);
    vector<double> bottomFluxVector = EulerStateVector::computeYFluxVector(bottomExtrapolatedValue, materialParameters);
    vector<double> evolutionVector = computeEvolutionVector(topFluxVector, bottomFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStep(topExtrapolatedValue, bottomExtrapolatedValue, evolutionVector, side, materialParameters);
}

EulerMultiphysicsStateVector Solvers::evolveStateByHalfYTimeStep(EulerMultiphysicsStateVector topStateVector, EulerMultiphysicsStateVector middleStateVector,
                                                                 EulerMultiphysicsStateVector bottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter, int side,
                                                                 EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(topStateVector, middleStateVector, bottomStateVector, bias, slopeLimiter, material1Parameters, material2Parameters);
    vector<double> topExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                         VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> bottomExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                       VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> topFluxVector = EulerMultiphysicsStateVector::computeYFluxVector(topExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> bottomFluxVector = EulerMultiphysicsStateVector::computeYFluxVector(bottomExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> evolutionVector = computeEvolutionVector(topFluxVector, bottomFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStep(topExtrapolatedValue, bottomExtrapolatedValue, evolutionVector, side, material1Parameters, material2Parameters);
}

EulerReducedStateVector Solvers::evolveStateByHalfYTimeStep(EulerReducedStateVector topStateVector, EulerReducedStateVector middleStateVector, EulerReducedStateVector bottomStateVector,
                                                            double cellSpacing, double timeStep, double bias, int slopeLimiter, int side, EulerMaterialParameters material1Parameters,
                                                            EulerMaterialParameters material2Parameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(topStateVector, middleStateVector, bottomStateVector, bias, slopeLimiter, material1Parameters, material2Parameters);
    vector<double> topExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                         VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> bottomExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                       VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> topFluxVector = EulerReducedStateVector::computeYFluxVector(topExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> bottomFluxVector = EulerReducedStateVector::computeYFluxVector(bottomExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> evolutionVector = computeEvolutionVector(topFluxVector, bottomFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStepReduced(topExtrapolatedValue, bottomExtrapolatedValue, evolutionVector, side, material1Parameters, material2Parameters);
}

EulerStateVector Solvers::evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                    EulerMaterialParameters materialParameters)
{
    EulerStateVector evolvedStateVector;

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

EulerMultiphysicsStateVector Solvers::evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                                EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    EulerMultiphysicsStateVector evolvedStateVector;

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

EulerReducedStateVector Solvers::evolveStateByHalfTimeStepReduced(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector, int side,
                                                                  EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    EulerReducedStateVector evolvedStateVector;

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

vector<double> Solvers::computeEvolutionVector(vector<double> leftFluxVector, vector<double> rightFluxVector, double cellSpacing, double timeStep)
{
    return computeFractionalEvolutionVector(0.5, leftFluxVector, rightFluxVector, cellSpacing, timeStep);
}

vector<double> Solvers::computeFractionalEvolutionVector(double stepFraction, vector<double> leftFluxVector, vector<double> rightFluxVector, double cellSpacing, double timeStep)
{
    return VectorAlgebra::multiplyVector(stepFraction * (timeStep / cellSpacing), VectorAlgebra::subtractVectors(leftFluxVector, rightFluxVector));
}

void Solvers::outputStatus(int currentIteration, double currentTime, double timeStep)
{
    cout << "Iteration = " << currentIteration << "; Time = " << currentTime << "; Timestep = " << timeStep << endl;
}

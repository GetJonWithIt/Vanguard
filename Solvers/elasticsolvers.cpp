#include "elasticsolvers.h"

ElasticSolvers::ElasticSolvers()
{
}

vector<ElasticStateVector> ElasticSolvers::insertBoundaryCells(vector<ElasticStateVector> & currentCells, int boundarySize)
{
    int cellCount = currentCells.size();
    vector<ElasticStateVector> currentCellsWithBoundary(cellCount + (2 * boundarySize));

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

vector<ElasticMultiphysicsStateVector> ElasticSolvers::insertBoundaryCells(vector<ElasticMultiphysicsStateVector> & currentCells, int boundarySize)
{
    int cellCount = currentCells.size();
    vector<ElasticMultiphysicsStateVector> currentCellsWithBoundary(cellCount + (2 * boundarySize));

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

vector<ElasticReducedStateVector> ElasticSolvers::insertBoundaryCells(vector<ElasticReducedStateVector> & currentCells, int boundarySize)
{
    int cellCount = currentCells.size();
    vector<ElasticReducedStateVector> currentCellsWithBoundary(cellCount + (2 * boundarySize));

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

vector<vector<ElasticStateVector> > ElasticSolvers::insertBoundaryCells2D(vector<vector<ElasticStateVector> > & currentCells, int boundarySize)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();
    vector<vector<ElasticStateVector> > currentCellsWithBoundary(rowCount + (2 * boundarySize), vector<ElasticStateVector>(columnCount + (2 * boundarySize)));

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

vector<vector<ElasticReducedStateVector> > ElasticSolvers::insertBoundaryCells2D(vector<vector<ElasticReducedStateVector> > & currentCells, int boundarySize)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();
    vector<vector<ElasticReducedStateVector> > currentCellsWithBoundary(rowCount + (2 * boundarySize), vector<ElasticReducedStateVector>(columnCount + (2 * boundarySize)));

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

double ElasticSolvers::computeMaximumWaveSpeed(vector<ElasticStateVector> & currentCells, HyperelasticMaterialParameters materialParameters)
{
    double maximumWaveSpeed = 0.0;
    int cellCount = currentCells.size();

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        double waveSpeed = abs(currentCells[i].getXVelocity()) + currentCells[i].computeSoundSpeed(materialParameters, 0);

        if (waveSpeed > maximumWaveSpeed)
        {
            maximumWaveSpeed = waveSpeed;
        }
    }

    return maximumWaveSpeed;
}

double ElasticSolvers::computeMaximumWaveSpeed(vector<ElasticMultiphysicsStateVector> & currentCells, HyperelasticMaterialParameters material1Parameters,
                                               HyperelasticMaterialParameters material2Parameters)
{
    double maximumWaveSpeed = 0.0;
    int cellCount = currentCells.size();

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        double waveSpeed = abs(currentCells[i].getInterfaceXVelocity()) + max(currentCells[i].computeMaterial1SoundSpeed(material1Parameters, 0),
                                                                              currentCells[i].computeMaterial2SoundSpeed(material2Parameters, 0));

        if (waveSpeed > maximumWaveSpeed)
        {
            maximumWaveSpeed = waveSpeed;
        }
    }

    return maximumWaveSpeed;
}

double ElasticSolvers::computeMaximumWaveSpeed(vector<ElasticReducedStateVector> & currentCells, HyperelasticMaterialParameters material1Parameters,
                                               HyperelasticMaterialParameters material2Parameters)
{
    double maximumWaveSpeed = 0.0;
    int cellCount = currentCells.size();

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        double waveSpeed = abs(currentCells[i].getInterfaceXVelocity()) + max(currentCells[i].computeMaterial1SoundSpeed(material1Parameters, 0),
                                                                              currentCells[i].computeMaterial2SoundSpeed(material2Parameters, 0));

        if (waveSpeed > maximumWaveSpeed)
        {
            maximumWaveSpeed = waveSpeed;
        }
    }

    return maximumWaveSpeed;
}

double ElasticSolvers::computeMaximumWaveSpeed2D(vector<vector<ElasticStateVector> > & currentCells, HyperelasticMaterialParameters materialParameters)
{
    double maximumWaveSpeed = 0.0;
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            double waveSpeed = max(abs(currentCells[i][j].getXVelocity()), abs(currentCells[i][j].getYVelocity())) + max(currentCells[i][j].computeSoundSpeed(materialParameters, 0),
                                                                                                                         currentCells[i][j].computeSoundSpeed(materialParameters, 1));

            if (waveSpeed > maximumWaveSpeed)
            {
                maximumWaveSpeed = waveSpeed;
            }
        }
    }

    return maximumWaveSpeed;
}

double ElasticSolvers::computeMaximumWaveSpeed2D(vector<vector<ElasticReducedStateVector> > & currentCells, HyperelasticMaterialParameters material1Parameters,
                                                 HyperelasticMaterialParameters material2Parameters)
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
                    max(max(currentCells[i][j].computeMaterial1SoundSpeed(material1Parameters, 0), currentCells[i][j].computeMaterial1SoundSpeed(material1Parameters, 1)),
                        max(currentCells[i][j].computeMaterial2SoundSpeed(material2Parameters, 0), currentCells[i][j].computeMaterial2SoundSpeed(material2Parameters, 1)));

            if (waveSpeed > maximumWaveSpeed)
            {
                maximumWaveSpeed = waveSpeed;
            }
        }
    }

    return maximumWaveSpeed;
}

double ElasticSolvers::computeStableTimeStep(vector<ElasticStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                             int currentIteration, HyperelasticMaterialParameters materialParameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed(currentCells, materialParameters));

    return Solvers::computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

double ElasticSolvers::computeStableTimeStep(vector<ElasticMultiphysicsStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime,
                                             double finalTime, int currentIteration, HyperelasticMaterialParameters material1Parameters,
                                             HyperelasticMaterialParameters material2Parameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed(currentCells, material1Parameters, material2Parameters));

    return Solvers::computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

double ElasticSolvers::computeStableTimeStep(vector<ElasticReducedStateVector> & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                             int currentIteration, HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed(currentCells, material1Parameters, material2Parameters));

    return Solvers::computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

double ElasticSolvers::computeStableTimeStep2D(vector<vector<ElasticStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                        int currentIteration, HyperelasticMaterialParameters materialParameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed2D(currentCells, materialParameters));

    return Solvers::computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

double ElasticSolvers::computeStableTimeStep2D(vector<vector<ElasticReducedStateVector> > & currentCells, double cellSpacing, double CFLCoefficient, double currentTime,
                                               double finalTime, int currentIteration, HyperelasticMaterialParameters material1Parameters,
                                               HyperelasticMaterialParameters material2Parameters)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed2D(currentCells, material1Parameters, material2Parameters));

    return Solvers::computeStableTimeStep(timeStep, currentTime, finalTime, currentIteration);
}

ElasticStateVector ElasticSolvers::evolveStateByHalfXTimeStep(ElasticStateVector leftStateVector, ElasticStateVector middleStateVector, ElasticStateVector rightStateVector,
                                                              double cellSpacing, double timeStep, double bias, int slopeLimiter, int side,
                                                              HyperelasticMaterialParameters materialParameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, slopeLimiter, materialParameters);
    vector<double> leftExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(materialParameters),
                                                                          VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightExtrapolatedValue  = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(materialParameters),
                                                                       VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> leftFluxVector = ElasticStateVector::computeXFluxVector(leftExtrapolatedValue, materialParameters);
    vector<double> rightFluxVector = ElasticStateVector::computeXFluxVector(rightExtrapolatedValue, materialParameters);
    vector<double> evolutionVector = Solvers::computeEvolutionVector(leftFluxVector, rightFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStep(leftExtrapolatedValue, rightExtrapolatedValue, evolutionVector, side, materialParameters);
}

ElasticMultiphysicsStateVector ElasticSolvers::evolveStateByHalfXTimeStep(ElasticMultiphysicsStateVector leftStateVector, ElasticMultiphysicsStateVector middleStateVector,
                                                                          ElasticMultiphysicsStateVector rightStateVector, double cellSpacing, double timeStep, double bias,
                                                                          int slopeLimiter, int side, HyperelasticMaterialParameters material1Parameters,
                                                                          HyperelasticMaterialParameters material2Parameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, slopeLimiter, material1Parameters,
                                                                   material2Parameters);
    vector<double> leftExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                          VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                      VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> leftFluxVector = ElasticMultiphysicsStateVector::computeXFluxVector(leftExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> rightFluxVector = ElasticMultiphysicsStateVector::computeXFluxVector(rightExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> evolutionVector = Solvers::computeEvolutionVector(leftFluxVector, rightFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStep(leftExtrapolatedValue, rightExtrapolatedValue, evolutionVector, side, material1Parameters, material2Parameters);
}

ElasticReducedStateVector ElasticSolvers::evolveStateByHalfXTimeStep(ElasticReducedStateVector leftStateVector, ElasticReducedStateVector middleStateVector,
                                                                     ElasticReducedStateVector rightStateVector, double cellSpacing, double timeStep, double bias,
                                                                     int slopeLimiter, int side, HyperelasticMaterialParameters material1Parameters,
                                                                     HyperelasticMaterialParameters material2Parameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, slopeLimiter, material1Parameters,
                                                                   material2Parameters);
    vector<double> leftExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                          VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                      VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> leftFluxVector = ElasticReducedStateVector::computeXFluxVector(leftExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> rightFluxVector = ElasticReducedStateVector::computeXFluxVector(rightExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> evolutionVector = Solvers::computeEvolutionVector(leftFluxVector, rightFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStepReduced(leftExtrapolatedValue, rightExtrapolatedValue, evolutionVector, side, material1Parameters, material2Parameters);
}

ElasticStateVector ElasticSolvers::evolveStateByHalfYTimeStep(ElasticStateVector topStateVector, ElasticStateVector middleStateVector, ElasticStateVector bottomStateVector,
                                                              double cellSpacing, double timeStep, double bias, int slopeLimiter, int side,
                                                              HyperelasticMaterialParameters materialParameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(topStateVector, middleStateVector, bottomStateVector, bias, slopeLimiter, materialParameters);
    vector<double> topExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(materialParameters),
                                                                         VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> bottomExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(materialParameters),
                                                                       VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> topFluxVector = ElasticStateVector::computeYFluxVector(topExtrapolatedValue, materialParameters);
    vector<double> bottomFluxVector = ElasticStateVector::computeYFluxVector(bottomExtrapolatedValue, materialParameters);
    vector<double> evolutionVector = Solvers::computeEvolutionVector(topFluxVector, bottomFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStep(topExtrapolatedValue, bottomExtrapolatedValue, evolutionVector, side, materialParameters);
}

ElasticReducedStateVector ElasticSolvers::evolveStateByHalfYTimeStep(ElasticReducedStateVector topStateVector, ElasticReducedStateVector middleStateVector,
                                                                     ElasticReducedStateVector bottomStateVector, double cellSpacing, double timeStep, double bias, int slopeLimiter,
                                                                     int side,
                                                              HyperelasticMaterialParameters material1Parameters, HyperelasticMaterialParameters material2Parameters)
{
    vector<double> slopeVector = SlopeLimiters::computeSlopeVector(topStateVector, middleStateVector, bottomStateVector, bias, slopeLimiter, material1Parameters,
                                                                   material2Parameters);
    vector<double> topExtrapolatedValue = VectorAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                         VectorAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> bottomExtrapolatedValue = VectorAlgebra::addVectors(middleStateVector.computeConservedVariableVector(material1Parameters, material2Parameters),
                                                                       VectorAlgebra::multiplyVector(0.5, slopeVector));

    vector<double> topFluxVector = ElasticReducedStateVector::computeYFluxVector(topExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> bottomFluxVector = ElasticReducedStateVector::computeYFluxVector(bottomExtrapolatedValue, material1Parameters, material2Parameters);
    vector<double> evolutionVector = Solvers::computeEvolutionVector(topFluxVector, bottomFluxVector, cellSpacing, timeStep);

    return evolveStateByHalfTimeStepReduced(topExtrapolatedValue, bottomExtrapolatedValue, evolutionVector, side, material1Parameters, material2Parameters);
}

ElasticStateVector ElasticSolvers::evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue, vector<double> evolutionVector,
                                                             int side, HyperelasticMaterialParameters materialParameters)
{
    ElasticStateVector evolvedStateVector;

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

ElasticMultiphysicsStateVector ElasticSolvers::evolveStateByHalfTimeStep(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue,
                                                                         vector<double> evolutionVector, int side, HyperelasticMaterialParameters material1Parameters,
                                                                         HyperelasticMaterialParameters material2Parameters)
{
    ElasticMultiphysicsStateVector evolvedStateVector;

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

ElasticReducedStateVector ElasticSolvers::evolveStateByHalfTimeStepReduced(vector<double> leftExtrapolatedValue, vector<double> rightExtrapolatedValue,
                                                                           vector<double> evolutionVector, int side, HyperelasticMaterialParameters material1Parameters,
                                                                           HyperelasticMaterialParameters material2Parameters)
{
    ElasticReducedStateVector evolvedStateVector;

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

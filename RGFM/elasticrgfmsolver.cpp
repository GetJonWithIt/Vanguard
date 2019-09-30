#include "elasticrgfmsolver.h"

ElasticRGFMSolver::ElasticRGFMSolver()
{
}

ElasticMultimaterialSystem ElasticRGFMSolver::applyTildeRGFMBoundaryConditions(ElasticMultimaterialSystem multimaterialSystem, HyperelasticMaterialParameters material1Parameters,
                                                                               HyperelasticMaterialParameters material2Parameters)
{
    vector<ElasticStateVector> material1Cells = multimaterialSystem.getMaterial1Cells();
    vector<ElasticStateVector> material2Cells = multimaterialSystem.getMaterial2Cells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    vector<ElasticStateVector> newMaterial1Cells = material1Cells;
    vector<ElasticStateVector> newMaterial2Cells = material2Cells;

    double levelSetValue = levelSetFunction[0];
    int cellCount = levelSetFunction.size();
    bool inGhostRegion = false;

    for (int i = 1; i < cellCount; i++)
    {
        if ((levelSetFunction[i] * levelSetValue) <= 0.0 && i > 2 && i < (cellCount - 3))
        {
            for (int j = 0; j < 4; j++)
            {
                if (!inGhostRegion)
                {
                    ElasticStateVector riemannProblemSolution = ElasticHLLCSolver::solveTildeX(material1Cells[i - 1], material2Cells[i], material1Parameters, material2Parameters);

                    double material1GhostEntropy = riemannProblemSolution.getEntropy();
                    vector<vector<double> > material1GhostDistortionTensor = riemannProblemSolution.getDistortionTensor();

                    double material1GhostXVelocity = riemannProblemSolution.getXVelocity();
                    double material1GhostYVelocity = riemannProblemSolution.getYVelocity();
                    double material1GhostZVelocity = riemannProblemSolution.getZVelocity();

                    double material1GhostDensity = ElasticHLLCSolver::computeLeftTildeRegionDensity(material1Cells[i - 1], material2Cells[i], material1Parameters, material2Parameters);

                    newMaterial1Cells[i + j - 1].setDensity(material1GhostDensity);

                    newMaterial1Cells[i + j - 1].setXVelocity(material1GhostXVelocity);
                    newMaterial1Cells[i + j - 1].setYVelocity(material1GhostYVelocity);
                    newMaterial1Cells[i + j - 1].setZVelocity(material1GhostZVelocity);

                    newMaterial1Cells[i + j - 1].setDistortionTensor(material1GhostDistortionTensor);
                    newMaterial1Cells[i + j - 1].setEntropy(material1GhostEntropy);

                    double material2GhostDensity = ElasticHLLCSolver::computeRightTildeRegionDensity(material1Cells[i - 1], material2Cells[i], material1Parameters, material2Parameters);

                    newMaterial2Cells[i - j].setDensity(material2GhostDensity);

                    newMaterial2Cells[i - j].setXVelocity(material1GhostXVelocity);
                    newMaterial2Cells[i - j].setYVelocity(material1GhostYVelocity);
                    newMaterial2Cells[i - j].setZVelocity(material1GhostZVelocity);

                    newMaterial2Cells[i - j].setDistortionTensor(material1GhostDistortionTensor);
                    newMaterial2Cells[i - j].setEntropy(material1GhostEntropy);
                }
                else
                {
                    ElasticStateVector riemannProblemSolution = ElasticHLLCSolver::solveTildeX(material2Cells[i - 1], material1Cells[i], material2Parameters, material1Parameters);

                    double material2GhostEntropy = riemannProblemSolution.getEntropy();
                    vector<vector<double> > material2GhostDistortionTensor = riemannProblemSolution.getDistortionTensor();

                    double material2GhostXVelocity = riemannProblemSolution.getXVelocity();
                    double material2GhostYVelocity = riemannProblemSolution.getYVelocity();
                    double material2GhostZVelocity = riemannProblemSolution.getZVelocity();

                    double material2GhostDensity = ElasticHLLCSolver::computeLeftTildeRegionDensity(material2Cells[i - 1], material1Cells[i], material2Parameters, material1Parameters);

                    newMaterial2Cells[i + j - 1].setDensity(material2GhostDensity);

                    newMaterial2Cells[i + j - 1].setXVelocity(material2GhostXVelocity);
                    newMaterial2Cells[i + j - 1].setYVelocity(material2GhostYVelocity);
                    newMaterial2Cells[i + j - 1].setZVelocity(material2GhostZVelocity);

                    newMaterial2Cells[i + j - 1].setDistortionTensor(material2GhostDistortionTensor);
                    newMaterial2Cells[i + j - 1].setEntropy(material2GhostEntropy);

                    double material1GhostDensity = ElasticHLLCSolver::computeRightTildeRegionDensity(material2Cells[i - 1], material1Cells[i], material2Parameters, material1Parameters);

                    newMaterial1Cells[i - j].setDensity(material1GhostDensity);

                    newMaterial1Cells[i - j].setXVelocity(material2GhostXVelocity);
                    newMaterial1Cells[i - j].setYVelocity(material2GhostYVelocity);
                    newMaterial1Cells[i - j].setZVelocity(material2GhostZVelocity);

                    newMaterial1Cells[i - j].setDistortionTensor(material2GhostDistortionTensor);
                    newMaterial1Cells[i - j].setEntropy(material2GhostEntropy);
                }
            }

            inGhostRegion = !inGhostRegion;
        }

        levelSetValue = levelSetFunction[i];
    }

    return ElasticMultimaterialSystem(newMaterial1Cells, newMaterial2Cells, levelSetFunction);
}

ElasticMultimaterialSystem ElasticRGFMSolver::applyStarRGFMBoundaryConditions(ElasticMultimaterialSystem multimaterialSystem, HyperelasticMaterialParameters material1Parameters,
                                                                              HyperelasticMaterialParameters material2Parameters)
{
    vector<ElasticStateVector> material1Cells = multimaterialSystem.getMaterial1Cells();
    vector<ElasticStateVector> material2Cells = multimaterialSystem.getMaterial2Cells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    vector<ElasticStateVector> newMaterial1Cells = material1Cells;
    vector<ElasticStateVector> newMaterial2Cells = material2Cells;

    double levelSetValue = levelSetFunction[0];
    int cellCount = levelSetFunction.size();
    bool inGhostRegion = false;

    for (int i = 1; i < cellCount; i++)
    {
        if ((levelSetFunction[i] * levelSetValue) <= 0.0 && i > 2 && i < (cellCount - 3))
        {
            for (int j = 0; j < 4; j++)
            {
                if (!inGhostRegion)
                {
                    ElasticStateVector riemannProblemSolution = ElasticHLLCSolver::solveStarX(material1Cells[i - 1], material2Cells[i], material1Parameters, material2Parameters);

                    double material1GhostEntropy = riemannProblemSolution.getEntropy();
                    vector<vector<double> > material1GhostDistortionTensor = riemannProblemSolution.getDistortionTensor();

                    double material1GhostXVelocity = riemannProblemSolution.getXVelocity();
                    double material1GhostYVelocity = riemannProblemSolution.getYVelocity();
                    double material1GhostZVelocity = riemannProblemSolution.getZVelocity();

                    double material1GhostDensity = ElasticHLLCSolver::computeLeftTildeRegionDensity(material1Cells[i - 1], material2Cells[i], material1Parameters, material2Parameters);

                    newMaterial1Cells[i + j - 1].setDensity(material1GhostDensity);

                    newMaterial1Cells[i + j - 1].setXVelocity(material1GhostXVelocity);
                    newMaterial1Cells[i + j - 1].setYVelocity(material1GhostYVelocity);
                    newMaterial1Cells[i + j - 1].setZVelocity(material1GhostZVelocity);

                    newMaterial1Cells[i + j - 1].setDistortionTensor(material1GhostDistortionTensor);
                    newMaterial1Cells[i + j - 1].setEntropy(material1GhostEntropy);

                    double material2GhostDensity = ElasticHLLCSolver::computeRightTildeRegionDensity(material1Cells[i - 1], material2Cells[i], material1Parameters, material2Parameters);

                    newMaterial2Cells[i - j].setDensity(material2GhostDensity);

                    newMaterial2Cells[i - j].setXVelocity(material1GhostXVelocity);
                    newMaterial2Cells[i - j].setYVelocity(material1GhostYVelocity);
                    newMaterial2Cells[i - j].setZVelocity(material1GhostZVelocity);

                    newMaterial2Cells[i - j].setDistortionTensor(material1GhostDistortionTensor);
                    newMaterial2Cells[i - j].setEntropy(material1GhostEntropy);
                }
                else
                {
                    ElasticStateVector riemannProblemSolution = ElasticHLLCSolver::solveStarX(material2Cells[i - 1], material1Cells[i], material2Parameters, material1Parameters);

                    double material2GhostEntropy = riemannProblemSolution.getEntropy();
                    vector<vector<double> > material2GhostDistortionTensor = riemannProblemSolution.getDistortionTensor();

                    double material2GhostXVelocity = riemannProblemSolution.getXVelocity();
                    double material2GhostYVelocity = riemannProblemSolution.getYVelocity();
                    double material2GhostZVelocity = riemannProblemSolution.getZVelocity();

                    double material2GhostDensity = ElasticHLLCSolver::computeLeftTildeRegionDensity(material2Cells[i - 1], material1Cells[i], material2Parameters, material1Parameters);

                    newMaterial2Cells[i + j - 1].setDensity(material2GhostDensity);

                    newMaterial2Cells[i + j - 1].setXVelocity(material2GhostXVelocity);
                    newMaterial2Cells[i + j - 1].setYVelocity(material2GhostYVelocity);
                    newMaterial2Cells[i + j - 1].setZVelocity(material2GhostZVelocity);

                    newMaterial2Cells[i + j - 1].setDistortionTensor(material2GhostDistortionTensor);
                    newMaterial2Cells[i + j - 1].setEntropy(material2GhostEntropy);

                    double material1GhostDensity = ElasticHLLCSolver::computeRightTildeRegionDensity(material2Cells[i - 1], material1Cells[i], material2Parameters, material1Parameters);

                    newMaterial1Cells[i - j].setDensity(material1GhostDensity);

                    newMaterial1Cells[i - j].setXVelocity(material2GhostXVelocity);
                    newMaterial1Cells[i - j].setYVelocity(material2GhostYVelocity);
                    newMaterial1Cells[i - j].setZVelocity(material2GhostZVelocity);

                    newMaterial1Cells[i - j].setDistortionTensor(material2GhostDistortionTensor);
                    newMaterial1Cells[i - j].setEntropy(material2GhostEntropy);
                }
            }

            inGhostRegion = !inGhostRegion;
        }

        levelSetValue = levelSetFunction[i];
    }

    return ElasticMultimaterialSystem(newMaterial1Cells, newMaterial2Cells, levelSetFunction);
}

vector<double> ElasticRGFMSolver::updateLevelSetFunction(vector<double> levelSetFunction, double cellSpacing, double timeStep, vector<ElasticStateVector> material1Cells,
                                                         vector<ElasticStateVector> material2Cells)
{
    int cellCount = levelSetFunction.size();
    vector<double> newLevelSetFunction = levelSetFunction;
    vector<double> levelSetFunctionWithBoundary = RGFMSolver::insertBoundaryCells(levelSetFunction, 2);

    double upwindApproximation;
    double velocity;

    double levelSetValue = levelSetFunction[0];
    double inGhostRegion = false;

    for (int i = 0; i < cellCount; i++)
    {
        if ((levelSetValue * levelSetFunction[i]) <= 0.0)
        {
            inGhostRegion = !inGhostRegion;
        }

        if (!inGhostRegion)
        {
            velocity = material1Cells[i].getXVelocity();
        }
        else
        {
            velocity = material2Cells[i].getXVelocity();
        }

        if (velocity < 0.0)
        {
            upwindApproximation = (-levelSetFunctionWithBoundary[i + 4] + (4.0 * levelSetFunctionWithBoundary[i + 3]) - (3.0 * levelSetFunctionWithBoundary[i + 2])) / (2.0 * cellSpacing);
        }
        else
        {
            upwindApproximation = ((3.0 * levelSetFunctionWithBoundary[i + 2]) - (4.0 * levelSetFunctionWithBoundary[i + 1]) + levelSetFunctionWithBoundary[i]) / (2.0 * cellSpacing);
        }
        newLevelSetFunction[i] = levelSetFunction[i] - (timeStep * (velocity * upwindApproximation));

        levelSetValue = levelSetFunction[i];
    }

    return newLevelSetFunction;
}

ElasticMultimaterialSystem ElasticRGFMSolver::solveTilde(ElasticMultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                         int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, HyperelasticMaterialParameters material1Parameters,
                                                         HyperelasticMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;

    vector<ElasticStateVector> material1Cells = multimaterialSystem.getMaterial1Cells();
    vector<ElasticStateVector> material2Cells = multimaterialSystem.getMaterial2Cells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    vector<ElasticStateVector> currentMaterial1Cells = material1Cells;
    vector<ElasticStateVector> currentMaterial2Cells = material2Cells;

    while (currentTime < finalTime)
    {
        ElasticMultimaterialSystem newMultimaterialSystem = applyTildeRGFMBoundaryConditions(ElasticMultimaterialSystem(currentMaterial1Cells, currentMaterial2Cells, levelSetFunction),
                                                                                             material1Parameters, material2Parameters);

        currentMaterial1Cells = newMultimaterialSystem.getMaterial1Cells();
        currentMaterial2Cells = newMultimaterialSystem.getMaterial2Cells();

        vector<ElasticStateVector> currentMaterial1CellsWithBoundary = ElasticSolvers::insertBoundaryCells(currentMaterial1Cells, 2);
        vector<ElasticStateVector> currentMaterial2CellsWithBoundary = ElasticSolvers::insertBoundaryCells(currentMaterial2Cells, 2);

        double timeStep = min(ElasticSolvers::computeStableTimeStep(currentMaterial1CellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters),
                              ElasticSolvers::computeStableTimeStep(currentMaterial2CellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material2Parameters));

        levelSetFunction = updateLevelSetFunction(levelSetFunction, cellSpacing, timeStep, currentMaterial1Cells, currentMaterial2Cells);

        ElasticSecondOrderSolver::computeSLICTimeStep(currentMaterial1Cells, currentMaterial1CellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters);
        ElasticSecondOrderSolver::computeSLICTimeStep(currentMaterial2Cells, currentMaterial2CellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material2Parameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return ElasticMultimaterialSystem(currentMaterial1Cells, currentMaterial2Cells, levelSetFunction);
}

ElasticMultimaterialSystem ElasticRGFMSolver::solveStar(ElasticMultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                        int slopeLimiter, int subcyclingIterations, int reinitialisationFrequency, HyperelasticMaterialParameters material1Parameters,
                                                        HyperelasticMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;

    vector<ElasticStateVector> material1Cells = multimaterialSystem.getMaterial1Cells();
    vector<ElasticStateVector> material2Cells = multimaterialSystem.getMaterial2Cells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    vector<ElasticStateVector> currentMaterial1Cells = material1Cells;
    vector<ElasticStateVector> currentMaterial2Cells = material2Cells;

    while (currentTime < finalTime)
    {
        ElasticMultimaterialSystem newMultimaterialSystem = applyStarRGFMBoundaryConditions(ElasticMultimaterialSystem(currentMaterial1Cells, currentMaterial2Cells, levelSetFunction),
                                                                                            material1Parameters, material2Parameters);

        currentMaterial1Cells = newMultimaterialSystem.getMaterial1Cells();
        currentMaterial2Cells = newMultimaterialSystem.getMaterial2Cells();

        vector<ElasticStateVector> currentMaterial1CellsWithBoundary = ElasticSolvers::insertBoundaryCells(currentMaterial1Cells, 2);
        vector<ElasticStateVector> currentMaterial2CellsWithBoundary = ElasticSolvers::insertBoundaryCells(currentMaterial2Cells, 2);

        double timeStep = min(ElasticSolvers::computeStableTimeStep(currentMaterial1CellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters),
                              ElasticSolvers::computeStableTimeStep(currentMaterial2CellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material2Parameters));

        levelSetFunction = updateLevelSetFunction(levelSetFunction, cellSpacing, timeStep, currentMaterial1Cells, currentMaterial2Cells);

        ElasticSecondOrderSolver::computeSLICTimeStep(currentMaterial1Cells, currentMaterial1CellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters);
        ElasticSecondOrderSolver::computeSLICTimeStep(currentMaterial2Cells, currentMaterial2CellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material2Parameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return ElasticMultimaterialSystem(currentMaterial1Cells, currentMaterial2Cells, levelSetFunction);
}

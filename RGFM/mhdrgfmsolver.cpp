#include "mhdrgfmsolver.h"

MHDRGFMSolver::MHDRGFMSolver()
{
}

MHDMultimaterialSystem MHDRGFMSolver::applyRGFMBoundaryConditions(MHDMultimaterialSystem multimaterialSystem, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    vector<MHDStateVector> material1Cells = multimaterialSystem.getMaterial1Cells();
    vector<MHDStateVector> material2Cells = multimaterialSystem.getMaterial2Cells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    vector<MHDStateVector> newMaterial1Cells = material1Cells;
    vector<MHDStateVector> newMaterial2Cells = material2Cells;

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
                    MHDStateVector riemannProblemSolution = MHDHLLCSolver::solveX(material1Cells[i - 1], material2Cells[i], material1Parameters, material2Parameters);

                    double material1GhostXMagneticField = riemannProblemSolution.getXMagneticField();
                    double material1GhostYMagneticField = riemannProblemSolution.getYMagneticField();
                    double material1GhostZMagneticField = riemannProblemSolution.getZMagneticField();

                    double material1GhostPressure = riemannProblemSolution.getPressure();
                    if (material1GhostPressure < pow(10.0, -8.0))
                    {
                        material1GhostPressure = pow(10.0, -8.0);
                    }

                    double material1GhostXVelocity = riemannProblemSolution.getXVelocity();
                    double material1GhostYVelocity = riemannProblemSolution.getYVelocity();
                    double material1GhostZVelocity = riemannProblemSolution.getZVelocity();

                    double material1GhostDensity = MHDHLLCSolver::computeLeftStarRegionDensity(material1Cells[i - 1], material2Cells[i], material1Parameters, material2Parameters);

                    newMaterial1Cells[i + j - 1].setDensity(material1GhostDensity);

                    newMaterial1Cells[i + j - 1].setXVelocity(material1GhostXVelocity);
                    newMaterial1Cells[i + j - 1].setYVelocity(material1GhostYVelocity);
                    newMaterial1Cells[i + j - 1].setZVelocity(material1GhostZVelocity);

                    newMaterial1Cells[i + j - 1].setPressure(material1GhostPressure);

                    newMaterial1Cells[i + j - 1].setXMagneticField(material1GhostXMagneticField);
                    newMaterial1Cells[i + j - 1].setYMagneticField(material1GhostYMagneticField);
                    newMaterial1Cells[i + j - 1].setZMagneticField(material1GhostZMagneticField);

                    double material2GhostDensity = MHDHLLCSolver::computeRightStarRegionDensity(material1Cells[i - 1], material2Cells[i], material1Parameters, material2Parameters);

                    newMaterial2Cells[i - j].setDensity(material2GhostDensity);

                    newMaterial2Cells[i - j].setXVelocity(material1GhostXVelocity);
                    newMaterial2Cells[i - j].setYVelocity(material1GhostYVelocity);
                    newMaterial2Cells[i - j].setZVelocity(material1GhostZVelocity);

                    newMaterial2Cells[i - j].setPressure(material1GhostPressure);

                    newMaterial2Cells[i - j].setXMagneticField(material1GhostXMagneticField);
                    newMaterial2Cells[i - j].setYMagneticField(material1GhostYMagneticField);
                    newMaterial2Cells[i - j].setZMagneticField(material1GhostZMagneticField);
                }
                else
                {
                    MHDStateVector riemannProblemSolution = MHDHLLCSolver::solveX(material2Cells[i - 1], material1Cells[i], material2Parameters, material1Parameters);

                    double material2GhostXMagneticField = riemannProblemSolution.getXMagneticField();
                    double material2GhostYMagneticField = riemannProblemSolution.getYMagneticField();
                    double material2GhostZMagneticField = riemannProblemSolution.getZMagneticField();

                    double material2GhostPressure = riemannProblemSolution.getPressure();
                    if (material2GhostPressure < pow(10.0, -8.0))
                    {
                        material2GhostPressure = pow(10.0, -8.0);
                    }

                    double material2GhostXVelocity = riemannProblemSolution.getXVelocity();
                    double material2GhostYVelocity = riemannProblemSolution.getYVelocity();
                    double material2GhostZVelocity = riemannProblemSolution.getZVelocity();

                    double material2GhostDensity = MHDHLLCSolver::computeLeftStarRegionDensity(material2Cells[i - 1], material1Cells[i], material2Parameters, material1Parameters);

                    newMaterial2Cells[i + j - 1].setDensity(material2GhostDensity);

                    newMaterial2Cells[i + j - 1].setXVelocity(material2GhostXVelocity);
                    newMaterial2Cells[i + j - 1].setYVelocity(material2GhostYVelocity);
                    newMaterial2Cells[i + j - 1].setZVelocity(material2GhostZVelocity);

                    newMaterial2Cells[i + j - 1].setPressure(material2GhostPressure);

                    newMaterial2Cells[i + j - 1].setXMagneticField(material2GhostXMagneticField);
                    newMaterial2Cells[i + j - 1].setYMagneticField(material2GhostYMagneticField);
                    newMaterial2Cells[i + j - 1].setZMagneticField(material2GhostZMagneticField);

                    double material1GhostDensity = MHDHLLCSolver::computeRightStarRegionDensity(material2Cells[i - 1], material1Cells[i], material2Parameters, material1Parameters);

                    newMaterial1Cells[i - j].setDensity(material1GhostDensity);

                    newMaterial1Cells[i - j].setXVelocity(material2GhostXVelocity);
                    newMaterial1Cells[i - j].setYVelocity(material2GhostYVelocity);
                    newMaterial1Cells[i - j].setZVelocity(material2GhostZVelocity);

                    newMaterial1Cells[i - j].setPressure(material2GhostPressure);

                    newMaterial1Cells[i - j].setXMagneticField(material2GhostXMagneticField);
                    newMaterial1Cells[i - j].setYMagneticField(material2GhostYMagneticField);
                    newMaterial1Cells[i - j].setZMagneticField(material2GhostZMagneticField);
                }
            }

            inGhostRegion = !inGhostRegion;
        }

        levelSetValue = levelSetFunction[i];
    }

    return MHDMultimaterialSystem(newMaterial1Cells, newMaterial2Cells, levelSetFunction);
}

vector<double> MHDRGFMSolver::updateLevelSetFunction(vector<double> levelSetFunction, double cellSpacing, double timeStep, vector<MHDStateVector> material1Cells,
                                                     vector<MHDStateVector> material2Cells)
{
    int cellCount = levelSetFunction.size();
    vector<double> newLevelSetFunction = levelSetFunction;
    vector<double> levelSetFunctionWithBoundary = RGFMSolver::insertBoundaryCells(levelSetFunction, 2);

    double upwindApproximation;
    double velocity;

    double levelSetValue = levelSetFunction[0];
    bool inGhostRegion = false;

    for (int i = 0; i < cellCount; i++)
    {
        /*
        if ((levelSetValue * levelSetFunction[i]) <= 0.0 && (levelSetValue * levelSetFunction[i - 1]) > 0.0)
        {
            inGhostRegion = !inGhostRegion;
        }
        */
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

        if (velocity < 0)
        {
            upwindApproximation = (-levelSetFunctionWithBoundary[i + 4] + (4 * levelSetFunctionWithBoundary[i + 3]) - (3 * levelSetFunctionWithBoundary[i + 2])) / (2 * cellSpacing);
        }
        else
        {
            upwindApproximation = ((3 * levelSetFunctionWithBoundary[i + 2]) - (4 * levelSetFunctionWithBoundary[i + 1]) + levelSetFunctionWithBoundary[i]) / (2 * cellSpacing);
        }
        newLevelSetFunction[i] = levelSetFunction[i] - (timeStep * (velocity * upwindApproximation));

        levelSetValue = levelSetFunction[i];
    }

    return newLevelSetFunction;
}

MHDMultimaterialSystem MHDRGFMSolver::solve(MHDMultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                            int subcyclingIterations, int reinitialisationFrequency, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;

    vector<MHDStateVector> material1Cells = multimaterialSystem.getMaterial1Cells();
    vector<MHDStateVector> material2Cells = multimaterialSystem.getMaterial2Cells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    vector<MHDStateVector> currentMaterial1Cells = material1Cells;
    vector<MHDStateVector> currentMaterial2Cells = material2Cells;

    while (currentTime < finalTime)
    {
        MHDMultimaterialSystem newMultimaterialSystem = applyRGFMBoundaryConditions(MHDMultimaterialSystem(currentMaterial1Cells, currentMaterial2Cells, levelSetFunction),
                                                                                    material1Parameters, material2Parameters);

        currentMaterial1Cells = newMultimaterialSystem.getMaterial1Cells();
        currentMaterial2Cells = newMultimaterialSystem.getMaterial2Cells();

        vector<MHDStateVector> currentMaterial1CellsWithBoundary = MHDSolvers::insertBoundaryCells(currentMaterial1Cells, 2);
        vector<MHDStateVector> currentMaterial2CellsWithBoundary = MHDSolvers::insertBoundaryCells(currentMaterial2Cells, 2);

        double timeStep = min(MHDSolvers::computeStableTimeStep(currentMaterial1CellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters),
                              MHDSolvers::computeStableTimeStep(currentMaterial2CellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material2Parameters));

        levelSetFunction = updateLevelSetFunction(levelSetFunction, cellSpacing, timeStep, currentMaterial1Cells, currentMaterial2Cells);

        MHDSecondOrderSolver::computeSLICTimeStep(currentMaterial1Cells, currentMaterial1CellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters);
        MHDSecondOrderSolver::computeSLICTimeStep(currentMaterial2Cells, currentMaterial2CellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material2Parameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return MHDMultimaterialSystem(currentMaterial1Cells, currentMaterial2Cells, levelSetFunction);
}

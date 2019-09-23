#include "rgfmsolver.h"

RGFMSolver::RGFMSolver()
{
}

vector<double> RGFMSolver::insertBoundaryCells(vector<double> & levelSetFunction, int boundarySize)
{
    int cellCount = levelSetFunction.size();
    vector<double> levelSetFunctionWithBoundary(cellCount + (2 * boundarySize));

    if (boundarySize == 1)
    {
        levelSetFunctionWithBoundary[0] = levelSetFunction[0];
        levelSetFunctionWithBoundary[cellCount + 1] = levelSetFunction[cellCount - 1];
    }
    else if (boundarySize == 2)
    {
        levelSetFunctionWithBoundary[0] = levelSetFunction[1];
        levelSetFunctionWithBoundary[1] = levelSetFunction[0];

        levelSetFunctionWithBoundary[cellCount + 2] = levelSetFunction[cellCount - 1];
        levelSetFunctionWithBoundary[cellCount + 3] = levelSetFunction[cellCount - 2];
    }

#pragma omp parallel for
    for (int i = 0; i < cellCount; i++)
    {
        levelSetFunctionWithBoundary[i + boundarySize] = levelSetFunction[i];
    }

    return levelSetFunctionWithBoundary;
}

vector<vector<double> > RGFMSolver::insertBoundaryCells2D(vector<vector<double> > & levelSetFunction, int boundarySize)
{
    int rowCount = levelSetFunction.size();
    int columnCount = levelSetFunction[0].size();
    vector<vector<double> > levelSetFunctionWithBoundary(rowCount + (2 * boundarySize), vector<double>(columnCount + (2 * boundarySize)));

    if (boundarySize == 1)
    {
#pragma omp parallel for
        for (int i = 0; i < rowCount; i++)
        {
            levelSetFunctionWithBoundary[i + 1][0] = levelSetFunction[i][0];
            levelSetFunctionWithBoundary[i + 1][columnCount + 1] = levelSetFunction[i][columnCount - 1];
        }

#pragma omp parallel for
        for (int i = 0; i < columnCount; i++)
        {
            levelSetFunctionWithBoundary[0][i + 1] = levelSetFunction[0][i];
            levelSetFunctionWithBoundary[rowCount + 1][i + 1] = levelSetFunction[rowCount - 1][i];
        }
    }
    else if (boundarySize == 2)
    {
#pragma omp parallel for
        for (int i = 0; i < rowCount; i++)
        {
            levelSetFunctionWithBoundary[i + 2][0] = levelSetFunction[i][1];
            levelSetFunctionWithBoundary[i + 2][1] = levelSetFunction[i][0];

            levelSetFunctionWithBoundary[i + 2][columnCount + 2] = levelSetFunction[i][columnCount - 1];
            levelSetFunctionWithBoundary[i + 2][columnCount + 3] = levelSetFunction[i][columnCount - 2];
        }

#pragma omp parallel for
        for (int i = 0; i < columnCount; i++)
        {
            levelSetFunctionWithBoundary[0][i + 2] = levelSetFunction[1][i];
            levelSetFunctionWithBoundary[1][i + 2] = levelSetFunction[0][i];

            levelSetFunctionWithBoundary[rowCount + 2][i + 2] = levelSetFunction[rowCount - 1][i];
            levelSetFunctionWithBoundary[rowCount + 3][i + 2] = levelSetFunction[rowCount - 2][i];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            levelSetFunctionWithBoundary[i + boundarySize][j + boundarySize] = levelSetFunction[i][j];
        }
    }

    return levelSetFunctionWithBoundary;
}

MultimaterialSystem RGFMSolver::applyRGFMBoundaryConditionsExact(MultimaterialSystem multimaterialSystem, EulerMaterialParameters material1Parameters,
                                                                 EulerMaterialParameters material2Parameters)
{
    vector<EulerStateVector> material1Cells = multimaterialSystem.getMaterial1Cells();
    vector<EulerStateVector> material2Cells = multimaterialSystem.getMaterial2Cells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    vector<EulerStateVector> newMaterial1Cells = material1Cells;
    vector<EulerStateVector> newMaterial2Cells = material2Cells;

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
                    EulerStateVector riemannProblemSolution = ExactSolver::solveX(0.0, 1.0, 0.0, material1Cells[i - 1], material2Cells[i], material1Parameters, material2Parameters);

                    double material1GhostPressure = riemannProblemSolution.getPressure();
                    double material1GhostVelocity = riemannProblemSolution.getXVelocity();
                    double material1GhostDensity = ExactSolver::computeStarRegionDensity(material1GhostPressure, material1Cells[i - 1], material1Parameters);

                    newMaterial1Cells[i + j - 1].setDensity(material1GhostDensity);
                    newMaterial1Cells[i + j - 1].setXVelocity(material1GhostVelocity);
                    newMaterial1Cells[i + j - 1].setPressure(material1GhostPressure);

                    double material2GhostDensity = ExactSolver::computeStarRegionDensity(material1GhostPressure, material2Cells[i], material2Parameters);

                    newMaterial2Cells[i - j].setDensity(material2GhostDensity);
                    newMaterial2Cells[i - j].setXVelocity(material1GhostVelocity);
                    newMaterial2Cells[i - j].setPressure(material1GhostPressure);
                }
                else
                {
                    EulerStateVector riemannProblemSolution = ExactSolver::solveX(0.0, 1.0, 0.0, material2Cells[i - 1], material1Cells[i], material2Parameters, material1Parameters);

                    double material2GhostPressure = riemannProblemSolution.getPressure();
                    double material2GhostVelocity = riemannProblemSolution.getXVelocity();
                    double material2GhostDensity = ExactSolver::computeStarRegionDensity(material2GhostPressure, material2Cells[i - 1], material2Parameters);

                    newMaterial2Cells[i + j - 1].setDensity(material2GhostDensity);
                    newMaterial2Cells[i + j - 1].setXVelocity(material2GhostVelocity);
                    newMaterial2Cells[i + j - 1].setPressure(material2GhostPressure);

                    double material1GhostDensity = ExactSolver::computeStarRegionDensity(material2GhostPressure, material1Cells[i], material1Parameters);

                    newMaterial1Cells[i - j].setDensity(material1GhostDensity);
                    newMaterial1Cells[i - j].setXVelocity(material2GhostVelocity);
                    newMaterial1Cells[i - j].setPressure(material2GhostPressure);
                }
            }

            inGhostRegion = !inGhostRegion;
        }

        levelSetValue = levelSetFunction[i];
    }

    return MultimaterialSystem(newMaterial1Cells, newMaterial2Cells, levelSetFunction);
}

MultimaterialSystem RGFMSolver::applyRGFMBoundaryConditionsExact2D(MultimaterialSystem multimaterialSystem, EulerMaterialParameters material1Parameters,
                                                                   EulerMaterialParameters material2Parameters)
{
    vector<vector<EulerStateVector> > material1Cells = multimaterialSystem.getMaterial1Cells2D();
    vector<vector<EulerStateVector> > material2Cells = multimaterialSystem.getMaterial2Cells2D();
    vector<vector<double> > levelSetFunction = multimaterialSystem.getLevelSetFunction2D();

    vector<vector<EulerStateVector> > newMaterial1Cells = material1Cells;
    vector<vector<EulerStateVector> > newMaterial2Cells = material2Cells;

    int rowCount = levelSetFunction.size();
    int columnCount = levelSetFunction[0].size();

    double levelSetValue;
    bool inGhostRegion;

    for (int i = 0; i < rowCount; i++)
    {
        levelSetValue = levelSetFunction[i][0];
        inGhostRegion = true;

        for (int j = 1; j < columnCount; j++)
        {
            if ((levelSetFunction[i][j] * levelSetValue) <= 0.0 && j > 2 && j < (columnCount - 3))
            {
                for (int k = 0; k < 4; k++)
                {
                    if (!inGhostRegion)
                    {
                        EulerStateVector riemannProblemSolution = ExactSolver::solveX(0.0, 1.0, 0.0, material1Cells[i][j - 1], material2Cells[i][j], material1Parameters, material2Parameters);

                        double material1GhostPressure = riemannProblemSolution.getPressure();
                        double material1GhostVelocity = riemannProblemSolution.getXVelocity();
                        double material1GhostDensity = ExactSolver::computeStarRegionDensity(material1GhostPressure, material1Cells[i][j - 1], material1Parameters);

                        newMaterial1Cells[i][j + k - 1].setDensity(material1GhostDensity);
                        newMaterial1Cells[i][j + k - 1].setXVelocity(material1GhostVelocity);
                        newMaterial1Cells[i][j + k - 1].setPressure(material1GhostPressure);

                        double material2GhostDensity = ExactSolver::computeStarRegionDensity(material1GhostPressure, material2Cells[i][j], material2Parameters);

                        newMaterial2Cells[i][j - k].setDensity(material2GhostDensity);
                        newMaterial2Cells[i][j - k].setXVelocity(material1GhostVelocity);
                        newMaterial2Cells[i][j - k].setPressure(material1GhostPressure);
                    }
                    else
                    {
                        EulerStateVector riemannProblemSolution = ExactSolver::solveX(0.0, 1.0, 0.0, material2Cells[i][j - 1], material1Cells[i][j], material2Parameters, material1Parameters);

                        double material2GhostPressure = riemannProblemSolution.getPressure();
                        double material2GhostVelocity = riemannProblemSolution.getXVelocity();
                        double material2GhostDensity = ExactSolver::computeStarRegionDensity(material2GhostPressure, material2Cells[i][j - 1], material2Parameters);

                        newMaterial2Cells[i][j + k - 1].setDensity(material2GhostDensity);
                        newMaterial2Cells[i][j + k - 1].setXVelocity(material2GhostVelocity);
                        newMaterial2Cells[i][j + k - 1].setPressure(material2GhostPressure);

                        double material1GhostDensity = ExactSolver::computeStarRegionDensity(material2GhostPressure, material1Cells[i][j], material1Parameters);

                        newMaterial1Cells[i][j - k].setDensity(material1GhostDensity);
                        newMaterial1Cells[i][j - k].setXVelocity(material2GhostVelocity);
                        newMaterial1Cells[i][j - k].setPressure(material2GhostPressure);
                    }
                }

                inGhostRegion = !inGhostRegion;
            }

            levelSetValue = levelSetFunction[i][j];
        }
    }

    material1Cells = newMaterial1Cells;
    material2Cells = newMaterial2Cells;

    for (int i = 0; i < columnCount; i++)
    {
        levelSetValue = levelSetFunction[0][i];
        inGhostRegion = true;

        for (int j = 1; j < rowCount; j++)
        {
            if ((levelSetFunction[j][i] * levelSetValue) <= 0.0 && j > 2 && j < (rowCount - 3))
            {
                for (int k = 0; k < 4; k++)
                {
                    if (!inGhostRegion)
                    {
                        EulerStateVector riemannProblemSolution = ExactSolver::solveY(0.0, 1.0, 0.0, material1Cells[j - 1][i], material2Cells[j][i], material1Parameters, material2Parameters);

                        double material1GhostPressure = riemannProblemSolution.getPressure();
                        double material1GhostVelocity = riemannProblemSolution.getYVelocity();
                        double material1GhostDensity = ExactSolver::computeStarRegionDensity(material1GhostPressure, material1Cells[j - 1][i], material1Parameters);

                        newMaterial1Cells[j + k - 1][i].setDensity(material1GhostDensity);
                        newMaterial1Cells[j + k - 1][i].setYVelocity(material1GhostVelocity);
                        newMaterial1Cells[j + k - 1][i].setPressure(material1GhostPressure);

                        double material2GhostDensity = ExactSolver::computeStarRegionDensity(material1GhostPressure, material2Cells[j][i], material2Parameters);

                        newMaterial2Cells[j - k][i].setDensity(material2GhostDensity);
                        newMaterial2Cells[j - k][i].setYVelocity(material1GhostVelocity);
                        newMaterial2Cells[j - k][i].setPressure(material1GhostPressure);
                    }
                    else
                    {
                        EulerStateVector riemannProblemSolution = ExactSolver::solveY(0.0, 1.0, 0.0, material2Cells[j - 1][i], material1Cells[j][i], material2Parameters, material1Parameters);

                        double material2GhostPressure = riemannProblemSolution.getPressure();
                        double material2GhostVelocity = riemannProblemSolution.getYVelocity();
                        double material2GhostDensity = ExactSolver::computeStarRegionDensity(material2GhostPressure, material2Cells[j - 1][i], material2Parameters);

                        newMaterial2Cells[j + k - 1][i].setDensity(material2GhostDensity);
                        newMaterial2Cells[j + k - 1][i].setYVelocity(material2GhostVelocity);
                        newMaterial2Cells[j + k - 1][i].setPressure(material2GhostPressure);

                        double material1GhostDensity = ExactSolver::computeStarRegionDensity(material2GhostPressure, material1Cells[j][i], material1Parameters);

                        newMaterial1Cells[j - k][i].setDensity(material1GhostDensity);
                        newMaterial1Cells[j - k][i].setYVelocity(material2GhostVelocity);
                        newMaterial1Cells[j - k][i].setPressure(material2GhostPressure);
                    }
                }

                inGhostRegion = !inGhostRegion;
            }

            levelSetValue = levelSetFunction[j][i];
        }
    }

    return MultimaterialSystem(newMaterial1Cells, newMaterial2Cells, levelSetFunction);
}

MultimaterialSystem RGFMSolver::applyRGFMBoundaryConditionsHLLC(MultimaterialSystem multimaterialSystem, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    vector<EulerStateVector> material1Cells = multimaterialSystem.getMaterial1Cells();
    vector<EulerStateVector> material2Cells = multimaterialSystem.getMaterial2Cells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    vector<EulerStateVector> newMaterial1Cells = material1Cells;
    vector<EulerStateVector> newMaterial2Cells = material2Cells;

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
                    EulerStateVector riemannProblemSolution = HLLCSolver::solveX(material1Cells[i - 1], material2Cells[i], material1Parameters, material2Parameters);

                    double material1GhostPressure = riemannProblemSolution.getPressure();
                    double material1GhostVelocity = riemannProblemSolution.getXVelocity();
                    double material1GhostDensity = HLLCSolver::computeLeftStarRegionDensity(material1Cells[i - 1], material2Cells[i], material1Parameters, material2Parameters);

                    newMaterial1Cells[i + j - 1].setDensity(material1GhostDensity);
                    newMaterial1Cells[i + j - 1].setXVelocity(material1GhostVelocity);
                    newMaterial1Cells[i + j - 1].setPressure(material1GhostPressure);

                    double material2GhostDensity = HLLCSolver::computeRightStarRegionDensity(material1Cells[i - 1], material2Cells[i], material1Parameters, material2Parameters);

                    newMaterial2Cells[i - j].setDensity(material2GhostDensity);
                    newMaterial2Cells[i - j].setXVelocity(material1GhostVelocity);
                    newMaterial2Cells[i - j].setPressure(material1GhostPressure);
                }
                else
                {
                    EulerStateVector riemannProblemSolution = HLLCSolver::solveX(material2Cells[i - 1], material1Cells[i], material2Parameters, material1Parameters);

                    double material2GhostPressure = riemannProblemSolution.getPressure();
                    double material2GhostVelocity = riemannProblemSolution.getXVelocity();
                    double material2GhostDensity = HLLCSolver::computeLeftStarRegionDensity(material2Cells[i - 1], material1Cells[i], material2Parameters, material1Parameters);

                    newMaterial2Cells[i + j - 1].setDensity(material2GhostDensity);
                    newMaterial2Cells[i + j - 1].setXVelocity(material2GhostVelocity);
                    newMaterial2Cells[i + j - 1].setPressure(material2GhostPressure);

                    double material1GhostDensity = HLLCSolver::computeRightStarRegionDensity(material2Cells[i - 1], material1Cells[i], material2Parameters, material1Parameters);

                    newMaterial1Cells[i - j].setDensity(material1GhostDensity);
                    newMaterial1Cells[i - j].setXVelocity(material2GhostVelocity);
                    newMaterial1Cells[i - j].setPressure(material2GhostPressure);
                }
            }

            inGhostRegion = !inGhostRegion;
        }

        levelSetValue = levelSetFunction[i];
    }

    return MultimaterialSystem(newMaterial1Cells, newMaterial2Cells, levelSetFunction);
}

MultimaterialSystem RGFMSolver::applyRGFMBoundaryConditionsHLLC2D(MultimaterialSystem multimaterialSystem, EulerMaterialParameters material1Parameters,
                                                                  EulerMaterialParameters material2Parameters)
{
    vector<vector<EulerStateVector> > material1Cells = multimaterialSystem.getMaterial1Cells2D();
    vector<vector<EulerStateVector> > material2Cells = multimaterialSystem.getMaterial2Cells2D();
    vector<vector<double> > levelSetFunction = multimaterialSystem.getLevelSetFunction2D();

    vector<vector<EulerStateVector> > newMaterial1Cells = material1Cells;
    vector<vector<EulerStateVector> > newMaterial2Cells = material2Cells;

    int rowCount = levelSetFunction.size();
    int columnCount = levelSetFunction[0].size();

    double levelSetValue;
    bool inGhostRegion;

    for (int i = 0; i < rowCount; i++)
    {
        levelSetValue = levelSetFunction[i][0];
        inGhostRegion = true;

        for (int j = 1; j < columnCount; j++)
        {
            if ((levelSetFunction[i][j] * levelSetValue) <= 0.0 && j > 2 && j < (columnCount - 3))
            {
                for (int k = 0; k < 4; k++)
                {
                    if (!inGhostRegion)
                    {
                        EulerStateVector riemannProblemSolution = HLLCSolver::solveX(material1Cells[i][j - 1], material2Cells[i][j], material1Parameters, material2Parameters);

                        double material1GhostPressure = riemannProblemSolution.getPressure();
                        double material1GhostVelocity = riemannProblemSolution.getXVelocity();
                        double material1GhostDensity = HLLCSolver::computeLeftStarRegionDensity(material1Cells[i][j - 1], material2Cells[i][j], material1Parameters, material2Parameters);

                        newMaterial1Cells[i][j + k - 1].setDensity(material1GhostDensity);
                        newMaterial1Cells[i][j + k - 1].setXVelocity(material1GhostVelocity);
                        newMaterial1Cells[i][j + k - 1].setPressure(material1GhostPressure);

                        double material2GhostDensity = HLLCSolver::computeRightStarRegionDensity(material1Cells[i][j - 1], material2Cells[i][j], material1Parameters, material2Parameters);

                        newMaterial2Cells[i][j - k].setDensity(material2GhostDensity);
                        newMaterial2Cells[i][j - k].setXVelocity(material1GhostVelocity);
                        newMaterial2Cells[i][j - k].setPressure(material1GhostPressure);
                    }
                    else
                    {
                        EulerStateVector riemannProblemSolution = HLLCSolver::solveX(material2Cells[i][j - 1], material1Cells[i][j], material2Parameters, material1Parameters);

                        double material2GhostPressure = riemannProblemSolution.getPressure();
                        double material2GhostVelocity = riemannProblemSolution.getXVelocity();
                        double material2GhostDensity = HLLCSolver::computeLeftStarRegionDensity(material2Cells[i][j - 1], material1Cells[i][j], material2Parameters, material1Parameters);

                        newMaterial2Cells[i][j + k - 1].setDensity(material2GhostDensity);
                        newMaterial2Cells[i][j + k - 1].setXVelocity(material2GhostVelocity);
                        newMaterial2Cells[i][j + k - 1].setPressure(material2GhostPressure);

                        double material1GhostDensity = HLLCSolver::computeRightStarRegionDensity(material2Cells[i][j - 1], material1Cells[i][j], material2Parameters, material1Parameters);

                        newMaterial1Cells[i][j - k].setDensity(material1GhostDensity);
                        newMaterial1Cells[i][j - k].setXVelocity(material2GhostVelocity);
                        newMaterial1Cells[i][j - k].setPressure(material2GhostPressure);
                    }
                }

                inGhostRegion = !inGhostRegion;
            }

            levelSetValue = levelSetFunction[i][j];
        }
    }

    material1Cells = newMaterial1Cells;
    material2Cells = newMaterial2Cells;

    for (int i = 0; i < columnCount; i++)
    {
        levelSetValue = levelSetFunction[0][i];
        inGhostRegion = true;

        for (int j = 1; j < rowCount; j++)
        {
            if ((levelSetFunction[j][i] * levelSetValue) <= 0.0 && j > 2 && j < (rowCount - 3))
            {
                for (int k = 0; k < 4; k++)
                {
                    if (!inGhostRegion)
                    {
                        EulerStateVector riemannProblemSolution = HLLCSolver::solveY(material1Cells[j - 1][i], material2Cells[j][i], material1Parameters, material2Parameters);

                        double material1GhostPressure = riemannProblemSolution.getPressure();
                        double material1GhostVelocity = riemannProblemSolution.getYVelocity();
                        double material1GhostDensity = HLLCSolver::computeTopStarRegionDensity(material1Cells[j - 1][i], material2Cells[j][i], material1Parameters, material2Parameters);

                        newMaterial1Cells[j + k - 1][i].setDensity(material1GhostDensity);
                        newMaterial1Cells[j + k - 1][i].setYVelocity(material1GhostVelocity);
                        newMaterial1Cells[j + k - 1][i].setPressure(material1GhostPressure);

                        double material2GhostDensity = HLLCSolver::computeBottomStarRegionDensity(material1Cells[j - 1][i], material2Cells[j][i], material1Parameters, material2Parameters);

                        newMaterial2Cells[j - k][i].setDensity(material2GhostDensity);
                        newMaterial2Cells[j - k][i].setYVelocity(material1GhostVelocity);
                        newMaterial2Cells[j - k][i].setPressure(material1GhostPressure);
                    }
                    else
                    {
                        EulerStateVector riemannProblemSolution = HLLCSolver::solveY(material2Cells[j - 1][i], material1Cells[j][i], material2Parameters, material1Parameters);

                        double material2GhostPressure = riemannProblemSolution.getPressure();
                        double material2GhostVelocity = riemannProblemSolution.getYVelocity();
                        double material2GhostDensity = HLLCSolver::computeTopStarRegionDensity(material2Cells[j - 1][i], material1Cells[j][i], material2Parameters, material1Parameters);

                        newMaterial2Cells[j + k - 1][i].setDensity(material2GhostDensity);
                        newMaterial2Cells[j + k - 1][i].setYVelocity(material2GhostVelocity);
                        newMaterial2Cells[j + k - 1][i].setPressure(material2GhostPressure);

                        double material1GhostDensity = HLLCSolver::computeBottomStarRegionDensity(material2Cells[j - 1][i], material1Cells[j][i], material2Parameters, material1Parameters);

                        newMaterial1Cells[j - k][i].setDensity(material1GhostDensity);
                        newMaterial1Cells[j - k][i].setYVelocity(material2GhostVelocity);
                        newMaterial1Cells[j - k][i].setPressure(material2GhostPressure);
                    }
                }

                inGhostRegion = !inGhostRegion;
            }

            levelSetValue = levelSetFunction[j][i];
        }
    }

    return MultimaterialSystem(newMaterial1Cells, newMaterial2Cells, levelSetFunction);
}

vector<double> RGFMSolver::updateLevelSetFunction(vector<double> levelSetFunction, double cellSpacing, double timeStep, vector<EulerStateVector> material1Cells,
                                                  vector<EulerStateVector> material2Cells)
{
    int cellCount = levelSetFunction.size();
    vector<double> newLevelSetFunction = levelSetFunction;
    vector<double> levelSetFunctionWithBoundary = insertBoundaryCells(levelSetFunction, 2);

    double upwindApproximation;
    double velocity;

    double levelSetValue = levelSetFunction[0];
    bool inGhostRegion = false;

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

vector<vector<double> > RGFMSolver::update2DLevelSetFunctionX(vector<vector<double> > levelSetFunction, double cellSpacing, double timeStep,
                                                              vector<vector<EulerStateVector> > material1Cells, vector<vector<EulerStateVector> > material2Cells)
{
    int rowCount = levelSetFunction.size();
    int columnCount = levelSetFunction[0].size();

    vector<vector<double> > newLevelSetFunction = levelSetFunction;
    vector<vector<double> > levelSetFunctionWithBoundary = insertBoundaryCells2D(levelSetFunction, 2);

    double upwindApproximation;
    double velocity;

    double levelSetValue;
    bool inGhostRegion;

    for (int i = 0; i < rowCount; i++)
    {
        levelSetValue = levelSetFunction[i][0];
        inGhostRegion = true;

        for (int j = 0; j < columnCount; j++)
        {
            if ((levelSetValue * levelSetFunction[i][j]) <= 0.0)
            {
                inGhostRegion = !inGhostRegion;
            }

            if (!inGhostRegion)
            {
                velocity = material1Cells[i][j].getXVelocity();
            }
            else
            {
                velocity = material2Cells[i][j].getXVelocity();
            }

            if (velocity < 0)
            {
                upwindApproximation = (-levelSetFunctionWithBoundary[i + 2][j + 4] + (4 * levelSetFunctionWithBoundary[i + 2][j + 3]) - (3 * levelSetFunctionWithBoundary[i + 2][j + 2])) /
                        (2 * cellSpacing);
            }
            else
            {
                upwindApproximation = ((3 * levelSetFunctionWithBoundary[i + 2][j + 2]) - (4 * levelSetFunctionWithBoundary[i + 2][j + 1]) + levelSetFunctionWithBoundary[i + 2][j]) /
                        (2 * cellSpacing);
            }
            newLevelSetFunction[i][j] = levelSetFunction[i][j] - (timeStep * (velocity * upwindApproximation));

            levelSetValue = levelSetFunction[i][j];
        }
    }

    return newLevelSetFunction;
}

vector<vector<double> > RGFMSolver::update2DLevelSetFunctionY(vector<vector<double> > levelSetFunction, double cellSpacing, double timeStep,
                                                              vector<vector<EulerStateVector> > material1Cells, vector<vector<EulerStateVector> > material2Cells)
{
    int rowCount = levelSetFunction.size();
    int columnCount = levelSetFunction[0].size();

    vector<vector<double> > newLevelSetFunction = levelSetFunction;
    vector<vector<double> > levelSetFunctionWithBoundary = insertBoundaryCells2D(levelSetFunction, 2);

    double upwindApproximation;
    double velocity;

    double levelSetValue;
    bool inGhostRegion;

    for (int i = 0; i < columnCount; i++)
    {
        levelSetValue = levelSetFunction[0][i];
        inGhostRegion = true;

        for (int j = 0; j < rowCount; j++)
        {
            if ((levelSetValue * levelSetFunction[j][i]) <= 0.0)
            {
                inGhostRegion = !inGhostRegion;
            }

            if (!inGhostRegion)
            {
                velocity = material1Cells[j][i].getYVelocity();
            }
            else
            {
                velocity = material2Cells[j][i].getYVelocity();
            }

            if (velocity < 0)
            {
                upwindApproximation = (-levelSetFunctionWithBoundary[j + 4][i + 2] + (4 * levelSetFunctionWithBoundary[j + 3][i + 2]) - (3 * levelSetFunctionWithBoundary[j + 2][i + 2])) /
                        (2 * cellSpacing);
            }
            else
            {
                upwindApproximation = ((3 * levelSetFunctionWithBoundary[j + 2][i + 2]) - (4 * levelSetFunctionWithBoundary[j + 1][i + 2]) + levelSetFunctionWithBoundary[j][i + 2]) /
                        (2 * cellSpacing);
            }
            newLevelSetFunction[j][i] = levelSetFunction[j][i] - (timeStep * (velocity * upwindApproximation));

            levelSetValue = levelSetFunction[j][i];
        }
    }

    return newLevelSetFunction;
}

MultimaterialSystem RGFMSolver::solveExact(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                           int subcyclingIterations, int reinitialisationFrequency, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;

    vector<EulerStateVector> material1Cells = multimaterialSystem.getMaterial1Cells();
    vector<EulerStateVector> material2Cells = multimaterialSystem.getMaterial2Cells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    vector<EulerStateVector> currentMaterial1Cells = material1Cells;
    vector<EulerStateVector> currentMaterial2Cells = material2Cells;

    while (currentTime < finalTime)
    {
        MultimaterialSystem newMultimaterialSystem = applyRGFMBoundaryConditionsExact(MultimaterialSystem(currentMaterial1Cells, currentMaterial2Cells, levelSetFunction), material1Parameters,
                                                                                      material2Parameters);

        currentMaterial1Cells = newMultimaterialSystem.getMaterial1Cells();
        currentMaterial2Cells = newMultimaterialSystem.getMaterial2Cells();

        vector<EulerStateVector> currentMaterial1CellsWithBoundary = Solvers::insertBoundaryCells(currentMaterial1Cells, 2);
        vector<EulerStateVector> currentMaterial2CellsWithBoundary = Solvers::insertBoundaryCells(currentMaterial2Cells, 2);

        double timeStep = min(Solvers::computeStableTimeStep(currentMaterial1CellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters),
                              Solvers::computeStableTimeStep(currentMaterial2CellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material2Parameters));

        levelSetFunction = updateLevelSetFunction(levelSetFunction, cellSpacing, timeStep, currentMaterial1Cells, currentMaterial2Cells);

        SecondOrderSolver::computeSLICTimeStep(currentMaterial1Cells, currentMaterial1CellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters);
        SecondOrderSolver::computeSLICTimeStep(currentMaterial2Cells, currentMaterial2CellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material2Parameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return MultimaterialSystem(currentMaterial1Cells, currentMaterial2Cells, levelSetFunction);
}

MultimaterialSystem RGFMSolver::solveExact2D(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                             int subcyclingIterations, int reinitialisationFrequency, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;

    vector<vector<EulerStateVector> > material1Cells = multimaterialSystem.getMaterial1Cells2D();
    vector<vector<EulerStateVector> > material2Cells = multimaterialSystem.getMaterial2Cells2D();
    vector<vector<double> > levelSetFunction = multimaterialSystem.getLevelSetFunction2D();

    vector<vector<EulerStateVector> > currentMaterial1Cells = material1Cells;
    vector<vector<EulerStateVector> > currentMaterial2Cells = material2Cells;

    while (currentTime < finalTime)
    {
        MultimaterialSystem newMultimaterialSystem = applyRGFMBoundaryConditionsExact2D(MultimaterialSystem(currentMaterial1Cells, currentMaterial2Cells, levelSetFunction), material1Parameters,
                                                                                        material2Parameters);

        currentMaterial1Cells = newMultimaterialSystem.getMaterial1Cells2D();
        currentMaterial2Cells = newMultimaterialSystem.getMaterial2Cells2D();

        vector<vector<EulerStateVector> > currentMaterial1CellsWithBoundary = Solvers::insertBoundaryCells2D(currentMaterial1Cells, 2);
        vector<vector<EulerStateVector> > currentMaterial2CellsWithBoundary = Solvers::insertBoundaryCells2D(currentMaterial2Cells, 2);

        double timeStep = min(Solvers::computeStableTimeStep2D(currentMaterial1CellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters),
                              Solvers::computeStableTimeStep2D(currentMaterial2CellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material2Parameters));

        levelSetFunction = update2DLevelSetFunctionX(levelSetFunction, cellSpacing, 0.5 * timeStep, currentMaterial1Cells, currentMaterial2Cells);
        levelSetFunction = update2DLevelSetFunctionY(levelSetFunction, cellSpacing, timeStep, currentMaterial1Cells, currentMaterial2Cells);
        levelSetFunction = update2DLevelSetFunctionX(levelSetFunction, cellSpacing, 0.5 * timeStep, currentMaterial1Cells, currentMaterial2Cells);

        SecondOrderSolver::computeXSLICTimeStep2D(currentMaterial1Cells, currentMaterial1CellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, material1Parameters);

        currentMaterial1CellsWithBoundary = Solvers::insertBoundaryCells2D(currentMaterial1Cells, 2);
        SecondOrderSolver::computeYSLICTimeStep2D(currentMaterial1Cells, currentMaterial1CellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters);

        currentMaterial1CellsWithBoundary = Solvers::insertBoundaryCells2D(currentMaterial1Cells, 2);
        SecondOrderSolver::computeXSLICTimeStep2D(currentMaterial1Cells, currentMaterial1CellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, material1Parameters);

        SecondOrderSolver::computeXSLICTimeStep2D(currentMaterial2Cells, currentMaterial2CellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, material2Parameters);

        currentMaterial2CellsWithBoundary = Solvers::insertBoundaryCells2D(currentMaterial2Cells, 2);
        SecondOrderSolver::computeYSLICTimeStep2D(currentMaterial2Cells, currentMaterial2CellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material2Parameters);

        currentMaterial2CellsWithBoundary = Solvers::insertBoundaryCells2D(currentMaterial2Cells, 2);
        SecondOrderSolver::computeXSLICTimeStep2D(currentMaterial2Cells, currentMaterial2CellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, material2Parameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return MultimaterialSystem(currentMaterial1Cells, currentMaterial2Cells, levelSetFunction);
}

MultimaterialSystem RGFMSolver::solveHLLC(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                          int subcyclingIterations, int reinitialisationFrequency, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;

    vector<EulerStateVector> material1Cells = multimaterialSystem.getMaterial1Cells();
    vector<EulerStateVector> material2Cells = multimaterialSystem.getMaterial2Cells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    vector<EulerStateVector> currentMaterial1Cells = material1Cells;
    vector<EulerStateVector> currentMaterial2Cells = material2Cells;

    while (currentTime < finalTime)
    {
        MultimaterialSystem newMultimaterialSystem = applyRGFMBoundaryConditionsHLLC(MultimaterialSystem(currentMaterial1Cells, currentMaterial2Cells, levelSetFunction), material1Parameters,
                                                                                     material2Parameters);

        currentMaterial1Cells = newMultimaterialSystem.getMaterial1Cells();
        currentMaterial2Cells = newMultimaterialSystem.getMaterial2Cells();

        vector<EulerStateVector> currentMaterial1CellsWithBoundary = Solvers::insertBoundaryCells(currentMaterial1Cells, 2);
        vector<EulerStateVector> currentMaterial2CellsWithBoundary = Solvers::insertBoundaryCells(currentMaterial2Cells, 2);

        double timeStep = min(Solvers::computeStableTimeStep(currentMaterial1CellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters),
                              Solvers::computeStableTimeStep(currentMaterial2CellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material2Parameters));

        levelSetFunction = updateLevelSetFunction(levelSetFunction, cellSpacing, timeStep, currentMaterial1Cells, currentMaterial2Cells);

        SecondOrderSolver::computeSLICTimeStep(currentMaterial1Cells, currentMaterial1CellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters);
        SecondOrderSolver::computeSLICTimeStep(currentMaterial2Cells, currentMaterial2CellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material2Parameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return MultimaterialSystem(currentMaterial1Cells, currentMaterial2Cells, levelSetFunction);
}

MultimaterialSystem RGFMSolver::solveHLLC2D(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                            int subcyclingIterations, int reinitialisationFrequency, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;

    vector<vector<EulerStateVector> > material1Cells = multimaterialSystem.getMaterial1Cells2D();
    vector<vector<EulerStateVector> > material2Cells = multimaterialSystem.getMaterial2Cells2D();
    vector<vector<double> > levelSetFunction = multimaterialSystem.getLevelSetFunction2D();

    vector<vector<EulerStateVector> > currentMaterial1Cells = material1Cells;
    vector<vector<EulerStateVector> > currentMaterial2Cells = material2Cells;

    while (currentTime < finalTime)
    {
        MultimaterialSystem newMultimaterialSystem = applyRGFMBoundaryConditionsHLLC2D(MultimaterialSystem(currentMaterial1Cells, currentMaterial2Cells, levelSetFunction), material1Parameters,
                                                                                       material2Parameters);

        currentMaterial1Cells = newMultimaterialSystem.getMaterial1Cells2D();
        currentMaterial2Cells = newMultimaterialSystem.getMaterial2Cells2D();

        vector<vector<EulerStateVector> > currentMaterial1CellsWithBoundary = Solvers::insertBoundaryCells2D(currentMaterial1Cells, 2);
        vector<vector<EulerStateVector> > currentMaterial2CellsWithBoundary = Solvers::insertBoundaryCells2D(currentMaterial2Cells, 2);

        double timeStep = min(Solvers::computeStableTimeStep2D(currentMaterial1CellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters),
                              Solvers::computeStableTimeStep2D(currentMaterial2CellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material2Parameters));

        levelSetFunction = update2DLevelSetFunctionX(levelSetFunction, cellSpacing, 0.5 * timeStep, currentMaterial1Cells, currentMaterial2Cells);
        levelSetFunction = update2DLevelSetFunctionY(levelSetFunction, cellSpacing, timeStep, currentMaterial1Cells, currentMaterial2Cells);
        levelSetFunction = update2DLevelSetFunctionX(levelSetFunction, cellSpacing, 0.5 * timeStep, currentMaterial1Cells, currentMaterial2Cells);

        SecondOrderSolver::computeXSLICTimeStep2D(currentMaterial1Cells, currentMaterial1CellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, material1Parameters);

        currentMaterial1CellsWithBoundary = Solvers::insertBoundaryCells2D(currentMaterial1Cells, 2);
        SecondOrderSolver::computeYSLICTimeStep2D(currentMaterial1Cells, currentMaterial1CellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters);

        currentMaterial1CellsWithBoundary = Solvers::insertBoundaryCells2D(currentMaterial1Cells, 2);
        SecondOrderSolver::computeXSLICTimeStep2D(currentMaterial1Cells, currentMaterial1CellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, material1Parameters);

        SecondOrderSolver::computeXSLICTimeStep2D(currentMaterial2Cells, currentMaterial2CellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, material2Parameters);

        currentMaterial2CellsWithBoundary = Solvers::insertBoundaryCells2D(currentMaterial2Cells, 2);
        SecondOrderSolver::computeYSLICTimeStep2D(currentMaterial2Cells, currentMaterial2CellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material2Parameters);

        currentMaterial2CellsWithBoundary = Solvers::insertBoundaryCells2D(currentMaterial2Cells, 2);
        SecondOrderSolver::computeXSLICTimeStep2D(currentMaterial2Cells, currentMaterial2CellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, material2Parameters);

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return MultimaterialSystem(currentMaterial1Cells, currentMaterial2Cells, levelSetFunction);
}

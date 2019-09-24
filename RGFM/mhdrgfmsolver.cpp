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

MHDMultimaterialSystem MHDRGFMSolver::applyRGFMBoundaryConditions2D(MHDMultimaterialSystem multimaterialSystem, MHDMaterialParameters material1Parameters,
                                                                    MHDMaterialParameters material2Parameters)
{
    vector<vector<MHDStateVector> > material1Cells = multimaterialSystem.getMaterial1Cells2D();
    vector<vector<MHDStateVector> > material2Cells = multimaterialSystem.getMaterial2Cells2D();
    vector<vector<double> > levelSetFunction = multimaterialSystem.getLevelSetFunction2D();

    vector<vector<MHDStateVector> > newMaterial1Cells = material1Cells;
    vector<vector<MHDStateVector> > newMaterial2Cells = material2Cells;

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
                        MHDStateVector riemannProblemSolution = MHDHLLCSolver::solveX(material1Cells[i][j - 1], material2Cells[i][j], material1Parameters, material2Parameters);

                        double material1GhostAuxiliaryField = riemannProblemSolution.getAuxiliaryField();
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

                        double material1GhostDensity = MHDHLLCSolver::computeLeftStarRegionDensity(material1Cells[i][j - 1], material2Cells[i][j], material1Parameters, material2Parameters);

                        newMaterial1Cells[i][j + k - 1].setDensity(material1GhostDensity);

                        newMaterial1Cells[i][j + k - 1].setXVelocity(material1GhostXVelocity);
                        newMaterial1Cells[i][j + k - 1].setYVelocity(material1GhostYVelocity);
                        newMaterial1Cells[i][j + k - 1].setZVelocity(material1GhostZVelocity);

                        newMaterial1Cells[i][j + k - 1].setPressure(material1GhostPressure);

                        newMaterial1Cells[i][j + k - 1].setXMagneticField(material1GhostXMagneticField);
                        newMaterial1Cells[i][j + k - 1].setYMagneticField(material1GhostYMagneticField);
                        newMaterial1Cells[i][j + k - 1].setZMagneticField(material1GhostZMagneticField);
                        newMaterial1Cells[i][j + k - 1].setAuxiliaryField(material1GhostAuxiliaryField);

                        double material2GhostDensity = MHDHLLCSolver::computeRightStarRegionDensity(material1Cells[i][j - 1], material2Cells[i][j], material1Parameters, material2Parameters);

                        newMaterial2Cells[i][j - k].setDensity(material2GhostDensity);

                        newMaterial2Cells[i][j - k].setXVelocity(material1GhostXVelocity);
                        newMaterial2Cells[i][j - k].setYVelocity(material1GhostYVelocity);
                        newMaterial2Cells[i][j - k].setZVelocity(material1GhostZVelocity);

                        newMaterial2Cells[i][j - k].setPressure(material1GhostPressure);

                        newMaterial2Cells[i][j - k].setXMagneticField(material1GhostXMagneticField);
                        newMaterial2Cells[i][j - k].setYMagneticField(material1GhostYMagneticField);
                        newMaterial2Cells[i][j - k].setZMagneticField(material1GhostZMagneticField);
                        newMaterial2Cells[i][j - k].setAuxiliaryField(material1GhostAuxiliaryField);
                    }
                    else
                    {
                        MHDStateVector riemannProblemSolution = MHDHLLCSolver::solveX(material2Cells[i][j - 1], material1Cells[i][j], material2Parameters, material1Parameters);

                        double material2GhostAuxiliaryField = riemannProblemSolution.getAuxiliaryField();
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

                        double material2GhostDensity = MHDHLLCSolver::computeLeftStarRegionDensity(material2Cells[i][j - 1], material1Cells[i][j], material2Parameters, material1Parameters);

                        newMaterial2Cells[i][j + k - 1].setDensity(material2GhostDensity);

                        newMaterial2Cells[i][j + k - 1].setXVelocity(material2GhostXVelocity);
                        newMaterial2Cells[i][j + k - 1].setYVelocity(material2GhostYVelocity);
                        newMaterial2Cells[i][j + k - 1].setZVelocity(material2GhostZVelocity);

                        newMaterial2Cells[i][j + k - 1].setPressure(material2GhostPressure);

                        newMaterial2Cells[i][j + k - 1].setXMagneticField(material2GhostXMagneticField);
                        newMaterial2Cells[i][j + k - 1].setYMagneticField(material2GhostYMagneticField);
                        newMaterial2Cells[i][j + k - 1].setZMagneticField(material2GhostZMagneticField);
                        newMaterial2Cells[i][j + k - 1].setAuxiliaryField(material2GhostAuxiliaryField);

                        double material1GhostDensity = MHDHLLCSolver::computeRightStarRegionDensity(material2Cells[i][j - 1], material1Cells[i][j], material2Parameters, material1Parameters);

                        newMaterial1Cells[i][j - k].setDensity(material1GhostDensity);

                        newMaterial1Cells[i][j - k].setXVelocity(material2GhostXVelocity);
                        newMaterial1Cells[i][j - k].setYVelocity(material2GhostYVelocity);
                        newMaterial1Cells[i][j - k].setZVelocity(material2GhostZVelocity);

                        newMaterial1Cells[i][j - k].setPressure(material2GhostPressure);

                        newMaterial1Cells[i][j - k].setXMagneticField(material2GhostXMagneticField);
                        newMaterial1Cells[i][j - k].setYMagneticField(material2GhostYMagneticField);
                        newMaterial1Cells[i][j - k].setZMagneticField(material2GhostZMagneticField);
                        newMaterial1Cells[i][j - k].setAuxiliaryField(material2GhostAuxiliaryField);
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
            if ((levelSetFunction[j][i]  * levelSetValue) <= 0.0 && j > 2 && j < (rowCount - 3))
            {
                for (int k = 0; k < 4; k++)
                {
                    if (!inGhostRegion)
                    {
                        MHDStateVector riemannProblemSolution = MHDHLLCSolver::solveY(material1Cells[j - 1][i], material2Cells[j][i], material1Parameters, material2Parameters);

                        double material1GhostAuxiliaryField = riemannProblemSolution.getAuxiliaryField();
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

                        double material1GhostDensity = MHDHLLCSolver::computeTopStarRegionDensity(material1Cells[j - 1][i], material2Cells[j][i], material1Parameters, material2Parameters);

                        newMaterial1Cells[j + k - 1][i].setDensity(material1GhostDensity);

                        newMaterial1Cells[j + k - 1][i].setXVelocity(material1GhostXVelocity);
                        newMaterial1Cells[j + k - 1][i].setYVelocity(material1GhostYVelocity);
                        newMaterial1Cells[j + k - 1][i].setZVelocity(material1GhostZVelocity);

                        newMaterial1Cells[j + k - 1][i].setPressure(material1GhostPressure);

                        newMaterial1Cells[j + k - 1][i].setXMagneticField(material1GhostXMagneticField);
                        newMaterial1Cells[j + k - 1][i].setYMagneticField(material1GhostYMagneticField);
                        newMaterial1Cells[j + k - 1][i].setZMagneticField(material1GhostZMagneticField);
                        newMaterial1Cells[j + k - 1][i].setAuxiliaryField(material1GhostAuxiliaryField);

                        double material2GhostDensity = MHDHLLCSolver::computeBottomStarRegionDensity(material1Cells[j - 1][i], material2Cells[j][i], material1Parameters, material2Parameters);

                        newMaterial2Cells[j - k][i].setDensity(material2GhostDensity);

                        newMaterial2Cells[j - k][i].setXVelocity(material1GhostXVelocity);
                        newMaterial2Cells[j - k][i].setYVelocity(material1GhostYVelocity);
                        newMaterial2Cells[j - k][i].setZVelocity(material1GhostZVelocity);

                        newMaterial2Cells[j - k][i].setPressure(material1GhostPressure);

                        newMaterial2Cells[j - k][i].setXMagneticField(material1GhostXMagneticField);
                        newMaterial2Cells[j - k][i].setYMagneticField(material1GhostYMagneticField);
                        newMaterial2Cells[j - k][i].setZMagneticField(material1GhostZMagneticField);
                        newMaterial2Cells[j - k][i].setAuxiliaryField(material1GhostAuxiliaryField);
                    }
                    else
                    {
                        MHDStateVector riemannProblemSolution = MHDHLLCSolver::solveY(material2Cells[j - 1][i], material1Cells[j][i], material2Parameters, material1Parameters);

                        double material2GhostAuxiliaryField = riemannProblemSolution.getAuxiliaryField();
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

                        double material2GhostDensity = MHDHLLCSolver::computeTopStarRegionDensity(material2Cells[j - 1][i], material1Cells[j][i], material2Parameters, material1Parameters);

                        newMaterial2Cells[j + k - 1][i].setDensity(material2GhostDensity);

                        newMaterial2Cells[j + k - 1][i].setXVelocity(material2GhostXVelocity);
                        newMaterial2Cells[j + k - 1][i].setYVelocity(material2GhostYVelocity);
                        newMaterial2Cells[j + k - 1][i].setZVelocity(material2GhostZVelocity);

                        newMaterial2Cells[j + k - 1][i].setPressure(material2GhostPressure);

                        newMaterial2Cells[j + k - 1][i].setXMagneticField(material2GhostXMagneticField);
                        newMaterial2Cells[j + k - 1][i].setYMagneticField(material2GhostYMagneticField);
                        newMaterial2Cells[j + k - 1][i].setZMagneticField(material2GhostZMagneticField);
                        newMaterial2Cells[j + k - 1][i].setAuxiliaryField(material2GhostAuxiliaryField);

                        double material1GhostDensity = MHDHLLCSolver::computeBottomStarRegionDensity(material2Cells[j - 1][i], material1Cells[j][i], material2Parameters, material1Parameters);

                        newMaterial1Cells[j - k][i].setDensity(material1GhostDensity);

                        newMaterial1Cells[j - k][i].setXVelocity(material2GhostXVelocity);
                        newMaterial1Cells[j - k][i].setYVelocity(material2GhostYVelocity);
                        newMaterial1Cells[j - k][i].setZVelocity(material2GhostZVelocity);

                        newMaterial1Cells[j - k][i].setPressure(material2GhostPressure);

                        newMaterial1Cells[j - k][i].setXMagneticField(material2GhostXMagneticField);
                        newMaterial1Cells[j - k][i].setYMagneticField(material2GhostYMagneticField);
                        newMaterial1Cells[j - k][i].setZMagneticField(material2GhostZMagneticField);
                        newMaterial1Cells[j - k][i].setAuxiliaryField(material2GhostAuxiliaryField);
                    }
                }

                inGhostRegion = !inGhostRegion;
            }

            levelSetValue = levelSetFunction[j][i];
        }
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

vector<vector<double> > MHDRGFMSolver::update2DLevelSetFunctionX(vector<vector<double> > levelSetFunction, double cellSpacing, double timeStep,
                                                                 vector<vector<MHDStateVector> > material1Cells, vector<vector<MHDStateVector> > material2Cells)
{
    int rowCount = levelSetFunction.size();
    int columnCount = levelSetFunction[0].size();

    vector<vector<double> > newLevelSetFunction = levelSetFunction;
    vector<vector<double> > levelSetFunctionWithBoundary = RGFMSolver::insertBoundaryCells2D(levelSetFunction, 2);

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
            if ((levelSetValue * levelSetFunction[i][j] <= 0.0))
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

            if (velocity < 0.0)
            {
                upwindApproximation = (-levelSetFunctionWithBoundary[i + 2][j + 4] + (4.0 * levelSetFunctionWithBoundary[i + 2][j + 3]) - (3.0 * levelSetFunctionWithBoundary[i + 2][j + 2])) /
                        (2.0 * cellSpacing);
            }
            else
            {
                upwindApproximation = ((3.0 * levelSetFunctionWithBoundary[i + 2][j + 2]) - (4.0 * levelSetFunctionWithBoundary[i + 2][j + 1]) + levelSetFunctionWithBoundary[i + 2][j]) /
                        (2.0 * cellSpacing);
            }
            newLevelSetFunction[i][j] = levelSetFunction[i][j] - (timeStep * (velocity * upwindApproximation));

            levelSetValue = levelSetFunction[i][j];
        }
    }

    return newLevelSetFunction;
}

vector<vector<double> > MHDRGFMSolver::update2DLevelSetFunctionY(vector<vector<double> > levelSetFunction, double cellSpacing, double timeStep,
                                                                 vector<vector<MHDStateVector> > material1Cells, vector<vector<MHDStateVector> > material2Cells)
{
    int rowCount = levelSetFunction.size();
    int columnCount = levelSetFunction[0].size();

    vector<vector<double> > newLevelSetFunction = levelSetFunction;
    vector<vector<double> > levelSetFunctionWithBoundary = RGFMSolver::insertBoundaryCells2D(levelSetFunction, 2);

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

            if (velocity < 0.0)

            {
                upwindApproximation = (-levelSetFunctionWithBoundary[j + 4][i + 2] + (4.0 * levelSetFunctionWithBoundary[j + 3][i + 2]) - (3.0 * levelSetFunctionWithBoundary[j + 2][i + 2])) /
                        (2.0 * cellSpacing);
            }
            else
            {
                upwindApproximation = ((3.0 * levelSetFunctionWithBoundary[j + 2][i + 2]) - (4.0 * levelSetFunctionWithBoundary[j + 1][i + 2]) + levelSetFunctionWithBoundary[j][i + 2]) /
                        (2.0 * cellSpacing);
            }
            newLevelSetFunction[j][i] = levelSetFunction[j][i] - (timeStep * (velocity * upwindApproximation));

            levelSetValue = levelSetFunction[j][i];
        }
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

MHDMultimaterialSystem MHDRGFMSolver::solve2D(MHDMultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int slopeLimiter,
                                              int subcyclingIterations, int reinitialisationFrequency, MHDMaterialParameters material1Parameters, MHDMaterialParameters material2Parameters)
{
    double currentTime = 0.0;
    int currentIteration = 0;

    vector<vector<MHDStateVector> > material1Cells = multimaterialSystem.getMaterial1Cells2D();
    vector<vector<MHDStateVector> > material2Cells = multimaterialSystem.getMaterial2Cells2D();
    vector<vector<double> > levelSetFunction = multimaterialSystem.getLevelSetFunction2D();

    vector<vector<MHDStateVector> > currentMaterial1Cells = material1Cells;
    vector<vector<MHDStateVector> > currentMaterial2Cells = material2Cells;

    while (currentTime < finalTime)
    {
        MHDMultimaterialSystem newMultimaterialSystem = applyRGFMBoundaryConditions2D(MHDMultimaterialSystem(currentMaterial1Cells, currentMaterial2Cells, levelSetFunction), material1Parameters,
                                                                                      material2Parameters);

        currentMaterial1Cells = newMultimaterialSystem.getMaterial1Cells2D();
        currentMaterial2Cells = newMultimaterialSystem.getMaterial2Cells2D();

        vector<vector<MHDStateVector> > currentMaterial1CellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentMaterial1Cells, 2);
        vector<vector<MHDStateVector> > currentMaterial2CellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentMaterial2Cells, 2);

        double timeStep = min(MHDSolvers::computeStableTimeStep2D(currentMaterial1CellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material1Parameters),
                              MHDSolvers::computeStableTimeStep2D(currentMaterial2CellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration, material2Parameters));

        levelSetFunction = update2DLevelSetFunctionX(levelSetFunction, cellSpacing, 0.5 * timeStep, currentMaterial1Cells, currentMaterial2Cells);
        levelSetFunction = update2DLevelSetFunctionY(levelSetFunction, cellSpacing, timeStep, currentMaterial1Cells, currentMaterial2Cells);
        levelSetFunction = update2DLevelSetFunctionX(levelSetFunction, cellSpacing, 0.5 * timeStep, currentMaterial1Cells, currentMaterial2Cells);

        double maximumWaveSpeed = min(MHDSolvers::computeMaximumWaveSpeed2D(currentMaterial1CellsWithBoundary, material1Parameters),
                                      MHDSolvers::computeMaximumWaveSpeed2D(currentMaterial2CellsWithBoundary, material2Parameters));
        material1Parameters.setHyperbolicWaveSpeed(maximumWaveSpeed);
        material2Parameters.setHyperbolicWaveSpeed(maximumWaveSpeed);

        material1Parameters.configureParabolicDamping();
        material2Parameters.configureParabolicDamping();
        /*
        double material1MaximumWaveSpeed = MHDSolvers::computeMaximumWaveSpeed2D(currentMaterial1CellsWithBoundary, material1Parameters);
        double material2MaximumWaveSpeed = MHDSolvers::computeMaximumWaveSpeed2D(currentMaterial2CellsWithBoundary, material2Parameters);

        material1Parameters.setHyperbolicWaveSpeed(material1MaximumWaveSpeed);
        material2Parameters.setHyperbolicWaveSpeed(material2MaximumWaveSpeed);

        material1Parameters.configureParabolicDamping();
        material2Parameters.configureParabolicDamping();
        */

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentMaterial1CellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentMaterial1Cells, 1);

            MHDForcingSolver::computeRungeKuttaTimeStep2D(currentMaterial1Cells, currentMaterial1CellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter,
                                                          material1Parameters);
        }

        currentMaterial1CellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentMaterial1Cells, 2);
        MHDSecondOrderSolver::computeXSLICTimeStep2D(currentMaterial1Cells, currentMaterial1CellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, material1Parameters);

        currentMaterial1CellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentMaterial1Cells, 2);
        MHDSecondOrderSolver::computeYSLICTimeStep2D(currentMaterial1Cells, currentMaterial1CellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material1Parameters);

        currentMaterial1CellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentMaterial1Cells, 2);
        MHDSecondOrderSolver::computeXSLICTimeStep2D(currentMaterial1Cells, currentMaterial1CellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, material1Parameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentMaterial1CellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentMaterial1Cells, 1);

            MHDForcingSolver::computeRungeKuttaTimeStep2D(currentMaterial1Cells, currentMaterial1CellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter,
                                                          material1Parameters);
        }

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentMaterial2CellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentMaterial2Cells, 1);

            MHDForcingSolver::computeRungeKuttaTimeStep2D(currentMaterial2Cells, currentMaterial2CellsWithBoundary, cellSpacing, 0.5 * (timeStep / subcyclingIterations), bias, slopeLimiter,
                                                          material2Parameters);
        }

        currentMaterial2CellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentMaterial2Cells, 2);
        MHDSecondOrderSolver::computeXSLICTimeStep2D(currentMaterial2Cells, currentMaterial2CellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, material2Parameters);

        currentMaterial2CellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentMaterial2Cells, 2);
        MHDSecondOrderSolver::computeYSLICTimeStep2D(currentMaterial2Cells, currentMaterial2CellsWithBoundary, cellSpacing, timeStep, bias, slopeLimiter, material2Parameters);

        currentMaterial2CellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentMaterial2Cells, 2);
        MHDSecondOrderSolver::computeXSLICTimeStep2D(currentMaterial2Cells, currentMaterial2CellsWithBoundary, cellSpacing, 0.5 * timeStep, bias, slopeLimiter, material2Parameters);

        for (int i = 0; i < subcyclingIterations; i++)
        {
            currentMaterial2CellsWithBoundary = MHDSolvers::insertBoundaryCells2D(currentMaterial2Cells, 1);

            MHDForcingSolver::computeRungeKuttaTimeStep2D(currentMaterial2Cells, currentMaterial2CellsWithBoundary, cellSpacing, 0.5 * (timeStep /subcyclingIterations), bias, slopeLimiter,
                                                          material2Parameters);
        }

        currentTime += timeStep;
        currentIteration += 1;

        Solvers::outputStatus(currentIteration, currentTime, timeStep);
    }

    return MHDMultimaterialSystem(currentMaterial1Cells, currentMaterial2Cells, levelSetFunction);
}

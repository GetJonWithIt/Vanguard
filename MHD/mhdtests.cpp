#include "mhdtests.h"

MHDTests::MHDTests()
{
}

void MHDTests::solveDumbserTest1(int cellCount, int subcyclingIterations)
{
    double cellSpacing = 1.0 / cellCount;

    vector<MHDStateVector> initialCells(cellCount);
    MHDMaterialParameters materialParameters(2.0);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = MHDStateVector(1.0, 0.0, 0.0, 0.0, 1.0, 0.75, 1.0, 0.0, 0.0);
        }
        else
        {
            initialCells[i] = MHDStateVector(0.125, 0.0, 0.0, 0.0, 0.1, 0.75, -1.0, 0.0, 0.0);
        }
    }

    outputSolution(MHDSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.1, 0.0, 0, subcyclingIterations, materialParameters));
}

void MHDTests::solveDumbserTest2(int cellCount, int subcyclingIterations)
{
    double cellSpacing = 1.0 / cellCount;

    vector<MHDStateVector> initialCells(cellCount);
    MHDMaterialParameters materialParameters(5.0 / 3.0);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.4 * cellCount))
        {
            initialCells[i] = MHDStateVector(1.08, 1.2, 0.01, 0.5, 0.95, 0.564189, 1.015541, 0.564189, 0.0);
        }
        else
        {
            initialCells[i] = MHDStateVector(0.9891, -0.0131, 0.0269, 0.010037, 0.97159, 0.564189, 1.135262, 0.564923, 0.0);
        }
    }

    outputSolution(MHDSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.2, 0.0, 0, subcyclingIterations, materialParameters));
}

void MHDTests::solveDumbserTest3(int cellCount, int subcyclingIterations)
{
    double cellSpacing = 1.0 / cellCount;

    vector<MHDStateVector> initialCells(cellCount);
    MHDMaterialParameters materialParameters(5.0 / 3.0);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = MHDStateVector(1.0, 0.0, 0.0, 0.0, 1.0, 1.3, 1.0, 0.0, 0.0);
        }
        else
        {
            initialCells[i] = MHDStateVector(0.4, 0.0, 0.0, 0.0, 0.4, 1.3, -1.0, 0.0, 0.0);
        }
    }

    outputSolution(MHDSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.16, 0.0, 0, subcyclingIterations, materialParameters));
}

void MHDTests::solve2DDumbserTest1(int cellCount, int subcyclingIterations)
{
    double cellSpacing = 1.0 / cellCount;

    vector<vector<MHDStateVector> > initialCells(cellCount, vector<MHDStateVector>(cellCount));
    MHDMaterialParameters materialParameters(2.0);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.2 * cellCount))
            {
                initialCells[i][j] = MHDStateVector(1.0, 0.0, 0.0, 0.0, 1.0, 0.75, 1.0, 0.0, 0.0);
            }
            else
            {
                initialCells[i][j] = MHDStateVector(0.125, 0.0, 0.0, 0.0, 0.1, 0.75, -1.0, 0.0, 0.0);
            }
        }
    }

    outputSolution2D(MHDSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.6, 0.05, 0.0, 0, subcyclingIterations, materialParameters));
}

void MHDTests::solve2DDumbserTest2(int cellCount, int subcyclingIterations)
{
    double cellSpacing = 1.0 / cellCount;

    vector<vector<MHDStateVector> > initialCells(cellCount, vector<MHDStateVector>(cellCount));
    MHDMaterialParameters materialParameters(5.0 / 3.0);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.16 * cellCount))
            {
                initialCells[i][j] = MHDStateVector(1.08, 1.2, 0.01, 0.5, 0.95, 0.564189, 1.015541, 0.564189, 0.0);
            }
            else
            {
                initialCells[i][j] = MHDStateVector(0.9891, -0.0131, 0.0269, 0.010037, 0.97159, 0.564189, 1.135262, 0.564923, 0.0);
            }
        }
    }

    outputSolution2D(MHDSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.8, 0.1, 0.0, 0, subcyclingIterations, materialParameters));
}

void MHDTests::outputSolution(vector<MHDStateVector> solution)
{
    int cellCount = solution.size();
    double cellSpacing = 1.0 / cellCount;

    ofstream densityFile("density.dat");
    ofstream yMagneticFieldFile("yMagneticField.dat");

    for (int i = 0; i < cellCount; i++)
    {
        densityFile << (cellSpacing * i) << " " << solution[i].getDensity() << endl;
        yMagneticFieldFile << (cellSpacing * i) << " " << solution[i].getYMagneticField() << endl;
    }

    densityFile.close();
    yMagneticFieldFile.close();
}

void MHDTests::outputSolution2D(vector<vector<MHDStateVector> > solution)
{
    int rowCount = solution.size();
    int columnCount = solution[0].size();
    double cellSpacing = 1.0 / rowCount;

    ofstream densityFile("density.dat");
    ofstream yMagneticFieldFile("yMagneticField.dat");

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            densityFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << solution[i][j].getDensity() << endl;
            yMagneticFieldFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << solution[i][j].getYMagneticField() << endl;
        }
    }

    densityFile.close();
    yMagneticFieldFile.close();
}

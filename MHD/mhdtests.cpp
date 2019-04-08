#include "mhdtests.h"

MHDTests::MHDTests()
{
}

void MHDTests::solveDumbserTest1(int cellCount)
{
    double cellSpacing = 1.0 / cellCount;

    vector<MHDStateVector> initialCells(cellCount);
    MHDMaterialParameters materialParameters(5.0 / 3.0);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = MHDStateVector(1.0, 0.0, 0.0, 0.0, 1.0, 3.0 / 4.0, 1.0, 0.0, 0.0);
        }
        else
        {
            initialCells[i] = MHDStateVector(0.125, 0.0, 0.0, 0.0, 0.1, 3.0 / 4.0, -1.0, 0.0, 0.0);
        }
    }

    outputSolution(MHDSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.1, 0.0, 0, 0, materialParameters));
}

void MHDTests::solveDumbserTest2(int cellCount)
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

    outputSolution(MHDSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.2, 0.0, 0, 0, materialParameters));
}

void MHDTests::solveDumbserTest3(int cellCount)
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

    outputSolution(MHDSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.16, 0.0, 0, 0, materialParameters));
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

#include "eulertests.h"

EulerTests::EulerTests()
{
}

void EulerTests::solveToroTest1(int cellCount)
{
    double cellSpacing = 1.0 / cellCount;

    vector<EulerStateVector> initialCells(cellCount);
    EulerMaterialParameters materialParameters(1.4);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = EulerStateVector(1.0, 0.0, 0.0, 0.0, 1.0);
        }
        else
        {
            initialCells[i] = EulerStateVector(0.125, 0.0, 0.0, 0.0, 0.1);
        }
    }

    outputSolution(SecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.25, 0.0, 0, 0, materialParameters));
}

void EulerTests::solve2DToroTest1(int cellCount)
{
    double cellSpacing = 1.0 / cellCount;

    vector<vector<EulerStateVector> > initialCells(cellCount, vector<EulerStateVector>(cellCount));
    EulerMaterialParameters materialParameters(1.4);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.2 * cellCount))
            {
                initialCells[i][j] = EulerStateVector(1.0, 0.0, 0.0, 0.0, 1.0);
            }
            else
            {
                initialCells[i][j] = EulerStateVector(0.125, 0.0, 0.0, 0.0, 0.1);
            }
        }
    }

    outputSolution2D(SecondOrderSolver::solve2D(initialCells, cellSpacing, 0.8, 0.125, 0.0, 0, 0, materialParameters));
}

void EulerTests::outputSolution(vector<EulerStateVector> solution)
{
    int cellCount = solution.size();
    double cellSpacing = 1.0 / cellCount;

    ofstream densityFile("density.dat");
    ofstream xVelocityFile("xVelocity.dat");
    ofstream yVelocityFile("yVelocity.dat");
    ofstream zVelocityFile("zVelocity.dat");
    ofstream pressureFile("pressure.dat");

    for (int i = 0; i < cellCount; i++)
    {
        densityFile << (cellSpacing * i) << " " << solution[i].getDensity() << endl;
        xVelocityFile << (cellSpacing * i) << " " << solution[i].getXVelocity() << endl;
        yVelocityFile << (cellSpacing * i) << " " << solution[i].getYVelocity() << endl;
        zVelocityFile << (cellSpacing * i) << " " << solution[i].getZVelocity() << endl;
        pressureFile << (cellSpacing * i) << " " << solution[i].getPressure() << endl;
    }

    densityFile.close();
    xVelocityFile.close();
    yVelocityFile.close();
    zVelocityFile.close();
    pressureFile.close();
}

void EulerTests::outputSolution2D(vector<vector<EulerStateVector> > solution)
{
    int rowCount = solution.size();
    int columnCount = solution[0].size();
    double cellSpacing = 1.0 / rowCount;

    ofstream densityFile("density.dat");
    ofstream xVelocityFile("xVelocity.dat");
    ofstream yVelocityFile("yVelocity.dat");
    ofstream zVelocityFile("zVelocity.dat");
    ofstream pressureFile("pressure.dat");

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            densityFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << solution[i][j].getDensity() << endl;
            pressureFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << solution[i][j].getPressure() << endl;
        }
    }

    densityFile.close();
    xVelocityFile.close();
    yVelocityFile.close();
    zVelocityFile.close();
    pressureFile.close();
}

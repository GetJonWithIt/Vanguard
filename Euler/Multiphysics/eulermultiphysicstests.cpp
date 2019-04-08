#include "eulermultiphysicstests.h"

EulerMultiphysicsTests::EulerMultiphysicsTests()
{
}

void EulerMultiphysicsTests::solveToroTest1(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<EulerMultiphysicsStateVector> initialCells(cellCount);
    EulerMaterialParameters materialParameters(1.4);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = EulerMultiphysicsStateVector(0.999, 0.0, 0.0, 0.0, 1.0, 1.0, 0.125, 0.1);
        }
        else
        {
            initialCells[i] = EulerMultiphysicsStateVector(0.001, 0.0, 0.0, 0.0, 1.0, 1.0, 0.125, 0.1);
        }
    }

    outputSolution(SecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.25, 0.0, 0, 0, reinitialisationFrequency, materialParameters, materialParameters));
}

void EulerMultiphysicsTests::solveFedkiwTest(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<EulerMultiphysicsStateVector> initialCells(cellCount);
    EulerMaterialParameters material1Parameters(1.4);
    EulerMaterialParameters material2Parameters(1.67);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            if (i <= (0.05 * cellCount))
            {
                initialCells[i] = EulerMultiphysicsStateVector(0.999, 0.3535 * sqrt(pow(10.0, 5.0)), 0.0, 0.0, 1.3333, 1.5 * pow(10.0, 5.0), 0.1379, 1.0 * pow(10.0, 5.0));
            }
            else
            {
                initialCells[i] = EulerMultiphysicsStateVector(0.999, 0.0, 0.0, 0.0, 1.0, 1.0 * pow(10.0, 5.0), 0.1379, 1.0 * pow(10.0, 5.0));
            }
        }
        else
        {
            initialCells[i] = EulerMultiphysicsStateVector(0.001, 0.0, 0.0, 0.0, 1.0, 1.0 * pow(10.0, 5.0), 0.1379, 1.0 * pow(10.0, 5.0));
        }
    }

    outputSolution(SecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.0012, 0.0, 0, 0, reinitialisationFrequency, material1Parameters, material2Parameters));
}

void EulerMultiphysicsTests::solve2DToroTest1(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<vector<EulerMultiphysicsStateVector> > initialCells(cellCount, vector<EulerMultiphysicsStateVector>(cellCount));
    EulerMaterialParameters materialParameters(1.4);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.2 * cellCount))
            {
                initialCells[i][j] = EulerMultiphysicsStateVector(0.999, 0.0, 0.0, 0.0, 1.0, 1.0, 0.125, 0.1);
            }
            else
            {
                initialCells[i][j] = EulerMultiphysicsStateVector(0.001, 0.0, 0.0, 0.0, 1.0, 1.0, 0.125, 0.1);
            }
        }
    }

    outputSolution2D(SecondOrderSolver::solve2D(initialCells, cellSpacing, 0.8, 0.125, 0.0, 0, 0, reinitialisationFrequency, materialParameters, materialParameters));
}

void EulerMultiphysicsTests::solve2DFedkiwTest(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<vector<EulerMultiphysicsStateVector> > initialCells(cellCount, vector<EulerMultiphysicsStateVector>(cellCount));
    EulerMaterialParameters material1Parameters(1.4);
    EulerMaterialParameters material2Parameters(1.67);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.2 * cellCount))
            {
                if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.02 * cellCount))
                {
                    initialCells[i][j] = EulerMultiphysicsStateVector(0.999, 0.3535 * sqrt(pow(10.0, 5.0)), 0.0, 0.0, 1.3333, 1.5 * pow(10.0, 5.0), 0.1379, 1.0 * pow(10.0, 5.0));
                }
                else
                {
                    initialCells[i][j] = EulerMultiphysicsStateVector(0.999, 0.0, 0.0, 0.0, 1.0, 1.0 * pow(10.0, 5.0), 0.1379, 1.0 * pow(10.0, 5.0));
                }
            }
            else
            {
                initialCells[i][j] = EulerMultiphysicsStateVector(0.001, 0.0, 0.0, 0.0, 1.0, 1.0 * pow(10.0, 5.0), 0.1379, 1.0 * pow(10.0, 5.0));
            }
        }
    }

    outputSolution2D(SecondOrderSolver::solve2D(initialCells, cellSpacing, 0.8, 0.0006, 0.0, 0, 0, reinitialisationFrequency, material1Parameters, material2Parameters));
}

void EulerMultiphysicsTests::outputSolution(vector<EulerMultiphysicsStateVector> solution)
{
    int cellCount = solution.size();
    double cellSpacing = 1.0 / cellCount;

    ofstream volumeFractionFile("volumeFraction.dat");
    ofstream densityFile("density.dat");

    for (int i = 0; i < cellCount; i++)
    {
        volumeFractionFile << (cellSpacing * i) << " " << solution[i].getMaterial1VolumeFraction() << endl;
        densityFile << (cellSpacing * i) << " " << solution[i].computeTotalDensity() << endl;
    }

    volumeFractionFile.close();
    densityFile.close();
}

void EulerMultiphysicsTests::outputSolution2D(vector<vector<EulerMultiphysicsStateVector> > solution)
{
    int rowCount = solution.size();
    int columnCount = solution[0].size();
    double cellSpacing = 1.0 / rowCount;

    ofstream volumeFractionFile("volumeFraction.dat");
    ofstream densityFile("density.dat");

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            volumeFractionFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << solution[i][j].getMaterial1VolumeFraction() << endl;
            densityFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << solution[i][j].computeTotalDensity() << endl;
        }
    }
}

#include "mhdreducedtests.h"

MHDReducedTests::MHDReducedTests()
{
}

void MHDReducedTests::solveDumbserTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<MHDReducedStateVector> initialCells(cellCount);
    MHDMaterialParameters materialParameters(2.0);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = MHDReducedStateVector(0.999, 0.0, 0.0, 0.0, 1.0, 0.125, 1.0, 0.75, 1.0, 0.0, 0.0);
        }
        else
        {
            initialCells[i] = MHDReducedStateVector(0.001, 0.0, 0.0, 0.0, 1.0, 0.125, 0.1, 0.75, -1.0, 0.0, 0.0);
        }
    }

    outputSolution(MHDSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.1, 0.0, 0, subcyclingIterations, reinitialisationFrequency, materialParameters, materialParameters));
}

void MHDReducedTests::solveDumbserMultimaterialTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<MHDReducedStateVector> initialCells(cellCount);
    MHDMaterialParameters material1Parameters(1.4);
    MHDMaterialParameters material2Parameters(1.67);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = MHDReducedStateVector(0.999, 0.0, 0.0, 0.0, 1.0, 0.125, 1.0, 0.75, 1.0, 0.0, 0.0);
        }
        else
        {
            initialCells[i] = MHDReducedStateVector(0.001, 0.0, 0.0, 0.0, 1.0, 0.125, 0.1, 0.75, -1.0, 0.0, 0.0);
        }
    }

    outputSolution(MHDSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.1, 0.0, 0, subcyclingIterations, reinitialisationFrequency, material1Parameters, material2Parameters));
}

void MHDReducedTests::solve2DDumbserTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<vector<MHDReducedStateVector> > initialCells(cellCount, vector<MHDReducedStateVector>(cellCount));
    MHDMaterialParameters materialParameters(2.0);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.2 * cellCount))
            {
                initialCells[i][j] = MHDReducedStateVector(0.999, 0.0, 0.0, 0.0, 1.0, 0.125, 1.0, 0.75, 1.0, 0.0, 0.0);
            }
            else
            {
                initialCells[i][j] = MHDReducedStateVector(0.001, 0.0, 0.0, 0.0, 1.0, 0.125, 0.1, 0.75, -1.0, 0.0, 0.0);
            }
        }
    }

    outputSolution2D(MHDSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.8, 0.05, 0.0, 0, subcyclingIterations, reinitialisationFrequency, materialParameters, materialParameters));
}

void MHDReducedTests::solve2DDumbserMultimaterialTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<vector<MHDReducedStateVector> > initialCells(cellCount, vector<MHDReducedStateVector>(cellCount));
    MHDMaterialParameters material1Parameters(1.4);
    MHDMaterialParameters material2Parameters(1.67);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.2 * cellCount))
            {
                initialCells[i][j] = MHDReducedStateVector(0.999, 0.0, 0.0, 0.0, 1.0, 0.125, 1.0, 0.75, 1.0, 0.0, 0.0);
            }
            else
            {
                initialCells[i][j] = MHDReducedStateVector(0.001, 0.0, 0.0, 0.0, 1.0, 0.125, 0.1, 0.75, -1.0, 0.0, 0.0);
            }
        }
    }

    outputSolution2D(MHDSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.8, 0.05, 0.0, 0, subcyclingIterations, reinitialisationFrequency, material1Parameters, material2Parameters));
}

void MHDReducedTests::outputSolution(vector<MHDReducedStateVector> solution)
{
    int cellCount = solution.size();
    double cellSpacing = 1.0 / cellCount;

    ofstream densityFile("density.dat");
    ofstream volumeFractionFile("volumeFraction.dat");

    for (int i = 0; i < cellCount; i++)
    {
        densityFile << (cellSpacing * i) << " " << solution[i].computeTotalDensity() << endl;
        volumeFractionFile << (cellSpacing * i) << " " << solution[i].getMaterial1VolumeFraction() << endl;
    }

    densityFile.close();
    volumeFractionFile.close();
}

void MHDReducedTests::outputSolution2D(vector<vector<MHDReducedStateVector> > solution)
{
    int rowCount = solution.size();
    int columnCount = solution[0].size();
    double cellSpacing = 1.0 / rowCount;

    ofstream densityFile("density.dat");
    ofstream volumeFractionFile("volumeFraction.dat");

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            densityFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << solution[i][j].computeTotalDensity() << endl;
            volumeFractionFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << solution[i][j].getMaterial1VolumeFraction() << endl;
        }
    }

    densityFile.close();
    volumeFractionFile.close();
}

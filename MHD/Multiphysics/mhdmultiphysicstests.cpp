#include "mhdmultiphysicstests.h"

MHDMultiphysicsTests::MHDMultiphysicsTests()
{
}

void MHDMultiphysicsTests::solveDumbserTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<MHDMultiphysicsStateVector> initialCells(cellCount);
    MHDMaterialParameters materialParameters(2.0);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = MHDMultiphysicsStateVector(0.999, 0.0, 0.0, 0.0, 1.0, 1.0, 0.75, 1.0, 0.0, 0.0, 0.125, 0.1, 0.75, -1.0, 0.0, 0.0);
        }
        else
        {
            initialCells[i] = MHDMultiphysicsStateVector(0.001, 0.0, 0.0, 0.0, 1.0, 1.0, 0.75, 1.0, 0.0, 0.0, 0.125, 0.1, 0.75, -1.0, 0.0, 0.0);
        }
    }

    outputSolution(MHDSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.1, 0.0, 0, subcyclingIterations, reinitialisationFrequency, materialParameters, materialParameters));
}

void MHDMultiphysicsTests::solveDumbserMultimaterialTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<MHDMultiphysicsStateVector> initialCells(cellCount);
    MHDMaterialParameters material1Parameters(1.4);
    MHDMaterialParameters material2Parameters(1.67);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = MHDMultiphysicsStateVector(0.999, 0.0, 0.0, 0.0, 1.0, 1.0, 0.75, 1.0, 0.0, 0.0, 0.125, 0.1, 0.75, -1.0, 0.0, 0.0);
        }
        else
        {
            initialCells[i] = MHDMultiphysicsStateVector(0.001, 0.0, 0.0, 0.0, 1.0, 1.0, 0.75, 1.0, 0.0, 0.0, 0.125, 0.1, 0.75, -1.0, 0.0, 0.0);
        }
    }

    outputSolution(MHDSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.2, 0.0, 0, subcyclingIterations, reinitialisationFrequency, material1Parameters, material2Parameters));
}

void MHDMultiphysicsTests::outputSolution(vector<MHDMultiphysicsStateVector> solution)
{
    int cellCount = solution.size();
    double cellSpacing = 1.0 / cellCount;

    ofstream densityFile("density.dat");
    ofstream volumeFractionFile("volumeFraction.dat");

    for (int i = 0; i < cellCount; i++)
    {
        densityFile << (cellSpacing * i) << " " << solution[i].computeTotalDensity() << endl;
        volumeFractionFile << (cellSpacing * i) <<  " " << solution[i].getMaterial1VolumeFraction() << endl;
    }

    densityFile.close();
    volumeFractionFile.close();
}

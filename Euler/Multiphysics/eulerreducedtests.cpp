#include "eulerreducedtests.h"

EulerReducedTests::EulerReducedTests()
{

}

void EulerReducedTests::solveToroTest1(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<EulerReducedStateVector> initialCells(cellCount);
    EulerMaterialParameters materialParameters = EulerMaterialParameters(1.4);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = EulerReducedStateVector(0.999, 0.0, 0.0, 0.0, 1.0, 0.125, 1.0);
        }
        else
        {
            initialCells[i] = EulerReducedStateVector(0.001, 0.0, 0.0, 0.0, 1.0, 0.125, 0.1);
        }
    }

    outputSolution(SecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.25, 0.0, 0, 0, reinitialisationFrequency, materialParameters, materialParameters));
}

void EulerReducedTests::solveFedkiwTest(int cellCount, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<EulerReducedStateVector> initialCells(cellCount);
    EulerMaterialParameters material1Parameters(1.4);
    EulerMaterialParameters material2Parameters(1.67);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            if (i <= (0.05 * cellCount))
            {
                initialCells[i] = EulerReducedStateVector(0.999, 0.3535 * sqrt(pow(10.0, 5.0)), 0.0, 0.0, 1.3333, 0.1379, 1.5 * pow(10.0, 5.0));
            }
            else
            {
                initialCells[i] = EulerReducedStateVector(0.999, 0.0, 0.0, 0.0, 1.0, 0.1379, 1.0 * pow(10.0, 5.0));
            }
        }
        else
        {
            initialCells[i] = EulerReducedStateVector(0.001, 0.0, 0.0, 0.0, 1.0, 0.1379, 1.0 * pow(10.0, 5.0));
        }
    }

    outputSolution(SecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.0012, 0.0, 0, 0, reinitialisationFrequency, material1Parameters, material2Parameters));
}

void EulerReducedTests::outputSolution(vector<EulerReducedStateVector> solution)
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

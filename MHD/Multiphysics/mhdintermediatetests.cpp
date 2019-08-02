#include "mhdintermediatetests.h"

MHDIntermediateTests::MHDIntermediateTests()
{
}

void MHDIntermediateTests::solveDumbserTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<MHDIntermediateStateVector> initialCells(cellCount);
    MHDMaterialParameters materialParameters(2.0);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = MHDIntermediateStateVector(0.999, 0.0, 0.0, 0.0, 1.0, 1.0, 0.125, 0.1, 0.75, 1.0, 0.0, 0.0);
        }
        else
        {
            initialCells[i] = MHDIntermediateStateVector(0.001, 0.0, 0.0, 0.0, 1.0, 1.0, 0.125, 0.1, 0.75, -1.0, 0.0, 0.0);
        }
    }

    outputSolution(MHDSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.1, 0.0, 0, subcyclingIterations, reinitialisationFrequency, materialParameters, materialParameters));
}

void MHDIntermediateTests::solveDumbserMultimaterialTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<MHDIntermediateStateVector> initialCells(cellCount);
    MHDMaterialParameters material1Parameters(1.4);
    MHDMaterialParameters material2Parameters(1.67);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (0.5 * cellCount))
        {
            initialCells[i] = MHDIntermediateStateVector(0.999, 0.0, 0.0, 0.0, 1.0, 1.0, 0.125, 0.1, 0.75, 1.0, 0.0, 0.0);
        }
        else
        {
            initialCells[i] = MHDIntermediateStateVector(0.001, 0.0, 0.0, 0.0, 1.0, 1.0, 0.125, 0.1, 0.75, -1.0, 0.0, 0.0);
        }
    }

    outputSolution(MHDSecondOrderSolver::solve(initialCells, cellSpacing, 0.8, 0.2, 0.0, 0, subcyclingIterations, reinitialisationFrequency, material1Parameters, material2Parameters));
}

void MHDIntermediateTests::solve2DDumbserTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<vector<MHDIntermediateStateVector> > initialCells(cellCount, vector<MHDIntermediateStateVector>(cellCount));
    MHDMaterialParameters materialParameters(2.0);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.2 * cellCount))
            {
                initialCells[i][j] = MHDIntermediateStateVector(0.999, 0.0, 0.0, 0.0, 1.0, 1.0, 0.125, 0.1, 0.75, 1.0, 0.0, 0.0);
            }
            else
            {
                initialCells[i][j] = MHDIntermediateStateVector(0.001, 0.0, 0.0, 0.0, 1.0, 1.0, 0.125, 0.1, 0.75, -1.0, 0.0, 0.0);
            }
        }
    }

    outputSolution2D(MHDSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.8, 0.05, 0.0, 0, subcyclingIterations, reinitialisationFrequency, materialParameters, materialParameters));
}

void MHDIntermediateTests::solve2DDumbserMultimaterialTest1(int cellCount, int subcyclingIterations, int reinitialisationFrequency)
{
    double cellSpacing = 1.0 / cellCount;

    vector<vector<MHDIntermediateStateVector> > initialCells(cellCount, vector<MHDIntermediateStateVector>(cellCount));
    MHDMaterialParameters material1Parameters(1.4);
    MHDMaterialParameters material2Parameters(1.67);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.2 * cellCount))
            {
                initialCells[i][j] = MHDIntermediateStateVector(0.999, 0.0, 0.0, 0.0, 1.0, 1.0, 0.125, 0.1, 0.75, 1.0, 0.0, 0.0);
            }
            else
            {
                initialCells[i][j] = MHDIntermediateStateVector(0.001, 0.0, 0.0, 0.0, 1.0, 1.0, 0.125, 0.1, 0.75, -1.0, 0.0, 0.0);
            }
        }
    }

    outputSolution2D(MHDSecondOrderSolver::solve2D(initialCells, cellSpacing, 0.8, 0.1, 0.0, 0, subcyclingIterations, reinitialisationFrequency, material1Parameters, material2Parameters));
}

void MHDIntermediateTests::outputSolution(vector<MHDIntermediateStateVector> solution)
{
    int cellCount = solution.size();
    double cellSpacing = 1.0 / cellCount;

    ofstream densityFile("density.dat");
    ofstream yMagneticFieldFile("yMagneticField.dat");
    ofstream volumeFractionFile("volumeFraction.dat");

    for (int i = 0; i < cellCount; i++)
    {
        densityFile << (cellSpacing * i) << " " << solution[i].computeTotalDensity() << endl;
        yMagneticFieldFile << (cellSpacing * i) << " " << solution[i].getInterfaceYMagneticField() << endl;
        volumeFractionFile << (cellSpacing * i) << " " << solution[i].getMaterial1VolumeFraction() << endl;
    }

    densityFile.close();
    yMagneticFieldFile.close();
    volumeFractionFile.close();
}

void MHDIntermediateTests::outputSolution2D(vector<vector<MHDIntermediateStateVector> > solution)
{
    int rowCount = solution.size();
    int columnCount = solution[0].size();
    double cellSpacing = 1.0 / rowCount;

    ofstream densityFile("density.dat");
    ofstream yMagneticFieldFile("yMagneticField.dat");
    ofstream volumeFractionFile("volumeFraction.dat");

    for (int i = 0; i < rowCount; i++)
    {
        for (int j = 0; j < columnCount; j++)
        {
            densityFile << (cellSpacing * i ) << " " << (cellSpacing * j) << " " << solution[i][j].computeTotalDensity() << endl;
            yMagneticFieldFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << solution[i][j].getInterfaceYMagneticField() << endl;
            volumeFractionFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << solution[i][j].getMaterial1VolumeFraction() << endl;
        }
    }

    densityFile.close();
    yMagneticFieldFile.close();
    volumeFractionFile.close();
}

#include "eulerrgfmtests.h"

EulerRGFMTests::EulerRGFMTests()
{
}

void EulerRGFMTests::solveToroTest1Exact(int cellCount, int reinitialisationFrequency)
{
    vector<EulerStateVector> material1Cells(cellCount);
    vector<EulerStateVector> material2Cells(cellCount);
    vector<double> levelSetFunction(cellCount);

    double cellSpacing = 1.0 / cellCount;
    EulerMaterialParameters materialParameters(1.4);

    for (int i = 0; i < cellCount; i++)
    {
        material1Cells[i] = EulerStateVector(1.0, 0.0, 0.0, 0.0, 1.0);
        material2Cells[i] = EulerStateVector(0.125, 0.0, 0.0, 0.0, 0.1);

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputSolution(RGFMSolver::solveExact(MultimaterialSystem(material1Cells, material2Cells, levelSetFunction), cellSpacing, 0.8, 0.25, 0.0, 0, 0, 0, materialParameters, materialParameters));
}

void EulerRGFMTests::solveToroTest1HLLC(int cellCount, int reinitialisationFrequency)
{
    vector<EulerStateVector> material1Cells(cellCount);
    vector<EulerStateVector> material2Cells(cellCount);
    vector<double> levelSetFunction(cellCount);

    double cellSpacing = 1.0 / cellCount;
    EulerMaterialParameters materialParameters(1.4);

    for (int i = 0; i < cellCount; i++)
    {
        material1Cells[i] = EulerStateVector(1.0, 0.0, 0.0, 0.0, 1.0);
        material2Cells[i] = EulerStateVector(0.125, 0.0, 0.0, 0.0, 0.1);

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputSolution(RGFMSolver::solveHLLC(MultimaterialSystem(material1Cells, material2Cells, levelSetFunction), cellSpacing, 0.8, 0.25, 0.0, 0, 0, 0, materialParameters, materialParameters));
}

void EulerRGFMTests::solveFedkiwTestExact(int cellCount, int reinitialisationFrequency)
{
    vector<EulerStateVector> material1Cells(cellCount);
    vector<EulerStateVector> material2Cells(cellCount);
    vector<double> levelSetFunction(cellCount);

    double cellSpacing = 1.0 / cellCount;
    EulerMaterialParameters material1Parameters(1.4);
    EulerMaterialParameters material2Parameters(1.67);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (cellCount * 0.05))
        {
            material1Cells[i] = EulerStateVector(1.3333, 0.3535 * sqrt(pow(10.0, 5.0)), 0.0, 0.0, 1.5 * pow(10.0, 5.0));
        }
        else
        {
            material1Cells[i] = EulerStateVector(1.0, 0.0, 0.0, 0.0, 1.0 * pow(10.0, 5.0));
        }

        material2Cells[i] = EulerStateVector(0.1379, 0.0, 0.0, 0.0, 1.0 * pow(10.0, 5.0));

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputSolution(RGFMSolver::solveExact(MultimaterialSystem(material1Cells, material2Cells, levelSetFunction), cellSpacing, 0.8, 0.0012, 0.0, 0, 0, 0, material1Parameters, material2Parameters));
}

void EulerRGFMTests::solveFedkiwTestHLLC(int cellCount, int reinitialisationFrequency)
{
    vector<EulerStateVector> material1Cells(cellCount);
    vector<EulerStateVector> material2Cells(cellCount);
    vector<double> levelSetFunction(cellCount);

    double cellSpacing = 1.0 / cellCount;
    EulerMaterialParameters material1Parameters(1.4);
    EulerMaterialParameters material2Parameters(1.67);

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (cellCount * 0.05))
        {
            material1Cells[i] = EulerStateVector(1.3333, 0.3535 * sqrt(pow(10.0, 5.0)), 0.0, 0.0, 1.5 * pow(10.0, 5.0));
        }
        else
        {
            material1Cells[i] = EulerStateVector(1.0, 0.0, 0.0, 0.0, 1.0 * pow(10.0, 5.0));
        }

        material2Cells[i] = EulerStateVector(0.1379, 0.0, 0.0, 0.0, 1.0 * pow(10.0, 5.0));

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputSolution(RGFMSolver::solveHLLC(MultimaterialSystem(material1Cells, material2Cells, levelSetFunction), cellSpacing, 0.8, 0.0012, 0.0, 0, 0, 0, material1Parameters, material2Parameters));
}

void EulerRGFMTests::solve2DToroTest1Exact(int cellCount, int reinitialisationFrequency)
{
    vector<vector<EulerStateVector> > material1Cells(cellCount, vector<EulerStateVector>(cellCount));
    vector<vector<EulerStateVector> > material2Cells(cellCount, vector<EulerStateVector>(cellCount));
    vector<vector<double> > levelSetFunction(cellCount, vector<double>(cellCount));

    double cellSpacing = 1.0 / cellCount;
    EulerMaterialParameters materialParameters(1.4);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            material1Cells[i][j] = EulerStateVector(1.0, 0.0, 0.0, 0.0, 1.0);
            material2Cells[i][j] = EulerStateVector(0.125, 0.0, 0.0, 0.0, 0.1);

            levelSetFunction[i][j] = (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) / cellCount) - 0.2;
        }
    }

    outputSolution2D(RGFMSolver::solveExact2D(MultimaterialSystem(material1Cells, material2Cells, levelSetFunction), cellSpacing, 0.8, 0.125, 0.0, 0, 0, 0, materialParameters, materialParameters));
}

void EulerRGFMTests::solve2DToroTest1HLLC(int cellCount, int reinitialisationFrequency)
{
    vector<vector<EulerStateVector> > material1Cells(cellCount, vector<EulerStateVector>(cellCount));
    vector<vector<EulerStateVector> > material2Cells(cellCount, vector<EulerStateVector>(cellCount));
    vector<vector<double> > levelSetFunction(cellCount, vector<double>(cellCount));

    double cellSpacing = 1.0 / cellCount;
    EulerMaterialParameters materialParameters(1.4);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            material1Cells[i][j] = EulerStateVector(1.0, 0.0, 0.0, 0.0, 1.0);
            material2Cells[i][j] = EulerStateVector(0.125, 0.0, 0.0, 0.0, 0.1);

            levelSetFunction[i][j] = (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) / cellCount) - 0.2;
        }
    }

    outputSolution2D(RGFMSolver::solveHLLC2D(MultimaterialSystem(material1Cells, material2Cells, levelSetFunction), cellSpacing, 0.8, 0.125, 0.0, 0, 0, 0, materialParameters, materialParameters));
}

void EulerRGFMTests::solve2DFedkiwTestExact(int cellCount, int reinitialisationFrequency)
{
    vector<vector<EulerStateVector> > material1Cells(cellCount, vector<EulerStateVector>(cellCount));
    vector<vector<EulerStateVector> > material2Cells(cellCount, vector<EulerStateVector>(cellCount));
    vector<vector<double> > levelSetFunction(cellCount, vector<double>(cellCount));

    double cellSpacing = 1.0 / cellCount;
    EulerMaterialParameters material1Parameters(1.4);
    EulerMaterialParameters material2Parameters(1.67);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.02 * cellCount))
            {
                material1Cells[i][j] = EulerStateVector(1.3333, 0.3535 * sqrt(pow(10.0, 5.0)), 0.0, 0.0, 1.5 * pow(10.0, 5.0));
            }
            else
            {
                material1Cells[i][j] = EulerStateVector(1.0, 0.0, 0.0, 0.0, 1.0 * pow(10.0, 5.0));
            }

            material2Cells[i][j] = EulerStateVector(0.1379, 0.0, 0.0, 0.0, 1.0 * pow(10.0, 5.0));

            levelSetFunction[i][j] = (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) / cellCount) - 0.2;
        }
    }

    outputSolution2D(RGFMSolver::solveExact2D(MultimaterialSystem(material1Cells, material2Cells, levelSetFunction), cellSpacing, 0.8, 0.0006, 0.0, 0, 0, 0, material1Parameters, material2Parameters));
}

void EulerRGFMTests::solve2DFedkiwTestHLLC(int cellCount, int reinitialisationFrequency)
{
    vector<vector<EulerStateVector> > material1Cells(cellCount, vector<EulerStateVector>(cellCount));
    vector<vector<EulerStateVector> > material2Cells(cellCount, vector<EulerStateVector>(cellCount));
    vector<vector<double> > levelSetFunction(cellCount, vector<double>(cellCount));

    double cellSpacing = 1.0 / cellCount;
    EulerMaterialParameters material1Parameters(1.4);
    EulerMaterialParameters material2Parameters(1.67);

    for (int i = 0; i < cellCount; i++)
    {
        for (int j = 0; j < cellCount; j++)
        {
            if (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) <= (0.02 * cellCount))
            {
                material1Cells[i][j] = EulerStateVector(1.3333, 0.3535 * sqrt(pow(10.0, 5.0)), 0.0, 0.0, 1.5 * pow(10.0, 5.0));
            }
            else
            {
                material1Cells[i][j] = EulerStateVector(1.0, 0.0, 0.0, 0.0, 1.0 * pow(10.0, 5.0));
            }

            material2Cells[i][j] = EulerStateVector(0.1379, 0.0, 0.0, 0.0, 1.0 * pow(10.0, 5.0));

            levelSetFunction[i][j] = (sqrt(((i - (0.5 * cellCount)) * (i - (0.5 * cellCount))) + ((j - (0.5 * cellCount)) * (j - (0.5 * cellCount)))) / cellCount) - 0.2;
        }
    }

    outputSolution2D(RGFMSolver::solveHLLC2D(MultimaterialSystem(material1Cells, material2Cells, levelSetFunction), cellSpacing, 0.8, 0.0006, 0.0, 0, 0, 0, material1Parameters, material2Parameters));
}

void EulerRGFMTests::outputSolution(MultimaterialSystem multimaterialSystem)
{
    vector<EulerStateVector> material1Cells = multimaterialSystem.getMaterial1Cells();
    vector<EulerStateVector> material2Cells = multimaterialSystem.getMaterial2Cells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    int cellCount = levelSetFunction.size();
    double cellSpacing = 1.0 / cellCount;

    ofstream levelSetFunctionFile("levelSetFunction.dat");
    ofstream densityFile("density.dat");
    ofstream pressureFile("pressure.dat");

    bool inGhostRegion = false;
    double levelSetValue = levelSetFunction[0];

    for (int i = 0; i < cellCount; i++)
    {
        if ((levelSetFunction[i] * levelSetValue) <= 0.0 && (levelSetFunction[i - 1] * levelSetValue) > 0.0)
        {
            inGhostRegion = !inGhostRegion;
        }

        if (!inGhostRegion)
        {
            densityFile << (cellSpacing * i) << " " << material1Cells[i].getDensity() << endl;
            pressureFile << (cellSpacing * i) << " " << material1Cells[i].getPressure() << endl;
        }
        else
        {
            densityFile << (cellSpacing * i) << " " << material2Cells[i].getDensity() << endl;
            pressureFile << (cellSpacing * i) << " " << material2Cells[i].getPressure() << endl;
        }

        levelSetFunctionFile << (cellSpacing * i) << " " << levelSetFunction[i] << endl;

        levelSetValue = levelSetFunction[i];
    }

    levelSetFunctionFile.close();
    densityFile.close();
    pressureFile.close();
}

void EulerRGFMTests::outputSolution2D(MultimaterialSystem multimaterialSystem)
{
    vector<vector<EulerStateVector> > material1Cells = multimaterialSystem.getMaterial1Cells2D();
    vector<vector<EulerStateVector> > material2Cells = multimaterialSystem.getMaterial2Cells2D();
    vector<vector<double> > levelSetFunction = multimaterialSystem.getLevelSetFunction2D();

    int rowCount = levelSetFunction.size();
    int columnCount = levelSetFunction[0].size();
    double cellSpacing = 1.0 / rowCount;

    ofstream levelSetFunctionFile("levelSetFunction.dat");
    ofstream densityFile("density.dat");
    ofstream pressureFile("pressure.dat");

    ofstream levelSetFunctionSliceFile("levelSetFunctionSlice.dat");
    ofstream densitySliceFile("densitySlice.dat");
    ofstream pressureSliceFile("pressureSlice.dat");

    bool inGhostRegion;
    double levelSetValue;

    for (int i = 0; i < rowCount; i++)
    {
        inGhostRegion = true;
        levelSetValue = levelSetFunction[i][0];

        for (int j = 0; j < columnCount; j++)
        {
            if ((levelSetFunction[i][j] * levelSetValue) <= 0.0 && (levelSetFunction[i][j - 1] * levelSetValue) > 0.0)
            {
                inGhostRegion = !inGhostRegion;
            }

            if (!inGhostRegion)
            {
                densityFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << material1Cells[i][j].getDensity() << endl;
                pressureFile << (cellSpacing * i) << " " << (cellSpacing * j) <<  " " << material1Cells[i][j].getPressure() << endl;

                if (i == floor(0.5 * rowCount))
                {
                    densitySliceFile << (cellSpacing * j) << " " << material1Cells[i][j].getDensity() << endl;
                    pressureSliceFile << (cellSpacing * j) << " " << material1Cells[i][j].getPressure() << endl;
                }
            }
            else
            {
                densityFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << material2Cells[i][j].getDensity() << endl;
                pressureFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << material2Cells[i][j].getPressure() << endl;

                if (i == floor(0.5 * rowCount))
                {
                    densitySliceFile << (cellSpacing * j) << " " << material2Cells[i][j].getDensity() << endl;
                    pressureSliceFile << (cellSpacing * j) << " " << material2Cells[i][j].getPressure() << endl;
                }
            }

            levelSetFunctionFile << (cellSpacing * i) << " " << (cellSpacing * j) << " " << levelSetFunction[i][j] << endl;

            if (i == floor(0.5 * rowCount))
            {
                levelSetFunctionSliceFile << (cellSpacing * j) << levelSetFunction[i][j] << endl;
            }

            levelSetValue = levelSetFunction[i][j];
        }
    }

    levelSetFunctionFile.close();
    densityFile.close();
    pressureFile.close();

    levelSetFunctionSliceFile.close();
    densitySliceFile.close();
    pressureSliceFile.close();

}

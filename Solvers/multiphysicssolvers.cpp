#include "multiphysicssolvers.h"

MultiphysicsSolvers::MultiphysicsSolvers()
{
}

void MultiphysicsSolvers::reinitialiseVolumeFraction(vector<EulerMultiphysicsStateVector> & currentCells, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    int interfaceLocation = 0;
    int cellCount = currentCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        if (currentCells[i].getMaterial1VolumeFraction() < 0.5 && currentCells[i - 1].getMaterial1VolumeFraction() >= 0.5)
        {
            interfaceLocation = i;
        }
    }

    for (int i = interfaceLocation - 2; i < interfaceLocation + 2; i++)
    {
        currentCells[i].relaxTotalDensity();
        currentCells[i].relaxTotalPressure(material1Parameters, material2Parameters);

        if (i < interfaceLocation)
        {
            currentCells[i].setMaterial1VolumeFraction(0.999);
        }
        else
        {
            currentCells[i].setMaterial1VolumeFraction(0.001);
        }
    }
}

void MultiphysicsSolvers::reinitialiseVolumeFraction(vector<vector<EulerMultiphysicsStateVector> > & currentCells, EulerMaterialParameters material1Parameters,
                                                     EulerMaterialParameters material2Parameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

    for (int i = 0; i < rowCount; i++)
    {
        int interfaceCount1 = 0;
        int interfaceCount2 = 0;
        vector<int> interfaceLocations1(100);
        vector<int> interfaceLocations2(100);

        for (int j = 1; j < columnCount; j++)
        {
            if (currentCells[i][j].getMaterial1VolumeFraction() < 0.5 && currentCells[i][j - 1].getMaterial1VolumeFraction() >= 0.5)
            {
                interfaceLocations1[interfaceCount1] = j;
                interfaceCount1 += 1;
            }

            if (currentCells[i][j].getMaterial1VolumeFraction() > 0.5 && currentCells[i][j - 1].getMaterial1VolumeFraction() <= 0.5)
            {
                interfaceLocations2[interfaceCount2] = j;
                interfaceCount2 += 1;
            }
        }

        for (int k = 0; k < interfaceCount1; k++)
        {
            for (int j = interfaceLocations1[k] - 2; j < interfaceLocations1[k] + 2; j++)
            {
                currentCells[i][j].relaxTotalDensity();
                currentCells[i][j].relaxTotalPressure(material1Parameters, material2Parameters);

                if (j < interfaceLocations1[k])
                {
                    currentCells[i][j].setMaterial1VolumeFraction(0.999);
                }
                else
                {
                    currentCells[i][j].setMaterial1VolumeFraction(0.001);
                }
            }
        }

        for (int k = 0; k < interfaceCount2; k++)
        {
            for (int j = interfaceLocations2[k] - 2; j < interfaceLocations2[k] + 2; j++)
            {
                currentCells[i][j].relaxTotalDensity();
                currentCells[i][j].relaxTotalPressure(material1Parameters, material2Parameters);

                if (j < interfaceLocations2[k])
                {
                    currentCells[i][j].setMaterial1VolumeFraction(0.001);
                }
                else
                {
                    currentCells[i][j].setMaterial1VolumeFraction(0.999);
                }
            }
        }
    }

    for (int i = 0; i < columnCount; i++)
    {
        int interfaceCount1 = 0;
        int interfaceCount2 = 0;
        vector<int> interfaceLocations1(100);
        vector<int> interfaceLocations2(100);

        for (int j = 1; j < rowCount; j++)
        {
            if (currentCells[j][i].getMaterial1VolumeFraction() < 0.5 && currentCells[j - 1][i].getMaterial1VolumeFraction() >= 0.5)
            {
                interfaceLocations1[interfaceCount1] = j;
                interfaceCount1 += 1;
            }

            if (currentCells[j][i].getMaterial1VolumeFraction() > 0.5 && currentCells[j - 1][i].getMaterial1VolumeFraction() <= 0.5)
            {
                interfaceLocations2[interfaceCount2] = j;
                interfaceCount2 += 1;
            }
        }

        for (int k = 0; k < interfaceCount1; k++)
        {
            for (int j = interfaceLocations1[k] - 2; j < interfaceLocations1[k] + 2; j++)
            {
                currentCells[j][i].relaxTotalDensity();
                currentCells[j][i].relaxTotalPressure(material1Parameters, material2Parameters);

                if (j < interfaceLocations1[k])
                {
                    currentCells[j][i].setMaterial1VolumeFraction(0.999);
                }
                else
                {
                    currentCells[j][i].setMaterial1VolumeFraction(0.001);
                }
            }
        }

        for (int k = 0; k < interfaceCount2; k++)
        {
            for (int j = interfaceLocations2[k] - 2; j < interfaceLocations2[k] + 2; j++)
            {
                currentCells[j][i].relaxTotalDensity();
                currentCells[j][i].relaxTotalPressure(material1Parameters, material2Parameters);

                if (j < interfaceLocations2[k])
                {
                    currentCells[j][i].setMaterial1VolumeFraction(0.001);
                }
                else
                {
                    currentCells[j][i].setMaterial1VolumeFraction(0.999);
                }
            }
        }
    }
}

void MultiphysicsSolvers::reinitialiseVolumeFraction(vector<EulerReducedStateVector> & currentCells, EulerMaterialParameters material1Parameters, EulerMaterialParameters material2Parameters)
{
    int interfaceLocation = 0;
    int cellCount = currentCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        if (currentCells[i].getMaterial1VolumeFraction() < 0.5 && currentCells[i - 1].getMaterial1VolumeFraction() >= 0.5)
        {
            interfaceLocation = i;
        }
    }

    for (int i = interfaceLocation - 2; i < interfaceLocation + 2; i++)
    {
        currentCells[i].relaxTotalDensity();

        if (i < interfaceLocation)
        {
            currentCells[i].setMaterial1VolumeFraction(0.999);
        }
        else
        {
            currentCells[i].setMaterial1VolumeFraction(0.001);
        }
    }
}

void MultiphysicsSolvers::reinitialiseVolumeFraction(vector<vector<EulerReducedStateVector> > & currentCells, EulerMaterialParameters material1Parameters,
                                                     EulerMaterialParameters material2Parameters)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

    for (int i = 0; i < rowCount; i++)
    {
        int interfaceCount1 = 0;
        int interfaceCount2 = 0;
        vector<int> interfaceLocations1(100);
        vector<int> interfaceLocations2(100);

        for (int j = 1; j < columnCount; j++)
        {
            if (currentCells[i][j].getMaterial1VolumeFraction() < 0.5 && currentCells[i][j - 1].getMaterial1VolumeFraction() >= 0.5)
            {
                interfaceLocations1[interfaceCount1] = j;
                interfaceCount1 += 1;
            }

            if (currentCells[i][j].getMaterial1VolumeFraction() > 0.5 && currentCells[i][j - 1].getMaterial1VolumeFraction() <= 0.5)
            {
                interfaceLocations2[interfaceCount2] = j;
                interfaceCount2 += 1;
            }
        }

        for (int k = 0; k < interfaceCount1; k++)
        {
            for (int j = interfaceLocations1[k] - 2; j < interfaceLocations1[k] + 2; j++)
            {
                currentCells[i][j].relaxTotalDensity();

                if (j < interfaceLocations1[k])
                {
                    currentCells[i][j].setMaterial1VolumeFraction(0.999);
                }
                else
                {
                    currentCells[i][j].setMaterial1VolumeFraction(0.001);
                }
            }
        }

        for (int k = 0; k < interfaceCount2; k++)
        {
            for (int j = interfaceLocations2[k] - 2; j < interfaceLocations2[k] + 2; j++)
            {
                currentCells[i][j].relaxTotalDensity();

                if (j < interfaceLocations2[k])
                {
                    currentCells[i][j].setMaterial1VolumeFraction(0.001);
                }
                else
                {
                    currentCells[i][j].setMaterial1VolumeFraction(0.999);
                }
            }
        }
    }

    for (int i = 0; i < columnCount; i++)
    {
        int interfaceCount1 = 0;
        int interfaceCount2 = 0;
        vector<int> interfaceLocations1(100);
        vector<int> interfaceLocations2(100);

        for (int j = 1; j < rowCount; j++)
        {
            if (currentCells[j][i].getMaterial1VolumeFraction() < 0.5 && currentCells[j - 1][i].getMaterial1VolumeFraction() >= 0.5)
            {
                interfaceLocations1[interfaceCount1] = j;
                interfaceCount1 += 1;
            }

            if (currentCells[j][i].getMaterial1VolumeFraction() > 0.5 && currentCells[j - 1][i].getMaterial1VolumeFraction() <= 0.5)
            {
                interfaceLocations2[interfaceCount2] = j;
                interfaceCount2 += 1;
            }
        }

        for (int k = 0; k < interfaceCount1; k++)
        {
            for (int j = interfaceLocations1[k] - 2; j < interfaceLocations1[k] + 2; j++)
            {
                currentCells[j][i].relaxTotalDensity();

                if (j < interfaceLocations1[k])
                {
                    currentCells[j][i].setMaterial1VolumeFraction(0.999);
                }
                else
                {
                    currentCells[j][i].setMaterial1VolumeFraction(0.001);
                }
            }
        }

        for (int k = 0; k < interfaceCount2; k++)
        {
            for (int j = interfaceLocations2[k] - 2; j < interfaceLocations2[k] + 2; j++)
            {
                currentCells[j][i].relaxTotalDensity();

                if (j < interfaceLocations2[k])
                {
                    currentCells[j][i].setMaterial1VolumeFraction(0.001);
                }
                else
                {
                    currentCells[j][i].setMaterial1VolumeFraction(0.999);
                }
            }
        }
    }
}

void MultiphysicsSolvers::reinitialiseVolumeFraction(vector<ElasticMultiphysicsStateVector> & currentCells, HyperelasticMaterialParameters material1Parameters,
                                                     HyperelasticMaterialParameters material2Parameters)
{
    int interfaceLocation = 0;
    int cellCount = currentCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        if (currentCells[i].getMaterial1VolumeFraction() < 0.5 && currentCells[i - 1].getMaterial1VolumeFraction() >= 0.5)
        {
            interfaceLocation = i;
        }
    }

    for (int i = interfaceLocation - 2; i < interfaceLocation + 2; i++)
    {
        currentCells[i].relaxTotalDensity();
        currentCells[i].relaxTotalDistortionTensor();
        currentCells[i].relaxTotalEntropy(material1Parameters, material2Parameters);

        if (i < interfaceLocation)
        {
            currentCells[i].setMaterial1VolumeFraction(0.999);
        }
        else
        {
            currentCells[i].setMaterial1VolumeFraction(0.001);
        }
    }
}

void MultiphysicsSolvers::reinitialiseVolumeFraction(vector<ElasticReducedStateVector> & currentCells)
{
    int interfaceLocation = 0;
    int cellCount = currentCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        if (currentCells[i].getMaterial1VolumeFraction() < 0.5 && currentCells[i - 1].getMaterial1VolumeFraction() >= 0.5)
        {
            interfaceLocation = i;
        }
    }

    for (int i = interfaceLocation - 2; i < interfaceLocation + 2; i++)
    {
        currentCells[i].relaxTotalDensity();

        if (i < interfaceLocation)
        {
            currentCells[i].setMaterial1VolumeFraction(0.999);
        }
        else
        {
            currentCells[i].setMaterial1VolumeFraction(0.001);
        }
    }
}

void MultiphysicsSolvers::reinitialiseVolumeFraction(vector<vector<ElasticReducedStateVector> > & currentCells)
{
    int rowCount = currentCells.size();
    int columnCount = currentCells[0].size();

    for (int i = 0; i < rowCount; i++)
    {
        int interfaceCount1 = 0;
        int interfaceCount2 = 0;
        vector<int> interfaceLocations1(100);
        vector<int> interfaceLocations2(100);

        for (int j = 1; j < columnCount; j++)
        {
            if (currentCells[i][j].getMaterial1VolumeFraction() < 0.5 && currentCells[i][j - 1].getMaterial1VolumeFraction() >= 0.5)
            {
                interfaceLocations1[interfaceCount1] = j;
                interfaceCount1 += 1;
            }

            if (currentCells[i][j].getMaterial1VolumeFraction() > 0.5 && currentCells[i][j - 1].getMaterial1VolumeFraction() <= 0.5)
            {
                interfaceLocations2[interfaceCount2] = j;
                interfaceCount2 += 1;
            }
        }

        for (int k = 0; k < interfaceCount1; k++)
        {
            for (int j = interfaceLocations1[k] - 2; j < interfaceLocations1[k] + 2; j++)
            {
                currentCells[i][j].relaxTotalDensity();

                if (j < interfaceLocations1[k])
                {
                    currentCells[i][j].setMaterial1VolumeFraction(0.999);
                }
                else
                {
                    currentCells[i][j].setMaterial1VolumeFraction(0.001);
                }
            }
        }

        for (int k = 0; k < interfaceCount2; k++)
        {
            for (int j = interfaceLocations2[k] - 2; j < interfaceLocations2[k] + 2; j++)
            {
                currentCells[i][j].relaxTotalDensity();

                if (j < interfaceLocations2[k])
                {
                    currentCells[i][j].setMaterial1VolumeFraction(0.001);
                }
                else
                {
                    currentCells[i][j].setMaterial1VolumeFraction(0.999);
                }
            }
        }
    }

    for (int i = 0; i < columnCount; i++)
    {
        int interfaceCount1 = 0;
        int interfaceCount2 = 0;
        vector<int> interfaceLocations1(100);
        vector<int> interfaceLocations2(100);

        for (int j = 1; j < rowCount; j++)
        {
            if (currentCells[j][i].getMaterial1VolumeFraction() < 0.5 && currentCells[j - 1][i].getMaterial1VolumeFraction() >= 0.5)
            {
                interfaceLocations1[interfaceCount1] = j;
                interfaceCount1 += 1;
            }

            if (currentCells[j][i].getMaterial1VolumeFraction() > 0.5 && currentCells[j - 1][i].getMaterial1VolumeFraction() <= 0.5)
            {
                interfaceLocations2[interfaceCount2] = j;
                interfaceCount2 += 1;
            }
        }

        for (int k = 0; k < interfaceCount1; k++)
        {
            for (int j = interfaceLocations1[k] - 2; j < interfaceLocations1[k] + 2; j++)
            {
                currentCells[j][i].relaxTotalDensity();

                if (j < interfaceLocations1[k])
                {
                    currentCells[j][i].setMaterial1VolumeFraction(0.999);
                }
                else
                {
                    currentCells[j][i].setMaterial1VolumeFraction(0.001);
                }
            }
        }

        for (int k = 0; k < interfaceCount2; k++)
        {
            for (int j = interfaceLocations2[k] - 2; j < interfaceLocations2[k] + 2; j++)
            {
                currentCells[j][i].relaxTotalDensity();

                if (j < interfaceLocations2[k])
                {
                    currentCells[j][i].setMaterial1VolumeFraction(0.001);
                }
                else
                {
                    currentCells[j][i].setMaterial1VolumeFraction(0.999);
                }
            }
        }
    }
}

void MultiphysicsSolvers::reinitialiseVolumeFraction(vector<HPRReducedStateVector> & currentCells, HPRMaterialParameters materialParameters, HPRMaterialParameters material2Parameters)
{
    int interfaceLocation = 0;
    int cellCount = currentCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        if (currentCells[i].getMaterial1VolumeFraction() < 0.5 && currentCells[i - 1].getMaterial1VolumeFraction() >= 0.5)
        {
            interfaceLocation = i;
        }
    }

    for (int i = interfaceLocation - 2; i < interfaceLocation + 2; i++)
    {
        currentCells[i].relaxTotalDensity();

        if (i < interfaceLocation)
        {
            currentCells[i].setMaterial1VolumeFraction(0.999);
        }
        else
        {
            currentCells[i].setMaterial1VolumeFraction(0.001);
        }
    }
}

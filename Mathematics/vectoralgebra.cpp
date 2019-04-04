#include "vectoralgebra.h"

VectorAlgebra::VectorAlgebra()
{
}

vector<double> VectorAlgebra::addVectors(vector<double> vector1, vector<double> vector2)
{
    int componentCount = vector1.size();
    vector<double> sumVector(componentCount);

    for (int i = 0; i < componentCount; i++)
    {
        sumVector[i] = vector1[i] + vector2[i];
    }

    return sumVector;
}

vector<double> VectorAlgebra::subtractVectors(vector<double> vector1, vector<double> vector2)
{
    return addVectors(vector1, multiplyVector(-1.0, vector2));
}

vector<double> VectorAlgebra::multiplyVector(double scalar1, vector<double> vector1)
{
    int componentCount = vector1.size();
    vector<double> productVector(componentCount);

    for (int i = 0; i < componentCount; i++)
    {
        productVector[i] = scalar1 * vector1[i];
    }

    return productVector;
}

double VectorAlgebra::computeDotProduct(vector<double> vector1, vector<double> vector2)
{
    double dotProduct = 0.0;
    int componentCount = vector1.size();

    for (int i = 0; i < componentCount; i++)
    {
        dotProduct += vector1[i] * vector2[i];
    }

    return dotProduct;
}

double VectorAlgebra::computeNorm(vector<double> vector1)
{
    return sqrt(computeDotProduct(vector1, vector1));
}

double VectorAlgebra::computeComponentSum(vector<double> vector1)
{
    double componentSum = 0.0;
    int componentCount = vector1.size();

    for (int i = 0; i < componentCount; i++)
    {
        componentSum += vector1[i];
    }

    return componentSum;
}

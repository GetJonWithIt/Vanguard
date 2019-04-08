#include <QCoreApplication>

#include "Euler/eulertests.h"
#include "Euler/Multiphysics/eulermultiphysicstests.h"
#include "Euler/Multiphysics/eulerreducedtests.h"
#include "Elasticity/elastictests.h"
#include "Elasticity/Multiphysics/elasticmultiphysicstests.h"
#include "Elasticity/Multiphysics/elasticreducedtests.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    //EulerTests::solveToroTest1(400);
    //EulerTests::solve2DToroTest1(100);
    //EulerMultiphysicsTests::solveToroTest1(200, 3);
    //EulerMultiphysicsTests::solveFedkiwTest(400, 3);
    //EulerMultiphysicsTests::solveChinnayyaTest(100);

    //EulerMultiphysicsTests::solve2DToroTest1(100, 3);
    //EulerMultiphysicsTests::solve2DFedkiwTest(100, 3);

    //EulerReducedTests::solveToroTest1(100, 0);
    //EulerReducedTests::solveFedkiwTest(100, 0);

    ElasticTests::solveBartonTest1(100);
    //ElasticTests::solveBartonTest2(200);
    //ElasticTests::solve2DBartonTest1(100);
    //ElasticTests::solve2DBartonTest2(20);

    //ElasticMultiphysicsTests::solveBartonTest1(100, 3);
    //ElasticMultiphysicsTests::solveBartonTest2(100, 3);
    //ElasticMultiphysicsTests::solveZhangTest(100, 0);

    //ElasticReducedTests::solveBartonTest1(100, 3);
    //ElasticReducedTests::solveZhangTest(20, 0);
    //ElasticReducedTests::solve2DBartonTest1(100, 3);
    //ElasticReducedTests::solve2DZhangTest(100, 3);

    return a.exec();
}

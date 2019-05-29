#include <QCoreApplication>

#include "Euler/eulertests.h"
#include "Euler/Multiphysics/eulermultiphysicstests.h"
#include "Euler/Multiphysics/eulerreducedtests.h"
#include "Elasticity/elastictests.h"
#include "Elasticity/Multiphysics/elasticmultiphysicstests.h"
#include "Elasticity/Multiphysics/elasticreducedtests.h"
#include "MHD/mhdtests.h"
#include "HPR/hprtests.h"
#include "HPR/Multiphysics/hprmultiphysicstests.h"
#include "HPR/Multiphysics/hprintermediatetests.h"
#include "HPR/Multiphysics/hprreducedtests.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    //EulerTests::solveToroTest1(1000);
    //EulerTests::solve2DToroTest1(100);
    //EulerMultiphysicsTests::solveToroTest1(400, 3);
    //EulerMultiphysicsTests::solveFedkiwTest(100, 3);
    //EulerMultiphysicsTests::solveChinnayyaTest(100);

    //EulerMultiphysicsTests::solve2DToroTest1(100, 3);
    //EulerMultiphysicsTests::solve2DFedkiwTest(100, 3);

    //EulerReducedTests::solveToroTest1(100, 3);
    //EulerReducedTests::solveFedkiwTest(200, 3);
    //EulerReducedTests::solve2DToroTest1(100, 3);
    //EulerReducedTests::solve2DFedkiwTest(100, 3);

    //ElasticTests::solveBartonTest1(100);
    //ElasticTests::solveBartonTest2(200);
    //ElasticTests::solve2DBartonTest1(100);
    //ElasticTests::solve2DBartonTest2(20);

    //ElasticMultiphysicsTests::solveBartonTest1(100, 3);
    //ElasticMultiphysicsTests::solveBartonTest2(100, 3);
    //ElasticMultiphysicsTests::solveZhangTest(100, 0);

    //ElasticReducedTests::solveBartonTest1(100, 3);
    //ElasticReducedTests::solveZhangTest(800, 3);
    //ElasticReducedTests::solve2DBartonTest1(100, 3);
    //ElasticReducedTests::solve2DZhangTest(100, 3);
    //ElasticReducedTests::solve2DZhangTest2(100, 3);

    //MHDTests::solveDumbserTest2(400);
    //MHDTests::solveDumbserTest1(800);
    //MHDTests::solveDumbserTest3(400);

    //HPRTests::solveStokesFirstProblem(200, 1);
    //HPRTests::solveHeatConductionProblem(100, 20);
    //HPRTests::solveBartonTest1(100, 0);
    //HPRTests::solveZhangElasticPlasticTest(400, 1);
    //HPRTests::solve2DBartonTest1(100, 0);
    //HPRTests::solve2DZhangElasticPlasticTest(100, 1);

    //HPRMultiphysicsTests::solveBartonTest1(100, 3);
    HPRMultiphysicsTests::solveZhangTest(100, 0);

    //HPRIntermediateTests::solveBartonTest1(400, 3);
    //HPRIntermediateTests::solveZhangTest(400, 0);
    //HPRIntermediateTests::solveZhangElasticPlasticTest(400, 1, 3);
    //HPRIntermediateTests::solve2DBartonTest1(10, 3);

    //HPRReducedTests::solveBartonTest1(200, 0);
    //HPRReducedTests::solveZhangTest(800, 0);
    //HPRReducedTests::solveZhangElasticPlasticTest(100, 1, 0);
    //HPRReducedTests::solve2DBartonTest1(100, 0);
    //HPRReducedTests::solve2DZhangTest(100, 0);
    //HPRReducedTests::solve2DZhangTest2(100, 3);
    //HPRReducedTests::solve2DZhangElasticPlasticTest(100, 1, 0);
    //HPRReducedTests::solve2DZhangElasticPlasticTest2(100, 1, 0);

    return a.exec();
}

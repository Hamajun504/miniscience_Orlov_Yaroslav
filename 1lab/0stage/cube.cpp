#include <gmsh.h>

int main()
{
    gmsh::initialize();

    gmsh::model::add("cube");

    double lc = 1e-1;
    gmsh::model::geo::addPoint(0, 0, 0, lc, 1);
    gmsh::model::geo::addPoint(0, 0, 1, lc, 2);
    gmsh::model::geo::addPoint(0, 1, 1, lc, 3);
    gmsh::model::geo::addPoint(0, 1, 0, lc, 4);
    gmsh::model::geo::addPoint(1, 0, 0, lc, 5);
    gmsh::model::geo::addPoint(1, 0, 1, lc, 6);
    gmsh::model::geo::addPoint(1, 1, 1, lc, 7);
    gmsh::model::geo::addPoint(1, 1, 0, lc, 8);

    gmsh::model::geo::addLine(1, 2, 1);
    gmsh::model::geo::addLine(2, 3, 2);
    gmsh::model::geo::addLine(3, 4, 3);
    gmsh::model::geo::addLine(4, 1, 4);

    gmsh::model::geo::addLine(1, 5, 5);
    gmsh::model::geo::addLine(2, 6, 6);
    gmsh::model::geo::addLine(3, 7, 7);
    gmsh::model::geo::addLine(4, 8, 8);

    gmsh::model::geo::addLine(5, 6, 9);
    gmsh::model::geo::addLine(6, 7, 10);
    gmsh::model::geo::addLine(7, 8, 11);
    gmsh::model::geo::addLine(8, 5, 12);

    gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);
    gmsh::model::geo::addPlaneSurface({1}, 1);

    gmsh::model::geo::addCurveLoop({1, 6, -9, -5}, 2);
    gmsh::model::geo::addPlaneSurface({2}, 2);

    gmsh::model::geo::addCurveLoop({2, 7, -10, -6}, 3);
    gmsh::model::geo::addPlaneSurface({3}, 3);

    gmsh::model::geo::addCurveLoop({3, 8, -11, -7}, 4);
    gmsh::model::geo::addPlaneSurface({4}, 4);

    gmsh::model::geo::addCurveLoop({4, 5, -12, -8}, 5);
    gmsh::model::geo::addPlaneSurface({5}, 5);

    gmsh::model::geo::addCurveLoop({9, 10, 11, 12}, 6);
    gmsh::model::geo::addPlaneSurface({6}, 6);

    gmsh::model::geo::addSurfaceLoop({1, 2, 3, 4, 5, 6}, 1);
    gmsh::model::geo::addVolume({1});

    gmsh::model::geo::synchronize();

    gmsh::model::mesh::generate(3);

    gmsh::write("cube.msh");

    gmsh::fltk::run();

    gmsh::finalize();
}
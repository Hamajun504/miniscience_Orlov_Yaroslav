#include <gmsh.h>
#include <cmath>
#include <iostream>

int main()
{
    gmsh::initialize();

    gmsh::merge("Goose_full.stl");

    gmsh::model::mesh::createTopology();

    gmsh::model::mesh::classifySurfaces(180 * M_PI / 180, false, true, M_PI / 2);

    gmsh::model::mesh::createGeometry();

     std::vector<std::pair<int, int> > s;
    gmsh::model::getEntities(s, 2);
    std::vector<int> sl;
    for(auto surf : s) sl.push_back(surf.second);
    int l = gmsh::model::geo::addSurfaceLoop(sl);

    gmsh::model::geo::addVolume({l});

    gmsh::model::geo::synchronize();

    gmsh::model::mesh::field::add("MathEval", 3);
    gmsh::model::mesh::field::setString(3, "F", "3");
    gmsh::model::mesh::field::setAsBackgroundMesh(3);


    gmsh::model::mesh::generate(3);

    gmsh::write("goose.msh");

    gmsh::fltk::run();

    gmsh::finalize();

    return 0;
}
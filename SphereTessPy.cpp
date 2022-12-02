// C++ includes
#include <iostream>
#include <boost/python.hpp>
// CLHEP includes
//#include "CLHEP/Vector/ThreeVector.h"
// FreeCAD includes
#include "Base/Vector3D.h"

include "sphere_tesselation.hh"

using namespace std;
using namespace boost::python;

// ====================================================================
// Function definitions Class then Function
// ====================================================================
class sphereTesselation : public SphereTesselation
{

public:

MySphereTesselation::SphereTesselation(TESSELATION initial, size_t numSubdivisions)

private:
};


// Constructors
//MyG4TriangularFacet::MyG4TriangularFacet() {cout << "Default";}

//MyG4TriangularFacet::MyG4TriangularFacet(Base::Vector3d v0,
//                          Base::Vector3d v1,
//                          Base::Vector3d v2)

// {
// G4cout << "Facet Constructor FreeCAD : ";
// G4TriangularFacet(G4ThreeVector(v0.x,v0.y,v0.z),
//                  G4ThreeVector(v1.x,v1.y,v1.z),
//                  G4ThreeVector(v2.x,v2.y,v2.z),
//              ABSOLUTE);
//}

//G4TriangularFacet(v0,v1,v2,ABSOLUTE);
// Try the long way
//G4TriangularFacet facet;
//facet.SetVertex(0,v0);
//facet.SetVertex(1,v1);
//facet.SetVertex(2,v2);
//}

//MyG4TriangularFacet::MyG4TriangularFacet(CLHEP::HepVector& v0,
//                        CLHEP::HepVector& v1,
//                        CLHEP::HepVector& v2)

//{
//G4TriangularFacet(v0,v1,v2,ABSOLUTE);
//}


// ====================================================================
// module definition
// ====================================================================
BOOST_PYTHON_MODULE(SphereTess)
{
    using namespace boost::python;
    def("sphereTesselation", sphereTesselation);

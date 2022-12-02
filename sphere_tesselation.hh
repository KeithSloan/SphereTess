#ifndef SPHERE_TESSELATION_HH
#define SPHERE_TESSELATION_HH
#include "Vector3D.h"
#include <valarray>
#include <vector>

using namespace std;

/**
 * Subdivides a sphere into triangles of (approximately) equal area
 * The inital division into a THETRAHEDORn, or OCTAHEDRON, or ICOSAHEDRON
 * is "exact". The starting faces can then be subdivided into four faces each
 * by bisecting each edge. If the number of subdivisions \a numSubdivisions is
 * 0, you get the initial polyhedron
 *
 * For reference, the number of faces are:
 * Tetrahedron = 4 + 4*4*numSubdivisions;
 * Octahedron  = 8 + 8*4*numSubdivisions;
 * Icosahedron = 20 + 20*4*numSubdivisions;
 */

class SphereTesselation
{
public:
  enum TESSELATION {TETRAHEDRON, OCTAHEDRON, ICOSAHEDRON};
  
  SphereTesselation(TESSELATION initial, size_t numSubdivisions);
  ~SphereTesselation(){};
  
  int numVertices() const {return n_vertices;}
  int numFaces() const {return n_faces;}
  int numEdges() const {return n_edges;}
  Vector3D<double> faceMean(int n) const;
  Vector3D<double> vertex(int i) const {
    return ((i < n_faces) ? Vector3D<double>(m_vertices[3*i], m_vertices[3*i+1], m_vertices[3*i+2]) 
	    : Vector3D<double>(0.,0.,0.));
  }

  // return the index of the three vertexes making up this face
  bool face(int i, int *vertexIndex) const {
    if (i < n_faces) {
      vertexIndex[0] = m_faces[3*i];
      vertexIndex[1] = m_faces[3*i+1];
      vertexIndex[2] = m_faces[3*i+2];
      return true;
    }
    return false;
  }
  

  // return the index of the face closest to this one
  int closestFace(Vector3D<double> v) const;


private:
  void init_tetrahedron();
  void init_octahedron();
  void init_icosahedron();
  int search_midpoint (int index_start, int index_end,
		       valarray<int> &start,
		       valarray<int> &end,
		       valarray<int> &midpoint) ;
  void subdivide();

  // to speed searching for the containing a given vector we keep an
  // array of faces closest to a grid in costheta, phi. We then test only among
  // the faces in the array for that [ithet][iphi] bin
  void make_thet_phi_array();

  int closestFace(Vector3D<double> v, const vector<int> &faces) const ;

  
private:
  int n_vertices;
  int n_faces;
  int n_edges;
  valarray<float> m_vertices;
  valarray<int> m_faces; // for the ith face faces[3*i] = index of first vertex,
                       // faces[3*i+1] = index of second vertex;
                       // faces[3*i+2] = index of third vertex;

  // to help in finding which face covers a vector in 3D
  // we maintain a list of the faces which overlap a theta-phi 
  // division
  vector<int> thet_phi_faces[21][21];
};

#endif

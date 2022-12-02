#include "sphere_tesselation.hh"

static int edge_walk = 0; 

/******************************************************************/
SphereTesselation::SphereTesselation(TESSELATION initial, size_t numSubdivisions)
 {
   switch(initial) {
   case TETRAHEDRON:
     init_tetrahedron();
     break;
   case OCTAHEDRON:
     init_octahedron();
     break;
   case ICOSAHEDRON:
     init_icosahedron();
     break;
   }

   for (int i=0; i<numSubdivisions; i++)  {
     subdivide (); 
   }
   
   make_thet_phi_array();
 }

/******************************************************************/
void 
SphereTesselation::init_tetrahedron () 
{ 
  float sqrt3 = 1 / sqrt(3.0);
  float tetrahedron_vertices[] = {sqrt3, sqrt3, sqrt3,
				  -sqrt3, -sqrt3, sqrt3,
				  -sqrt3, sqrt3, -sqrt3,
				  sqrt3, -sqrt3, -sqrt3}; 
  int tetrahedron_faces[] = {0, 2, 1, 0, 1, 3, 2, 3, 1, 3, 2, 0};

  n_vertices = 4; 
  n_faces = 4; 
  n_edges = 6; 
  m_vertices.resize(3*n_vertices);
  m_vertices = valarray<float>(tetrahedron_vertices, 3*n_vertices);
  m_faces.resize(3*n_faces);
  m_faces =  valarray<int>(tetrahedron_faces,3*n_faces); 
} 

/******************************************************************/
void 
SphereTesselation::init_octahedron () 
{ 
  float octahedron_vertices[] = {0.0, 0.0, -1.0,
				 1.0, 0.0, 0.0,
				 0.0, -1.0, 0.0,
				 -1.0, 0.0, 0.0,
				 0.0, 1.0, 0.0,
				 0.0, 0.0, 1.0}; 
  int octahedron_faces[] = {0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 1, 5, 2, 1, 5, 3, 2, 5, 4, 3, 5, 1, 4}; 

  n_vertices = 6; 
  n_faces = 8;
  n_edges = 12; 
  m_vertices.resize(3*n_vertices);
  m_vertices = valarray<float>(octahedron_vertices, 3*n_vertices);
  m_faces.resize(3*n_faces);
  m_faces =  valarray<int>(octahedron_faces, 3*n_faces); 
} 

/******************************************************************/
void
SphereTesselation::init_icosahedron () 
{ 
  float t = (1+sqrt(5))/2;
  float tau = t/sqrt(1+t*t);
  float one = 1/sqrt(1+t*t);

  /*
  float icosahedron_vertices[] = {tau, one, 0.0,
				  -tau, one, 0.0,
				  -tau, -one, 0.0,
				  tau, -one, 0.0,
				  one, 0.0 ,  tau,
				  one, 0.0 , -tau,
				  -one, 0.0 , -tau,
				  -one, 0.0 , tau,
				  0.0 , tau, one,
				  0.0 , -tau, one,
				  0.0 , -tau, -one,
				  0.0 , tau, -one};
 int icosahedron_faces[] = {4, 8, 7,
			    4, 7, 9,
			    5, 6, 11,
			    5, 10, 6,
			    0, 4, 3,
			    0, 3, 5,
			    2, 7, 1,
			    2, 1, 6,
			    8, 0, 11,
			    8, 11, 1,
			    9, 10, 3,
			    9, 2, 10,
			    8, 4, 0,
			    11, 0, 5,
			    4, 9, 3,
			    5, 3, 10,
			    7, 8, 1,
			    6, 1, 11,
			    7, 2, 9,
			    6, 10, 2};
  */

  float icosahedron_vertices[] = {0, 0, 1/one,
				  t, one, tau,
				  0, 2*tau, (-1+t*t)*one,
				  -t, one, tau,
				  -1, -t*tau ,  tau,
				  1, -t*t*one, tau,
				  -1, t*t*one, -tau,
				  1, t*tau , -tau,
				  t , -one, -tau,
				  0.0 , -2*tau, (1-t*t)*one,
				  -t, -one, -tau,
				  0.0 , 0.0, -1/one};
 int icosahedron_faces[] = {0, 1, 2,
			    0, 2, 3,
			    0, 3, 4,
			    0, 4, 5,
			    0, 5, 1,
			    5, 8, 1,
			    1, 7, 2,
			    2, 6, 3,
			    3, 10, 4,
			    4, 9, 5,

			    7, 2, 6,
			    8, 1, 7,
			    9, 5, 8,
			    10, 4, 9,
			    6, 3, 10,
			    11, 6, 10,
			    11, 7, 6,
			    11, 8, 7,
			    11, 9, 8,
			    11, 10, 9};

 
  n_vertices = 12; 
  n_faces = 20;
  n_edges = 30;
  m_vertices.resize(3*n_vertices);
  m_vertices = valarray<float>(icosahedron_vertices, 3*n_vertices);
  m_faces.resize(3*n_faces);
  m_faces =  valarray<int>(icosahedron_faces, 3*n_faces); 
} 

/******************************************************************/
int 
SphereTesselation::search_midpoint (int index_start, int index_end,
				    valarray<int> &start,
				    valarray<int> &end,
				    valarray<int> &midpoint) 
{ 
  int i;
  for (i=0; i<edge_walk; i++) 
    if ((start[i] == index_start && end[i] == index_end) || 
	(start[i] == index_end && end[i] == index_start)) 
      {
	int res = midpoint[i];

	/* update the arrays */
	start[i]    = start[edge_walk-1];
	end[i]      = end[edge_walk-1];
	midpoint[i] = midpoint[edge_walk-1];
	edge_walk--;
	
	return res; 
      }

  /* vertex not in the list, so we add it */
  start[edge_walk] = index_start;
  end[edge_walk] = index_end; 
  midpoint[edge_walk] = n_vertices; 
  
  /* create new vertex */ 
  m_vertices[3*n_vertices]   = (m_vertices[3*index_start] + m_vertices[3*index_end]) / 2.0;
  m_vertices[3*n_vertices+1] = (m_vertices[3*index_start+1] + m_vertices[3*index_end+1]) / 2.0;
  m_vertices[3*n_vertices+2] = (m_vertices[3*index_start+2] + m_vertices[3*index_end+2]) / 2.0;
  
  /* normalize the new vertex */ 
  float length = sqrt (m_vertices[3*n_vertices] * m_vertices[3*n_vertices] +
		       m_vertices[3*n_vertices+1] * m_vertices[3*n_vertices+1] +
		       m_vertices[3*n_vertices+2] * m_vertices[3*n_vertices+2]);
  length = 1/length;
  m_vertices[3*n_vertices] *= length;
  m_vertices[3*n_vertices+1] *= length;
  m_vertices[3*n_vertices+2] *= length;
  
  n_vertices++;
  edge_walk++;
  return midpoint[edge_walk-1];
} 
/******************************************************************/
void 
SphereTesselation::subdivide() 
{ 
  int n_vertices_new = n_vertices+2*n_edges; 
  int n_faces_new = 4*n_faces; 
  int i; 

  n_edges = 2*n_vertices + 3*n_faces; 
  valarray<int> start(n_edges); 
  valarray<int> end(n_edges); 
  valarray<int> midpoint(n_edges); 

  valarray<int> faces_old(m_faces);
  valarray<float> temp_vertices(m_vertices);
  m_vertices.resize(3*n_vertices_new);
  for(size_t i=0; i < temp_vertices.size(); i++) 
    m_vertices[i] = temp_vertices[i];

  m_faces.resize(3*n_faces_new);

  n_faces_new = 0; 

  for (i=0; i<n_faces; i++) 
    { 
      int a = faces_old[3*i]; 
      int b = faces_old[3*i+1]; 
      int c = faces_old[3*i+2]; 

      int ab_midpoint = search_midpoint (b, a, start, end, midpoint); 
      int bc_midpoint = search_midpoint (c, b, start, end, midpoint); 
      int ca_midpoint = search_midpoint (a, c, start, end, midpoint); 

      m_faces[3*n_faces_new] = a; 
      m_faces[3*n_faces_new+1] = ab_midpoint; 
      m_faces[3*n_faces_new+2] = ca_midpoint; 
      n_faces_new++; 
      m_faces[3*n_faces_new] = ca_midpoint; 
      m_faces[3*n_faces_new+1] = ab_midpoint; 
      m_faces[3*n_faces_new+2] = bc_midpoint; 
      n_faces_new++; 
      m_faces[3*n_faces_new] = ca_midpoint; 
      m_faces[3*n_faces_new+1] = bc_midpoint; 
      m_faces[3*n_faces_new+2] = c; 
      n_faces_new++; 
      m_faces[3*n_faces_new] = ab_midpoint; 
      m_faces[3*n_faces_new+1] = b; 
      m_faces[3*n_faces_new+2] = bc_midpoint; 
      n_faces_new++; 
    } 
  n_faces = n_faces_new; 
} 

/******************************************************************/
Vector3D<double>
SphereTesselation::faceMean(int n) const
{ 
  Vector3D<double> defaultvector(0.,0.,0.);

  if(n < n_faces) {
    int i1 = m_faces[3*n];
    int i2 = m_faces[3*n+1];
    int i3 = m_faces[3*n+2]; 
    
    Vector3D<double> v1(m_vertices[3*i1],m_vertices[3*i1+1],m_vertices[3*i1+2]); 
    Vector3D<double> v2(m_vertices[3*i2],m_vertices[3*i2+1],m_vertices[3*i2+2]); 
    Vector3D<double> v3(m_vertices[3*i3],m_vertices[3*i3+1],m_vertices[3*i3+2]); 
    return (v1+v2+v3)/3;
  }

  return defaultvector;
}

/******************************************************************/
// Helper routine to make the theta_phi array of faces
void
SphereTesselation::make_thet_phi_array() 
{
  double dcos = 2./20;
  double dphi = 2*M_PI/20;

  vector<int> faces;
  for(int i=0; i < n_faces; i++) 
    faces.push_back(i);

  for(int i=0; i < n_faces; i++) {
    Vector3D<double> v = faceMean(i);
    int iTheta = (v.cosTheta() +1)/dcos;
    int iPhi = (v.phi()+M_PI)/dphi; // CLHEP's phi is from -Pi to +PI, not 0 to 2*pi
    thet_phi_faces[iTheta][iPhi].push_back(i);
  }

  // there may be bins without faces assigned to them
  // for those assign the closest face
  for(int i=0; i < 20; i++) {
    for(int j=0; j < 20; j++) {
      if(thet_phi_faces[i][j].size() == 0) {	
	// search for face closest to this bin
	double costhet = -1 + dcos*(i+0.5);
	double sinthet = costhet >= 1 ? 0 : sqrt(1.-costhet*costhet);
	double phi = (j+0.5)*dphi;
	Vector3D<double> v(sinthet*cos(phi), sinthet*sin(phi), costhet);
	thet_phi_faces[i][j].push_back(closestFace(v,faces)); 
      }
    }
  }
}

/******************************************************************/
int 
SphereTesselation::closestFace(Vector3D<double> v, const vector<int> &faces) const
{
  double cosmax = -3;
  int nface = 0;
  for(int i=0; i < faces.size(); i++) {
    int k = faces[i];
    Vector3D<double> vmid = faceMean(k).unit();
    double cosine = vmid.dot(v);
    if(cosine > cosmax) {
      cosmax = cosine;
      nface = k;
    } 
  }
  
  return nface;
}
/******************************************************************/
int 
SphereTesselation::closestFace(Vector3D<double> v) const
{
  double dcos = 2./20;
  double dphi = 2*M_PI/20;
  
  int iTheta = (v.cosTheta() +1)/dcos;
  int iPhi = (v.phi()+M_PI)/dphi;
  const vector<int> &faces = thet_phi_faces[iTheta][iPhi];

  return closestFace(v, faces);
}
#if 0
/******************************************************************/
int 
SphereTesselation::closestFace(Vector3D<double> v) const
{
  vector<int> faces;
  for(int i=0; i < n_faces; i++) 
    faces.push_back(i);

  return closestFace(v, faces);
}
#endif


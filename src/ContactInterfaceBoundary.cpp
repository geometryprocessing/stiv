#include "ContactInterfaceBoundary.h"

ContactInterfaceBoundary::ContactInterfaceBoundary()
{
  m_pos = NULL;
  mesh_manager = NULL;
  m_P = NULL;
}


ContactInterfaceBoundary::~ContactInterfaceBoundary()
{
  if(m_pos)
      delete[] m_pos;
  //mesh_manager->clear();
  if(mesh_manager)
      delete mesh_manager;
  //delete m_P;
}

void ContactInterfaceBoundary::getVolumeAndGradient(std::vector<double> &IV, std::vector<double> &IS_BDV, int &IV_size, 
        std::vector<double> &Grad, std::vector<int> &Grad_index,
        std::vector<double> &pos_s, std::vector<double> &pos_e, double *pos_s_pole, double *pos_e_pole,
        double minSep, double* pos_bd, int nv_resident, double periodicLength)
{
  mesh_manager->getMesh()->setThickness(minSep);
  mesh_manager->getMesh()->setNvResident(nv_resident);
  mesh_manager->getMesh()->setPeriodicLength(periodicLength);
  IV.clear();
  IV_size = 0;

  int vcount = 0;
  fillMPos(pos_s, pos_s_pole);
  for (vector<set<unsigned int> >::iterator vItr = mesh_manager->getMesh()->m_subMeshVertices.begin(); 
      vItr != mesh_manager->getMesh()->m_subMeshVertices.end(); ++vItr)
  {
      if(vcount < m_N)
      {
          int i = 0;
          for (set<unsigned int>::iterator sItr=vItr->begin(); sItr!=vItr->end(); ++sItr)
          {
              mesh_manager->getMesh()->getSlaveMesh()->getVertex(*sItr)->Pori() = Point3d(m_pos[i + vItr->size()*0 + 3*vcount],
                                                                                          m_pos[i + vItr->size()*1 + 3*vcount],
                                                                                          m_pos[i + vItr->size()*2 + 3*vcount]);
              i++;
          }
      }
      else
      {
          int i = 0;
          for (set<unsigned int>::iterator sItr=vItr->begin(); sItr!=vItr->end(); ++sItr)
          {
              mesh_manager->getMesh()->getSlaveMesh()->getVertex(*sItr)->Pori() = Point3d(pos_bd[(i+vcount-m_N)*3 + 0],
                                                                                          pos_bd[(i+vcount-m_N)*3 + 1],
                                                                                          pos_bd[(i+vcount-m_N)*3 + 2]);
              i++;
          }
      }
      vcount = vcount + vItr->size();
  }

  /*
  vcount = 0;
  fillMPos(pos_s, pos_s_pole);
  for (vector<set<unsigned int> >::iterator vItr = mesh_manager->getMesh()->m_subMeshVertices.begin(); 
      vItr != mesh_manager->getMesh()->m_subMeshVertices.end(); ++vItr)
  {
    int i = 0;
    for (set<unsigned int>::iterator sItr=vItr->begin(); sItr!=vItr->end(); ++sItr)
    {
      mesh_manager->getMesh()->getSlaveMesh()->getVertex(*sItr)->P() = Point3d(m_pos[i + vItr->size()*0 + 3*vcount],
                                                                               m_pos[i + vItr->size()*1 + 3*vcount],
                                                                               m_pos[i + vItr->size()*2 + 3*vcount]);

      i++;
    }
    vcount = vcount + vItr->size();
  }

  std::cout<< "Begin to calculate minDistance "<< std::endl;
  double t_md_start = omp_get_wtime();
  double minD = mesh_manager->getMesh()->minDistance();
  std::cout<< "start minDistance: " << minD << std::endl;
  double t_md_end = omp_get_wtime();
  std::cout<< "End to calculate minDistance, time used: "<< (t_md_end-t_md_start) << std::endl;
  assert(minD > minSep);
  */
  
  vcount = 0;
  fillMPos(pos_e, pos_e_pole);
  for (vector<set<unsigned int> >::iterator vItr = mesh_manager->getMesh()->m_subMeshVertices.begin(); 
      vItr != mesh_manager->getMesh()->m_subMeshVertices.end(); ++vItr)
  {
      if(vcount < m_N)
      {
          int i = 0;
          for (set<unsigned int>::iterator sItr=vItr->begin(); sItr!=vItr->end(); ++sItr)
          {
              mesh_manager->getMesh()->getSlaveMesh()->getVertex(*sItr)->P() = Point3d(m_pos[i + vItr->size()*0 + 3*vcount],
                                                                                       m_pos[i + vItr->size()*1 + 3*vcount],
                                                                                       m_pos[i + vItr->size()*2 + 3*vcount]);
              i++;
          }
      }
      else
      {
          int i = 0;
          for (set<unsigned int>::iterator sItr=vItr->begin(); sItr!=vItr->end(); ++sItr)
          {
              mesh_manager->getMesh()->getSlaveMesh()->getVertex(*sItr)->P() = Point3d(pos_bd[(i+vcount-m_N)*3 + 0],
                                                                                       pos_bd[(i+vcount-m_N)*3 + 1],
                                                                                       pos_bd[(i+vcount-m_N)*3 + 2]);
              i++;
          }
      }
      vcount = vcount + vItr->size();
  }

  /*
  std::cout<< "Begin to calculate minDistance "<< std::endl;
  t_md_start = omp_get_wtime();
  minD = mesh_manager->getMesh()->minDistance();
  std::cout<< "end minDistance: " << minD << std::endl;
  t_md_end = omp_get_wtime();
  std::cout<< "End to calculate minDistance, time used: "<< (t_md_end-t_md_start) << std::endl;
  */
  //double minDPori = mesh_manager->getMesh()->minDistanceBD(true);
  //double minDP = mesh_manager->getMesh()->minDistanceBD(false);
  
  std::cout<<std::scientific<<std::setprecision(16);
  //std::cout<< "start minDistance: " << minDPori << std::endl;
  //std::cout<< "end minDistance: " << minDP << std::endl;
  //assert(minDPori > minSep);

  // TODO: add to memember and constructor
  // set all vertices as seleted
  set<unsigned int> selected;
  /*
  for(int i=0; i<mesh_manager->getMesh()->getSlaveMesh()->size_of_vertices();++i)
  {
    selected.insert(i);
  }
  */
  
  // TODO: add to memember and constructor, just need to call setzeros
  // allocate gradient
  Eigen::SparseVector<double> gradient(3*mesh_manager->getMesh()->getSlaveMesh()
      ->size_of_vertices());
  // allocate gradient index
  Eigen::SparseVector<int> gradient_index(3*mesh_manager->getMesh()->getSlaveMesh()
      ->size_of_vertices());
  
  // set parameters for get gradient and volume
  std::map<unsigned int, unsigned int> m_fadeCounter;
  Vector_3 mov(0,0,0);
  bool handleCollisions = true;
  bool handleElasticity = false;
  bool rigid = false;
  bool moveAll = true;
  int limiter = -1;
  
  ////std::cout<< "Begin to calculate gradient "<< std::endl;
  // get volume gradient and intersection volume
  double t_check_start = omp_get_wtime();
  mesh_manager->getMesh()->controlMeshUpdate(selected, m_fadeCounter, handleCollisions, 
      handleElasticity, mov, rigid, moveAll, limiter, 
      gradient, IV, gradient_index, IV_size, IS_BDV);
  double t_check_end = omp_get_wtime();
  ////std::cout<< "End to calculate gradient, time used: "<< (t_check_end-t_check_start) << std::endl;

  // TODO: iterator nonzero element
  // map back to position without two poles
  //vcount = 0;
  int Ni = m_Ni-2;
  int Dim = 3;
  int m = 2*m_p;
  #pragma omp parallel for
  for(int i=0; i<m_nv;++i)
  {
    //#pragma omp parallel for
    for(int j=0; j<m_Ni; ++j)
    {
      if( (j != 0) && (j != (m_Ni-1)) )
      {
        int vcounti = i*m_Ni + j;
        Grad[j-1 + Ni*0 + i*Ni*Dim] = gradient.coeff(3*vcounti + 0);
        Grad[j-1 + Ni*1 + i*Ni*Dim] = gradient.coeff(3*vcounti + 1);
        Grad[j-1 + Ni*2 + i*Ni*Dim] = gradient.coeff(3*vcounti + 2);
        
        Grad_index[j-1 + Ni*0 + i*Ni*Dim] = gradient_index.coeff(3*vcounti + 0);
        Grad_index[j-1 + Ni*1 + i*Ni*Dim] = gradient_index.coeff(3*vcounti + 1);
        Grad_index[j-1 + Ni*2 + i*Ni*Dim] = gradient_index.coeff(3*vcounti + 2);
      }
      //vcount++;
    }

    // north pole
    int vcounti = i*m_Ni + 0;
    if(gradient_index.coeff(3*vcounti + 0) > 0 || 
            gradient_index.coeff(3*vcounti + 1) > 0 ||
            gradient_index.coeff(3*vcounti + 2) > 0 )
    {
        for(int ii=0;ii<m;++ii)
        {
            Grad[ii + 0*Ni + i*Ni*Dim] += gradient.coeff(3*vcounti + 0)/m;
            Grad[ii + 1*Ni + i*Ni*Dim] += gradient.coeff(3*vcounti + 1)/m;
            Grad[ii + 2*Ni + i*Ni*Dim] += gradient.coeff(3*vcounti + 2)/m;
                  
            if(Grad_index[ii + 0*Ni + i*Ni*Dim] > 0)
                continue;

            Grad_index[ii + 0*Ni + i*Ni*Dim] = gradient_index.coeff(3*vcounti + 0);
            Grad_index[ii + 1*Ni + i*Ni*Dim] = gradient_index.coeff(3*vcounti + 1);
            Grad_index[ii + 2*Ni + i*Ni*Dim] = gradient_index.coeff(3*vcounti + 2);
        }
    }

    // south pole
    vcounti = i*m_Ni + m_Ni -1;
    if(gradient_index.coeff(3*vcounti + 0) > 0 || 
            gradient_index.coeff(3*vcounti + 1) > 0 ||
            gradient_index.coeff(3*vcounti + 2) > 0 )
    {
        for(int ii=0;ii<m;++ii)
        {
            Grad[Ni-m + ii + 0*Ni + i*Ni*Dim] += gradient.coeff(3*vcounti + 0)/m;
            Grad[Ni-m + ii + 1*Ni + i*Ni*Dim] += gradient.coeff(3*vcounti + 1)/m;
            Grad[Ni-m + ii + 2*Ni + i*Ni*Dim] += gradient.coeff(3*vcounti + 2)/m;
            
            if(Grad_index[Ni-m + ii + 0*Ni + i*Ni*Dim] > 0)
                continue;

            Grad_index[Ni-m + ii + 0*Ni + i*Ni*Dim] = gradient_index.coeff(3*vcounti + 0);
            Grad_index[Ni-m + ii + 1*Ni + i*Ni*Dim] = gradient_index.coeff(3*vcounti + 1);
            Grad_index[Ni-m + ii + 2*Ni + i*Ni*Dim] = gradient_index.coeff(3*vcounti + 2);
        }
    }

  }
  
  std::cout<<"IV size is: "<<IV.size()<<std::endl;
  //if(IV.size() > 0)
  if(false)
  {
    static int count_i = 0;
    std::cout<< "Begin write to vtk "<< std::endl;
    // write vtk file
    std::ofstream myfile;
    std::stringstream fname1;
    fname1<<"data/testP_";
    fname1<<std::setfill('0')<<std::setw(5)<<count_i;
    fname1<<".vtk";
    myfile.open(fname1.str().c_str());
    myfile << "# vtk DataFile Version 2.0\n";
    myfile << "Unstructured Grid Example\n";
    myfile << "ASCII\n";
    myfile << "DATASET UNSTRUCTURED_GRID\n";
    
    myfile << "POINTS ";
    myfile << m_N+m_Np;
    myfile << " double\n";
  
    //myfile << m_N << " " << m_faces.size() << " " << m_edges.size() <<std::endl;
    //for(int iv=0;iv<m_nv;++iv)
    //{
    //  for(int i=0; i<m_Ni; ++i)
    //  {
    //    myfile << m_pos[i + m_Ni*0 + 3*iv*m_Ni] << " " << 
    //              m_pos[i + m_Ni*1 + 3*iv*m_Ni] << " " << 
    //              m_pos[i + m_Ni*2 + 3*iv*m_Ni] << std::endl;;
    //  }
    //}
    
    Vertex_iterator vit = mesh_manager->getMesh()->getControlMesh()->vertices_begin();
    while (vit != mesh_manager->getMesh()->getControlMesh()->vertices_end()) 
    {
      Point3d pt = vit->P();
      myfile << pt[0]<< " " << pt[1] << " " << pt[2] << std::endl;;
      vit++;
    }
    
    //for(int i=0; i<m_N; ++i)
    //{
    //  myfile << m_pos[i] << " " << m_pos[i + m_Ni*1] << " " << m_pos[i + m_Ni*2] << std::endl;;
    //}
    
    myfile << "CELLS ";
    myfile << m_faces.size();
    myfile << " ";
    myfile << 4*m_faces.size() << std::endl;
    for(int i=0; i<m_faces.size(); ++i)
    {
      myfile << "3 " << m_faces[i].ai << " " << m_faces[i].bi << " " <<
        m_faces[i].ci << std::endl;
    }
    
    myfile << "CELL_TYPES ";
    myfile << m_faces.size() << std::endl;
    for(int i=0; i<m_faces.size(); ++i)
      myfile << "5\n";
  
    myfile << "POINT_DATA ";
    myfile << m_N+m_Np << std::endl;
    myfile << "VECTORS gradient double\n";
  
    for(int i=0; i<m_N+m_Np; ++i)
    {
      myfile << gradient.coeffRef(3*i + 0) << " " << gradient.coeffRef(3*i + 1) << 
        " " << gradient.coeffRef(3*i + 2) << std::endl;;
      //myfile << 1 << " " << 1 << 
        //" " << 1 << std::endl;;
    }
    myfile.close();
    // end of write vtk file

    // write vtk file
    std::stringstream fname2;
    fname2<<"data/testPori_";
    fname2<<std::setfill('0')<<std::setw(5)<<count_i;
    fname2<<".vtk";
    myfile.open(fname2.str().c_str());
    myfile << "# vtk DataFile Version 2.0\n";
    myfile << "Unstructured Grid Example\n";
    myfile << "ASCII\n";
    myfile << "DATASET UNSTRUCTURED_GRID\n";
    
    myfile << "POINTS ";
    myfile << m_N+m_Np;
    myfile << " double\n";
  
    //myfile << m_N << " " << m_faces.size() << " " << m_edges.size() <<std::endl;
    //for(int iv=0;iv<m_nv;++iv)
    //{
    //  for(int i=0; i<m_Ni; ++i)
    //  {
    //    myfile << m_pos[i + m_Ni*0 + 3*iv*m_Ni] << " " << 
    //              m_pos[i + m_Ni*1 + 3*iv*m_Ni] << " " << 
    //              m_pos[i + m_Ni*2 + 3*iv*m_Ni] << std::endl;;
    //  }
    //}
    
    vit = mesh_manager->getMesh()->getControlMesh()->vertices_begin();
    while (vit != mesh_manager->getMesh()->getControlMesh()->vertices_end()) 
    {
      Point3d pt = vit->Pori();
      myfile << pt[0]<< " " << pt[1] << " " << pt[2] << std::endl;;
      vit++;
    }
    
    //for(int i=0; i<m_N; ++i)
    //{
    //  myfile << m_pos[i] << " " << m_pos[i + m_Ni*1] << " " << m_pos[i + m_Ni*2] << std::endl;;
    //}
    
    myfile << "CELLS ";
    myfile << m_faces.size();
    myfile << " ";
    myfile << 4*m_faces.size() << std::endl;
    for(int i=0; i<m_faces.size(); ++i)
    {
      myfile << "3 " << m_faces[i].ai << " " << m_faces[i].bi << " " <<
        m_faces[i].ci << std::endl;
    }
    
    myfile << "CELL_TYPES ";
    myfile << m_faces.size() << std::endl;
    for(int i=0; i<m_faces.size(); ++i)
      myfile << "5\n";
  
    myfile << "POINT_DATA ";
    myfile << m_N+m_Np << std::endl;
    myfile << "VECTORS gradient double\n";
  
    for(int i=0; i<m_N+m_Np; ++i)
    {
      myfile << gradient.coeffRef(3*i + 0) << " " << gradient.coeffRef(3*i + 1) << 
        " " << gradient.coeffRef(3*i + 2) << std::endl;;
      //myfile << 1 << " " << 1 << 
        //" " << 1 << std::endl;;
    }
    myfile.close();
    // end of write vtk file
    std::cout<< "End write to vtk "<< std::endl;
    //count_i++;
  }
}

void ContactInterfaceBoundary::init(SURF_TYPE st)
{
  mesh_manager = new MeshManager();
  m_P = new Polyhedron;
  //std::ifstream fin("data/tmp.off", std::ios_base::in | std::ios_base::binary);
  //fin >> *m_P;
  Build_surface_boundary<Polyhedron::HalfedgeDS> triangles(this);
  double t_con_start = omp_get_wtime();
  m_P->delegate(triangles);
  double t_con_end = omp_get_wtime();
  std::cout<<"gen con time: "<<t_con_end-t_con_start<<"\n";

  m_P->setFilename("");
  m_P->setName("");
  m_P->init();
  mesh_manager->setMesh(m_P,st, m_nv, m_p, m_np, m_pp);
  
  //fin.close();
 }

/*
void ContactInterface::writeOFF()
{
  std::ofstream myfile;
  myfile.open("data/tmp.off");
  myfile << "OFF\n";
  myfile << m_N << " " << m_faces.size() << " " << m_edges.size() <<std::endl;

  for(int iv=0;iv<m_nv;++iv)
  {
    for(int i=0; i<m_Ni; ++i)
    {
      myfile << m_pos[i + m_Ni*0 + 3*iv*m_Ni] << " " << 
                m_pos[i + m_Ni*1 + 3*iv*m_Ni] << " " << 
                m_pos[i + m_Ni*2 + 3*iv*m_Ni] << std::endl;;
    }
  }
  

  for(int i=0; i<m_faces.size(); ++i)
  {
    myfile << "3 " << m_faces[i].ai << " " << m_faces[i].bi << " " <<
      m_faces[i].ci << std::endl;
  }

  myfile.close();
}
*/

void ContactInterfaceBoundary::fillMPos(std::vector<double> &pos, double *pos_pole)
{
  
  int n = m_p+1; // number of points in latitude
  int m = 2*m_p; // number of points in longitude
  
  // number of points in pos per object
  int Ni = n*m;
  int Dim = 3;

  //begin to fill m_pos
  #pragma omp parallel for
  for(int iv=0;iv<m_nv;++iv)
  {
    double xn = 0, yn = 0, zn = 0;
    double xs = 0, ys = 0, zs = 0;

    xn = pos_pole[iv*3*2+2*0+0];
    xs = pos_pole[iv*3*2+2*0+1];

    yn = pos_pole[iv*3*2+2*1+0];
    ys = pos_pole[iv*3*2+2*1+1];

    zn = pos_pole[iv*3*2+2*2+0];
    zs = pos_pole[iv*3*2+2*2+1];

    /*
    for(int i=0;i<m;++i)
    {
      xn += pos[i + 0*Ni + iv*Dim*Ni];
      yn += pos[i + 1*Ni + iv*Dim*Ni];
      zn += pos[i + 2*Ni + iv*Dim*Ni];

      xs += pos[Ni-m + i + 0*Ni + iv*Dim*Ni];
      ys += pos[Ni-m + i + 1*Ni + iv*Dim*Ni];
      zs += pos[Ni-m + i + 2*Ni + iv*Dim*Ni];
    }
     xn = xn/m;
     yn = yn/m;
     zn = zn/m;

     xs = xs/m;
     ys = ys/m;
     zs = zs/m;
     */
     
     m_pos[0 + 0*m_Ni + iv*Dim*m_Ni] = xn;
     m_pos[0 + 1*m_Ni + iv*Dim*m_Ni] = yn;
     m_pos[0 + 2*m_Ni + iv*Dim*m_Ni] = zn;
     
     m_pos[m_Ni-1 + 0*m_Ni + iv*Dim*m_Ni] = xs;
     m_pos[m_Ni-1 + 1*m_Ni + iv*Dim*m_Ni] = ys;
     m_pos[m_Ni-1 + 2*m_Ni + iv*Dim*m_Ni] = zs;

     //#pragma omp parallel for
     for(int j=0;j<Ni;++j)
     {
       m_pos[j+1 + 0*m_Ni + iv*Dim*m_Ni] =
         pos[j   + 0*Ni   + iv*Dim*Ni];
       
       m_pos[j+1 + 1*m_Ni + iv*Dim*m_Ni] =
         pos[j   + 1*Ni   + iv*Dim*Ni];

       m_pos[j+1 + 2*m_Ni + iv*Dim*m_Ni] =
         pos[j   + 2*Ni   + iv*Dim*Ni];
     }
  }
  //end of fill m_pos
}

void ContactInterfaceBoundary::generateMesh(std::vector<double> pos, double *pos_pole, int p, int nv, double* pos_bd, int pp, int np)
{
  // boundary dimension
  m_pp = pp;
  m_np = np;
  m_Nip = m_pp*m_pp;
  m_Np = m_np*m_Nip;
  m_pos_bd = pos_bd;

  // p+1 points in latitude
  // (2p) points in longitude
  // total 2+(2p)*(p+1) points(two poles are added)
  m_p = p;
  m_nv = nv;
  
  int n = m_p+1; // number of points in latitude
  int m = 2*m_p; // number of points in longitude
 
  // number of points in pos per object
  int Ni = n*m;
  int Dim = 3;

  m_Ni = Ni+2;
  m_N = m_Ni*m_nv; 

  m_pos = new double[m_N*Dim];
  // number of points in m_pos per object

  fillMPos(pos, pos_pole);

  // fill out index
  int *index = new int[m_N];
  for(int i=0;i<m_N;++i)
    index[i] = i;

#ifdef VERBOSE
  std::cout<<"total num of points: "<<m_N<<std::endl;
  for(int i=0;i<m_N;i++)
    std::cout<<"index:"<<index[i]<<std::endl;
#endif

  // begin to fill m_edges and m_faces
  for(int iv=0;iv<m_nv;++iv)
  {
    /* An example of p=4, with 5 points in latitude
     * and 10 points in longitude.
     * each '.' represents a discretization point,
     * each '-', '|' or '\' represents an edge,
     * and any three closed edges form a triangle face. 
     *
     *
     *             .
     *   | | | | | | | | | |
     *   .-.-.-.-.-.-.-.-.-.-
     *   |\|\|\|\|\|\|\|\|\|\
     *   .-.-.-.-.-.-.-.-.-.-
     *   |\|\|\|\|\|\|\|\|\|\
     *   .-.-.-.-.-.-.-.-.-.-
     *   | | | | | | | | | |
     *             .
     *
     */ 

    /*
     * Suppose points are sorted in the order as following
     * for p=4
     *
     *                  0
     *   1  2  3  4  5  6  7  8  9  10
     *   11 12 13 14 15 16 17 18 19 20
     *   21 22 23 24 25 26 27 28 29 30
     *                  31
     */

    // construct edges [i-j] i,j are point index

    // vertical edges
    // with north pole
    for(int i=1;i<=m;++i)
      m_edges.push_back(std::make_pair(0 + iv*m_Ni, i + iv*m_Ni));
    // with south pole
    for(int i=1;i<=m;++i)
      m_edges.push_back(std::make_pair(m_Ni-1 + iv*m_Ni,m_Ni-1-i + iv*m_Ni));
    // middle 
    for(int i=1;i<n;++i)
      for(int j=0;j<m;++j)
        m_edges.push_back(std::make_pair(1+(i-1)*m+j + iv*m_Ni,1+i*m+j + iv*m_Ni));
  
#ifdef VERBOSE
    for(std::vector< std::pair<int,int> >::iterator 
    it = m_edges.begin() ; it != m_edges.end(); ++it)
      std::cout <<"vertial edge pair: "<< it->first<<"-"<<it->second<<std::endl;
    m_edges.clear();
#endif

    // horizontal edges
    for(int i=1;i<=n;++i)
    {
      for(int j=0;j<m-1;++j)
        m_edges.push_back(std::make_pair( (i-1)*m+j+1 + iv*m_Ni, (i-1)*m+j+2 + iv*m_Ni));
      int j = m-1;
      m_edges.push_back(std::make_pair( (i-1)*m+j+1 + iv*m_Ni, (i-1)*m+(j+2)%m + iv*m_Ni));
    }

#ifdef VERBOSE
    for(std::vector< std::pair<int,int> >::iterator 
    it = m_edges.begin() ; it != m_edges.end(); ++it)
      std::cout <<"horizontal edge pair: "<< it->first<<"-"<<it->second<<std::endl;
    m_edges.clear();
#endif

    // diagonal edges
    for(int i=1;i<n;++i)
    {
      for(int j=0;j<m-1;++j)
        m_edges.push_back(std::make_pair(1+(i-1)*m+j + iv*m_Ni, 2+i*m+j + iv*m_Ni));
      int j = m-1;  
      m_edges.push_back(std::make_pair( 1+(i-1)*m+j + iv*m_Ni, i*m+(j+2)%m + iv*m_Ni));
    }

#ifdef VERBOSE
    for(std::vector< std::pair<int,int> >::iterator 
    it = m_edges.begin() ; it != m_edges.end(); ++it)
      std::cout <<"diagnoal edge pair: "<< it->first<<"-"<<it->second<<std::endl;
    m_edges.clear();
#endif

    // construct faces [i-j-k-] i,j,k are point index
    // original tri-faces
    // with north pole
    for(int i=1;i<m;++i)
      m_faces.push_back(TrigBD<int>(0 + iv*m_Ni, i + iv*m_Ni , i+1 + iv*m_Ni));
    m_faces.push_back(TrigBD<int>(0 + iv*m_Ni, m + iv*m_Ni, 1 + iv*m_Ni));
    // with south pole
    for(int i=1;i<m;++i)
      //faces.push_back(TrigBD<int>(N-1,N-1-m+i-1,N-1-m+i));
      m_faces.push_back(TrigBD<int>(m_Ni-1 + iv*m_Ni, m_Ni-1-m+i + iv*m_Ni, 
            m_Ni-1-m+i-1 + iv*m_Ni));
    //faces.push_back(TrigBD<int>(N-1,N-2,N-1-m));
    m_faces.push_back(TrigBD<int>(m_Ni-1 + iv*m_Ni, m_Ni-1-m + iv*m_Ni,
          m_Ni-2 + iv*m_Ni));

#ifdef VERBOSE
    for(std::vector< TrigBD<int> >::iterator 
    it = m_faces.begin() ; it != m_faces.end(); ++it)
      std::cout <<"original tri-face: "<< it->ai <<"-"<< it->bi <<
        "-"<< it->ci <<std::endl;
    m_faces.clear();
#endif

    // splitted tri-faces
    for(int i=1;i<n;++i)
    {
      for(int j=0;j<m-1;++j)
      {
        m_faces.push_back(
            TrigBD<int>( 1+(i-1)*m+j + iv*m_Ni, 1+i*m+j + iv*m_Ni, 
              1+i*m+j+1 + iv*m_Ni) );
        //faces.push_back(
        //    TrigBD<int>( 1+(i-1)*m+j,1+(i-1)*m+j+1,1+i*m+j+1) );
        m_faces.push_back(
            TrigBD<int>( 1+(i-1)*m+j + iv*m_Ni ,1+i*m+j+1 + iv*m_Ni,
              1+(i-1)*m+j+1 + iv*m_Ni) );
      }
      int j = m-1;
      m_faces.push_back(
          TrigBD<int>( 1+(i-1)*m+j + iv*m_Ni, 1+i*m+j + iv*m_Ni,
            i*m+(j+2)%m + iv*m_Ni) );
      //faces.push_back(
      //    TrigBD<int>( 1+(i-1)*m+j, (i-1)*m+(j+2)%m, i*m+(j+2)%m ) );
      m_faces.push_back(
          TrigBD<int>( 1+(i-1)*m+j + iv*m_Ni, i*m+(j+2)%m + iv*m_Ni, 
            (i-1)*m+(j+2)%m + iv*m_Ni) );
    }

#ifdef VERBOSE
    for(std::vector< TrigBD<int> >::iterator 
    it = m_faces.begin() ; it != m_faces.end(); ++it)
      std::cout <<"splitted tri-face: "<< it->ai <<"-"<< it->bi <<
        "-"<< it->ci <<std::endl;
    //m_faces.clear();
#endif
  }
  std::cout<<"faces: "<<m_faces.size()<<"\n";
  // end of fill out m_edges and m_faces

  // begin to fill m_faces of boundary
  for(int ip=0; ip<m_np;++ip)
  {
      int off_set = m_N + m_Nip*ip;
      for(int i=0; i<m_pp-1; ++i)
          for(int j=0; j<m_pp-1; ++j)
          {
              int vindex = off_set + m_pp*i + j;
              m_faces.push_back(TrigBD<int>(vindex, vindex+1, vindex+m_pp+1));
              m_faces.push_back(TrigBD<int>(vindex, vindex+m_pp+1, vindex+m_pp));
          }

  }
  std::cout<<"faces: "<<m_faces.size()<<"\n";
  // end of fill m_faces of boundary

  delete[] index;
}

/*
int main()
{
  ContactInterface ci;
  
  std::cout<<"haha1"<<std::endl;
  
  std::string filename = "shape.txt";
  std::ifstream in(filename,std::ios::in);
  std::cout<<"haha2"<<std::endl;
  double number;
  std::vector<double> pos;
  while (in >> number)
  {
    pos.push_back(number);
  }
  in.close();

  std::cout<<"haha3"<<std::endl;

  std::ofstream myfile;
  myfile.open("shape2dumbbells.txt");
  for(int i=0;i<pos.size();++i)
    myfile << std::scientific << std::setprecision(16) << pos[i] << std::endl;

  
  
  for(int i=0;i<pos.size();++i)
  {
    if(i<pos.size()/3)
      myfile << std::scientific << std::setprecision(16) << pos[i]+2 << std::endl;
    else
      myfile << std::scientific << std::setprecision(16) << pos[i] << std::endl;
  }
  
  myfile.close();

  
  std::string filename2 = "shape2dumbbells.txt";
  std::ifstream in2(filename2,std::ios::in);
  double number2;
  std::vector<double> pos2;
  while (in2 >> number2)
  {
    pos2.push_back(number2);
  }
  in2.close();

  std::cout<<"haha4"<<std::endl;
  ci.generateMesh(pos2,8,2);
  std::cout<<"haha5"<<std::endl;
  ci.writeOFF();
  ci.init(SURF_SUBDIVISION);
  
  //ci.loadData("shape.txt");
  //ci.p = 8;
  //ci.nv = 1;
  //cout<<"order: "<<ci.p<<endl;
  //cout<<"nv: "<<ci.nv<<endl;
  //cout<<"np: "<<ci.pos.size()/3<<endl;

  //ci.generateMesh();

  //cout<<"number of edges: "<<ci.edges.size()<<endl;
  //cout<<"number of faces: "<<ci.faces.size()<<endl;

  //ci.writeOFF();

  return 0;
}
*/

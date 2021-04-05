#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include "MeshManager.h"
#include "Structures.h"

//using namespace std;

template<typename T>
class Trig
{
  public:
    
    Trig(T a, T b, T c):ai(a),bi(b),ci(c){};
    ~Trig(){};
    
    T ai,bi,ci;
};

class ContactInterface
{
  public:
    ContactInterface();
    ~ContactInterface();


    double *m_pos;
    int m_p;
    int m_nv;
    int m_N;
    int m_Ni;
  
    Polyhedron *m_P;
      
    std::vector< std::pair<int,int> > m_edges;
    std::vector< Trig<int> > m_faces;

    MeshManager* mesh_manager;

    void writeOFF();
    void generateMesh(std::vector<double> pos, double *pos_pole, int p, int nv);
    void init(SURF_TYPE);
    void getVolumeAndGradient(std::vector<double> &IV, std::vector<double> &IS_BDV, int &IV_size, std::vector<double> &Grad, 
            std::vector<int> &Grad_index, std::vector<double> &pos_s, std::vector<double> &pos_e, 
            double *pos_s_pole, double *pos_e_pole, double minSep, double periodicLength = 0);
    void fillMPos(std::vector<double> &pos, double *pos_pole);

  private:
};

template<typename HDS>
class Build_surface : public CGAL::Modifier_base<HDS> {
public:
    Build_surface(ContactInterface *ci):ci_(ci) {}
    void operator()( HDS& hds) {
        // Postcondition: hds is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        B.begin_surface(ci_->m_N, ci_->m_faces.size(), ci_->m_faces.size()*3);
        
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;
        int nv = ci_->m_nv;
        int Ni = ci_->m_Ni;
        const double *pos = ci_->m_pos;
        for(int iv=0; iv<nv; ++iv)
            for(int i=0; i<Ni; ++i)
                B.add_vertex(Point(pos[i+Ni*0+3*iv*Ni], pos[i+Ni*1+3*iv*Ni], pos[i+Ni*2+3*iv*Ni]));
        for(int i=0; i<ci_->m_faces.size(); i++)
        {
            B.begin_facet();
            B.add_vertex_to_facet( ci_->m_faces[i].ai);
            B.add_vertex_to_facet( ci_->m_faces[i].bi);
            B.add_vertex_to_facet( ci_->m_faces[i].ci);
            B.end_facet();
        }
        
        B.end_surface();
    }
    ContactInterface *ci_;
};

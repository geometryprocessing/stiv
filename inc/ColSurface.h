// Surface.h
//

#ifndef SURFACE_H_
#define SURFACE_H_

#include <string>
#include <set>
#include <map>

#include "Common.h"
#include "Octree.h"
#include "Collision.h"
#include "IntersectionVolume.h"
#include <Eigen/Core>
#include <PQP.h>
#include "Hash.h"

typedef Intersection::AABB<double, Eigen::Vector3d> AABB;
typedef Intersection::Hash MyHash;
typedef CGAL::Box_intersection_d::Box_with_handle_d<
            double, 3, Facet_const_handle> Box;
typedef CGAL::Bbox_3 Bbox;
typedef std::vector<Triangle_3> Triangles_IS;

class ColSurface
{
public:
        ColSurface(Polyhedron *p);
        ~ColSurface();
	
        void controlMeshUpdateIteration(set<unsigned int> &selectedVertices, map<unsigned int, unsigned int> &colliding, 
            bool handleCollisions, bool handleElasticity, Vector_3 &translate, bool rigid, bool moveAll, int limiter, 
            Eigen::SparseVector<double> &gradient, std::vector<double> &IV, Eigen::SparseVector<int> &gradient_index, int &IV_size,
            std::vector<double> &IS_BDV);

	void controlMeshUpdate(std::set<unsigned int> &selectedVertices, map<unsigned int, unsigned int> &colliding, 
            bool handleCollisions, bool handleElasticity, Vector_3 &translate, bool rigid, bool moveAll, int limiter, 
            Eigen::SparseVector<double> &gradient, std::vector<double> &IV, Eigen::SparseVector<int> &gradient_index, int &IV_size,
            std::vector<double> &IS_BDV);

        Polyhedron* getControlMesh() const
        { return m_controlMesh; }
    
        Polyhedron* getSlaveMesh() const
        { return m_slaveMesh; }
	
	AABB* getAABB()
	{
		return m_aabb;
	}
	
	double getEdgeLength()
	{
		return m_edgeLength;
	}
	
	vector<set<unsigned int> > m_subMeshVertices;
	std::vector<std::set<unsigned int> > m_subMeshes;

        bool intersection(const Polyhedron &P, Triangles_IS &triangles_is);
        double minDistance();
        double minDistanceBD(bool flag_orig);
        void loadModel(std::set<unsigned int> subMesh, PQP_Model *b, bool flag_orig=false);
        bool subMeshCheck(int i, int j, bool flag_orig);
        void setThickness(double thickness)
        {
            m_thickness = thickness;
        }
        
        void setNvResident(int nv_resident)
        {
            m_nv_resident = nv_resident;
        }

        void setPeriodicLength(double periodicLength)
        {
            m_periodicLength = periodicLength;
        }

        MyHash hash_grid;

        int m_nv; int m_p;
        int m_np; int m_pp;
        int m_nv_resident;
protected:

        virtual void initialize(bool precRigidGradients = true);

	virtual bool getContinuousIntersections(map<unsigned int, unsigned int> &colliding, vector<Collisions> &collisions);

	void groupCollisions(Collisions &inputColls, map<unsigned int, unsigned int> &colliding, vector<Collisions> &collisions);

	void precomputeRigidGradients();
	
	static void computeNormals(Polyhedron* s);

        Polyhedron *m_controlMesh;
        Polyhedron *m_slaveMesh;
	AABB       *m_aabb;
	
	double     m_edgeLength;
	double	m_thickness;
        double m_periodicLength;

	bool m_rigid;
	bool m_moveAll;

        struct Intersect_facets {
            Triangles_IS &triangles_is;
            Intersect_facets(Triangles_IS &triangles_is)
              : triangles_is(triangles_is)
            {}

	    void operator()( const Box* b, const Box* c) const {
		Halfedge_const_handle h = b->handle()->halfedge();
		// check for shared egde --> no intersection
		if ( h->opposite()->facet() == c->handle()
		     || h->next()->opposite()->facet() == c->handle()
		     || h->next()->next()->opposite()->facet() == c->handle())
		    return;
		// check for shared vertex --> maybe intersection, maybe not
		Halfedge_const_handle g = c->handle()->halfedge();
		Halfedge_const_handle v;
		if ( h->vertex() == g->vertex())
		    v = g;
		if ( h->vertex() == g->next()->vertex())
		    v = g->next();
		if ( h->vertex() == g->next()->next()->vertex())
		    v = g->next()->next();
		if ( v == Halfedge_const_handle()) {
		    h = h->next();
		    if ( h->vertex() == g->vertex())
			v = g;
		    if ( h->vertex() == g->next()->vertex())
			v = g->next();
		    if ( h->vertex() == g->next()->next()->vertex())
			v = g->next()->next();
		    if ( v == Halfedge_const_handle()) {
			h = h->next();
			if ( h->vertex() == g->vertex())
			    v = g;
			if ( h->vertex() == g->next()->vertex())
			    v = g->next();
			if ( h->vertex() == g->next()->next()->vertex())
			    v = g->next()->next();
		    }
		}
		if ( v != Halfedge_const_handle()) {
		    // found shared vertex:
		    CGAL_assertion( h->vertex() == v->vertex());
		    // geomtric check if the opposite segments intersect the triangles
		    Triangle_3 t1( h->vertex()->point(),
				 h->next()->vertex()->point(),
				 h->next()->next()->vertex()->point());
		    Triangle_3 t2( v->vertex()->point(),
				 v->next()->vertex()->point(),
				 v->next()->next()->vertex()->point());
		    Segment_3  s1( h->next()->vertex()->point(),
				 h->next()->next()->vertex()->point());
		    Segment_3  s2( v->next()->vertex()->point(),
				 v->next()->next()->vertex()->point());
		    if ( CGAL::do_intersect( t1, s2)) {
			//cerr << "Triangle_3s intersect (t1,s2):\n    T1: " << t1
			//     << "\n    T2 :" << t2 << endl;
			triangles_is.push_back(t1);
			triangles_is.push_back(t2);
		    } else if ( CGAL::do_intersect( t2, s1)) {
			//cerr << "Triangle_3s intersect (t2,s1):\n    T1: " << t1
			//     << "\n    T2 :" << t2 << endl;
			triangles_is.push_back(t1);
			triangles_is.push_back(t2);
		    }
		    return;
		}
		// check for geometric intersection
		Triangle_3 t1( h->vertex()->point(),
			     h->next()->vertex()->point(),
			     h->next()->next()->vertex()->point());
		Triangle_3 t2( g->vertex()->point(),
			     g->next()->vertex()->point(),
			     g->next()->next()->vertex()->point());
		if ( CGAL::do_intersect( t1, t2)) {
		    //cerr << "Triangle_3s intersect:\n    T1: " << t1 << "\n    T2 :"
		    //     << t2 << endl;
		    triangles_is.push_back(t1);
		    triangles_is.push_back(t2);
		}
	    }
	};

private:

};

#endif


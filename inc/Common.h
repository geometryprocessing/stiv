#pragma once

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

#define DETECTEE // Enable/disable edge detection in the hash grid
//#define DEBUGVG

#include <CGAL/Aff_transformation_3.h>
#include "Structures.h"
#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <limits.h>

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

using namespace std;

typedef ExtendedMesh                   Polyhedron;

typedef Polyhedron::Facet              Facet;

typedef Polyhedron::Vertex             Vertex;

typedef Polyhedron::Halfedge           Halfedge;


/* Iterators */

typedef Polyhedron::Vertex_iterator    Vertex_iterator;

typedef Polyhedron::Vertex_const_iterator
Vertex_const_iterator;

typedef Polyhedron::Halfedge_iterator  Halfedge_iterator;

typedef Polyhedron::Halfedge_const_iterator
Halfedge_const_iterator;

typedef Polyhedron::Facet_iterator     Facet_iterator;

typedef Polyhedron::Facet_const_iterator
Facet_const_iterator;

typedef Polyhedron::Edge_iterator      Edge_iterator;
typedef Polyhedron::Edge_const_iterator Edge_const_iterator;


/* Handles */

typedef Polyhedron::Vertex_handle          Vertex_handle;
typedef Polyhedron::Vertex_const_handle    Vertex_const_handle;

typedef Polyhedron::Halfedge_handle        Halfedge_handle;
typedef Polyhedron::Halfedge_const_handle  Halfedge_const_handle;

typedef Polyhedron::Facet_handle           Facet_handle;
typedef Polyhedron::Facet_const_handle     Facet_const_handle;


/* Circulators */

typedef Facet::Halfedge_around_facet_circulator
Halfedge_Facet_circulator;

typedef Facet::Halfedge_around_facet_const_circulator
Halfedge_Facet_const_circulator;

typedef Vertex::Halfedge_around_vertex_circulator 
Halfedge_Vertex_circulator;

typedef Vertex::Halfedge_around_vertex_const_circulator 
Halfedge_Vertex_const_circulator;


/* Misc */

typedef Polyhedron::size_type          size_type;

typedef CGAL::Color                    Color;

typedef CGAL::Aff_transformation_3<CGALKernel>
Transformation;

typedef CGAL::Triangle_3<CGALKernel>       Triangle_3;
typedef CGAL::Plane_3<CGALKernel>       Plane_3;
typedef CGAL::Segment_3<CGALKernel>     Segment_3;

typedef Eigen::Matrix<int,3,1> Point3i;
typedef Eigen::Matrix<double,3,1> Point3d;

template<class T>
class AABBT
{
public:
	
	AABBT(const T& min, const T& max)
	{
		m_min = min;
		m_max = max;
	}
	
	AABBT()
	{
		m_min = T(0,0,0);
		m_max = T(0,0,0);
	}
	
	~AABBT()
	{
	}
	
	
	T& getMin()
	{
		return m_min;
	}
	
	T& getMax()
	{
		return m_max;
	}
	
	friend std::ostream& operator<<(std::ostream& os, T& p)
	{
		cerr << "AABB: " << &p;
		os << "min: " << p.m_min << " ";
		os << "max: " << p.m_max << endl;
		return os;
	}
	
	T m_min; 
	T m_max; 
};

typedef AABBT<Point3i> AABBi;
typedef AABBT<Point3d> AABBd;

class Utils
{
public:
	static Vector_3 p2v(const Point_3& p)
	{
		return Vector_3(p.x(),p.y(),p.z());
	}

	static Point3d p2e(const Point_3& p)
	{
		return Point3d(p.x(),p.y(),p.z());
	}
	
	static Point_3 v2p(const Vector_3& v)
	{
		return Point_3(v.x(),v.y(),v.z());
	}

	static Vector_3 e2v(const Eigen::Vector3d& p)
	{
		return Vector_3(p.x(),p.y(),p.z());
	}

	static Eigen::Vector3d v2e(const Point_3& p)
	{
		return Eigen::Vector3d(p.x(),p.y(),p.z());
	}
	
	static Point_3 e2p(const Eigen::Vector3d& p)
	{
		return Point_3(p.x(),p.y(),p.z());
	}
	
	static Eigen::Vector3d v2e(const Vector_3& p)
	{
		return Eigen::Vector3d(p.x(),p.y(),p.z());
	}
	
	
	static Point_3 centroid(Triangle_3& t)
	{
		return v2p( 
				   (p2v(t.vertex(0)) + p2v(t.vertex(1)) + p2v(t.vertex(2))) / 3.0
				   );
	}
	
	static vector<Vertex_handle> get1RingV(Vertex_handle& v)
	{
		vector<Vertex_handle> ret;
		Halfedge_Vertex_circulator hvcit = v->vertex_begin();
		do 
		{
			Vertex_handle neighbor_v = hvcit->opposite()->vertex();
			ret.push_back(neighbor_v);
			hvcit++;
		} while (hvcit != v->vertex_begin());
		return ret;
	}

	static vector<Facet_handle> get1RingF(Vertex_handle& v)
	{
		vector<Facet_handle> ret;
		Halfedge_Vertex_circulator hvcit = v->vertex_begin();
		do 
		{
			if (!hvcit->is_border())
				ret.push_back(hvcit->facet());
			hvcit++;
		} while (hvcit != v->vertex_begin());
		return ret;
	}

	static vector<Facet*> get1RingFp(Vertex*& v)
	{
		vector<Facet*> ret;
		Halfedge_Vertex_circulator hvcit = v->vertex_begin();
		do 
		{
			if (!hvcit->is_border())
				ret.push_back(&*hvcit->facet());
			hvcit++;
		} while (hvcit != v->vertex_begin());
		return ret;
	}
	
	static std::vector<Vertex_const_handle> getVertices(const Facet_const_handle& t)
	{
		vector<Vertex_const_handle> vs;
		vs.reserve(3);
		Halfedge_Facet_const_circulator hfc2 = t->facet_begin();
		vs.push_back(hfc2->vertex());
		++hfc2;
		vs.push_back(hfc2->vertex());
		++hfc2;
		vs.push_back(hfc2->vertex());
		return vs;
	}	
	

};

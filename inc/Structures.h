#pragma once

#include <string>
#include <limits>
#include <vector>

//#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/Color.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

//#include "StructuresAPI.h"
//#include "Options.h"
#include <Eigen/Core>

#include <CGAL/box_intersection_d.h>

const double EPSILON      = std::numeric_limits<double>::min();
const int    LARGE_INT    = std::numeric_limits<int>::max();
const double LARGE_DOUBLE = std::numeric_limits<double>::max();

const CGAL::Color default_vertex_color = CGAL::Color(85, 156, 200, 255);
const CGAL::Color colliding_vertex_color = CGAL::Color(205, 92, 92, 255);
const CGAL::Color default_control_vertex_color = CGAL::Color(255,255,128, 255);
const CGAL::Color selected_vertex_color = CGAL::Color(72,61,139, 255);

class WeightInfo
{
public:
	int vis[4];
	
	bool is_border;
	bool opposite_is_border;
	
	double w;
	
};

//typedef CGAL::Simple_cartesian<double> CGALKernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel  CGALKernel;
//typedef CGAL::Cartesian<double> CGALKernel;

typedef CGAL::Point_3<CGALKernel>          Point_3;
typedef CGAL::Point_2<CGALKernel>          Point_2;

typedef CGAL::Vector_3<CGALKernel>         Vector_3;
typedef CGAL::Vector_2<CGALKernel>         Vector_2;

typedef CGAL::Polygon_2<CGALKernel>        Polygon_2;

class IndexedItem {
    public:
        IndexedItem() : index(-1) {}
        virtual ~IndexedItem() {}
    public:
        int getIndex() const { return index; }
        void setIndex(int v) { index = v; }
    private:
        int index;
};

///@TODO: Do we need to export this class with STRUCTURES_API ?
template <typename Refs, typename T, typename Point>
class /*STRUCTURES_API*/ MeshVertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, Point> , public IndexedItem 
{
public:
	MeshVertex() 
	{
		init();
	}
	
	MeshVertex(const Point& p) : CGAL::HalfedgeDS_vertex_base<Refs, T, Point>(p) 
	{
		init();
	}
	
	void reset() { init(); }
	
protected:
	void init() 
	{
          /*
	       CGAL::Color c(Options::getIntValue("color_r", 85),
				      Options::getIntValue("color_g", 156),
				      Options::getIntValue("color_b", 200),
					  Options::getIntValue("color_a", 255));
          */
                CGAL::Color c(85,156,200,255);
		color = c;
		default_color = c;
		normal  = Vector_3(0.0, 0.0, 0.0);
		//new_pos = Point_3(0.0, 0.0, 0.0);
	
		laplacian = Vector_3(0.0, 0.0, 0.0);
		marked  = true;
		handle = 0;
		handleCanMove = false;
		skip = false;
		substep = false;
		
		hide = false;
	}
	
public:
        int meshTag;
	CGAL::Color color;
	CGAL::Color default_color;
	Vector_3    normal;
//	Point_3     ori_pos;
//	Point_3     new_pos;
	bool        marked;
	short       handle;
	Vector_3    laplacian;
	bool handleCanMove;
	bool skip;
	bool substep;
	
	double min[3];
	double max[3];
	
	bool hide;
	
	std::vector<int> ns;
	std::vector<WeightInfo> nsW;
	
	Eigen::Vector3d m_pos;
	Eigen::Vector3d m_oripos;
	Eigen::Vector3d target;
	Eigen::Vector3d origin;
	
	Eigen::Vector3d& P()
	{
		return m_pos;
	}

	Eigen::Vector3d& Pori()
	{
		return m_oripos;
	}
	
	Point_3 Pcgal() const
	{
		return Point_3(m_pos[0],m_pos[1],m_pos[2]);
	}
};


template <typename Refs>
class /*STRUCTURES_API*/ MeshHalfedge : public CGAL::HalfedgeDS_halfedge_base<Refs>, public IndexedItem 
{
public:
	MeshHalfedge() 
	{
		color  = CGAL::BLACK;
		length = 0;
		weight = 0;
		marked = true;
	}
	
public:
	CGAL::Color color;
	CGALKernel::FT  length;
	double weight;
	bool        marked;
	
	double min[3];
	double max[3];
};

template <typename Refs, typename T, typename Plane>
class /*STRUCTURES_API*/ MeshFace : public CGAL::HalfedgeDS_face_base<Refs, T, Plane>, public IndexedItem {
public:
	MeshFace() {
		init();
	}
	
public:
	const Vector_3& getUnitNormal() const {
		return normal;
	}
	
	void setNormal(Vector_3 n) {
		normal = n;
	}
	
	void reset() { init(); }
	
	void init() {
		normal = Vector_3(0.0, 0.0, 0.0);
		marked = true;
		ignore = false;
		hide = false;
	}
	
public:
	Vector_3    normal;
	bool        marked; // marked faces have been moved in the current time step
	int vi[3];
	bool ignore;
	
	double min[3];
	double max[3];
	
	bool hide;
	
};

struct /*STRUCTURES_API*/ MeshItems : public CGAL::Polyhedron_items_3 {
    template < typename Refs, typename Traits>
    struct Vertex_wrapper {
        typedef typename Traits::Point_3 Point;
        typedef MeshVertex<Refs, CGAL::Tag_true, Point> Vertex;
    };
    template < typename Refs, typename Traits>
    struct Halfedge_wrapper {
        typedef MeshHalfedge<Refs> Halfedge;
    };
    template < typename Refs, typename Traits>
    struct Face_wrapper {
        typedef typename Traits::Plane_3 Plane;
        typedef MeshFace<Refs, CGAL::Tag_true, Plane> Face;
    };
};


typedef CGAL::Polyhedron_3<CGALKernel, MeshItems> ListBasedMesh;
typedef CGAL::Polyhedron_3<CGALKernel, MeshItems, CGAL::HalfedgeDS_vector>
                                              VectorBasedMesh;
class /*STRUCTURES_API*/ ExtendedMesh : public VectorBasedMesh {
    public:
        ExtendedMesh();

        /**
         * Reset all flags
         */
        void reset();

        /**
         * Initialized vertex/facet indices
         */
        void init();
	void updateWeight();

	    void initAABB(double h);

    public:
        /**
         * Member variable get/set methods.
         */
	bool getVisible() const { return visible; }
	void setVisible(bool v) { visible = v; }
	bool getActive() const { return active; }
	void setActive(bool v) { active = v; }
	std::string getName() const { return name; }
	void setName(const std::string& v) { name = v; }
	std::string getFilename() const { return filename; }
	void setFilename(const std::string& v) { filename = v; }
	Vertex* getVertex(int i) const { return v_array.at(i); }
	Facet*  getFacet (int i) const { return f_array.at(i); }
	
	double getArea(unsigned int idx);
	std::vector<std::pair<int, double> > get1Ring(int i);
	std::vector<std::pair<int, Halfedge*> > get1RingH(int i);
	
	std::vector<Vertex*> v_array;
	std::vector<Halfedge*> e_array; // one per edge, so half of the he are contained in this vector
	std::vector<Facet*>  f_array;
	
	void updateAllCotanPori();
	
	double computeWeight(Halfedge* h);
	double uniformWeight(Halfedge* h) const ;
	double cotangentWeight(Halfedge* h) const ;
	double cotangentWeightPori(Halfedge* h) const ;
	/// Return cotangent of (P,Q,R) corner (i.e. cotan of QP,QR angle).
	double cotangent(const Eigen::Vector3d& P,
					 const Eigen::Vector3d& Q,
					 const Eigen::Vector3d& R) const;
	
	WeightInfo getWeightVs(Halfedge* h);
	double computeWeightVs(std::vector<Eigen::Vector3d>& ps, WeightInfo& wi);
	
	void copyCgalToEigen();
	void copyEigenToCgal();
	
	Eigen::VectorXd& getAllPos()
	{
		return m_pos;
	}
	
	Eigen::VectorXd& getAllOriPos()
	{
		return m_oripos;
	}
	
    private:
        bool visible;
        bool active;
        std::string name;
        std::string filename; // including path
	
	Eigen::VectorXd m_pos;
	Eigen::VectorXd m_oripos;
};


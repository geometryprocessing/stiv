#ifndef PRISMSURFACE_H_
#define PRISMSURFACE_H_

#include "Common.h"
#include "ColSurface.h"
#include "Collision.h"

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/intersections.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include "Octree.h"

#include <Eigen/Core>


using namespace std;

typedef Intersection::AABB<double, Eigen::Vector3d> AABB;
typedef CGAL::Bbox_3                                     Bbox;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel  CGALKernel;
typedef CGALKernel::Triangle_3                               Triangle;
typedef CGALKernel::Segment_3                                Segment;
typedef CGAL::Box_intersection_d::Box_with_handle_d<
double, 3, Facet_const_handle>               Box;

using Eigen::Vector3d;

class PrismVolume
{
public:
			
	static bool validateCollision(Vertex*& vh, Vector3d& before0, Vector3d& before1, Vector3d& before2,
											   Vector3d& after0, Vector3d& after1, Vector3d& after2,
								  double h, double& tret);
	static void validateCollisionAndAdd(Facet* fit, Vertex* vit, double h, vector<Collision>& collisions, Polyhedron *p);
	static bool validateCollisionAndReplace(Facet* fit, Vertex* vit, double h, Collision& c, Polyhedron *p, double periodicLength = 0);
	static bool TestSegmentAABB(Vector3d p0, Vector3d p1, AABB& b);
	static bool TestSegmentAABB(Vector3d &p0, Vector3d &p1, Vector3d &min, Vector3d &max);
	static bool validateCollisionFastFast(Facet* fit, Vertex* vit);
	static bool validateCollisionFast(Facet* fit, Vertex* vit, Polyhedron *p, double periodicLength = 0);
	
        static void getTranslationVector(Vector3d &transV, Vertex* vit1, Vertex* vit2, double periodicLength);
	
	static bool validateCollisionEEFastFast(Halfedge* eit1, Halfedge* eit2);
	static bool validateCollisionEEFast(Halfedge* eit1, Halfedge* eit2, Polyhedron *p, double periodicLength = 0, double h = 0);
	static bool validateCollisionEEAndReplace(Halfedge* eit1, Halfedge* eit2, double h, Collision& c, Polyhedron *p, double periodicLength=0);

	static bool validateCollision(Eigen::Vector3d &s11, Eigen::Vector3d &s12, Eigen::Vector3d &s21, Eigen::Vector3d &s22,
								  Eigen::Vector3d &e11, Eigen::Vector3d &e12, Eigen::Vector3d &e21, Eigen::Vector3d &e22, double h, double& tret);

	static void computeVolumeWithHash(Polyhedron& m, AABB& aabb, double thickness, double edgeLength, std::vector<std::set<unsigned int> > &subMeshes, vector<Collision> &collisions, vector<Collision> &collisionsEE, MyHash &h, int nv_resident, double periodicLength = 0);

	static bool EdgeEdgeIntersection(Eigen::Vector3d &x0,Eigen::Vector3d &x1,Eigen::Vector3d &x2,Eigen::Vector3d &x3,
							  Eigen::Vector3d &v0,Eigen::Vector3d &v1,Eigen::Vector3d &v2,Eigen::Vector3d &v3, double &tRet);

	static bool VertexFaceIntersection(Eigen::Vector3d &x0,Eigen::Vector3d &x1,Eigen::Vector3d &x2,Eigen::Vector3d &x3,
									   Eigen::Vector3d &v0,Eigen::Vector3d &v1,Eigen::Vector3d &v2,Eigen::Vector3d &v3, double &tRet);
	
	static bool VertexFaceProximity(Eigen::Vector3d& x0, Eigen::Vector3d& x1, Eigen::Vector3d& x2, Eigen::Vector3d& x3,
									Eigen::Vector3d& v0, Eigen::Vector3d& v1, Eigen::Vector3d& v2, Eigen::Vector3d& v3,
									double h, double& tRet);
       
        static bool EdgeEdgeProximity(Eigen::Vector3d& x0, Eigen::Vector3d& x1, Eigen::Vector3d& x2, Eigen::Vector3d& x3,
									Eigen::Vector3d& v0, Eigen::Vector3d& v1, Eigen::Vector3d& v2, Eigen::Vector3d& v3,
									double h, double& tRet);


	static void getBarycentricCoordinates(Vector3d &x0, Vector3d &x1, Vector3d &x2, Vector3d &x3, bool vt, double w[4]);
	
	static void getTimeGradient(Vector3d &x1, Vector3d &x2, Vector3d &x3, Vector3d &x4,
								Vector3d &v1, Vector3d &v2, Vector3d &v3, Vector3d &v4, double w[4], double t, double h,
								Vector3d tGrad[4]);

	static void getEETimeGradient(Vector3d &x1, Vector3d &x2, Vector3d &x3, Vector3d &x4,
								  Vector3d &v1, Vector3d &v2, Vector3d &v3, Vector3d &v4, double w[4], double t, double h,
								  Vector3d tGrad[4]);
	
	static void getBarycentricGradients(Vector3d &x0, Vector3d &x1, Vector3d &x2, Vector3d &x3, Vector3d wGrad[4][4]);
	static void getEEBarycentricGradients(Vector3d &x0, Vector3d &x1, Vector3d &x2, Vector3d &x3, Vector3d wGrad[4][4]);

	static int areCoplanar(Eigen::Vector3d& x0, Eigen::Vector3d& x1, Eigen::Vector3d& x2, Eigen::Vector3d& x3,
									   Eigen::Vector3d& x0v, Eigen::Vector3d& x1v, Eigen::Vector3d& x2v, Eigen::Vector3d& x3v, double *dt);

	static int areProximate(Eigen::Vector3d& x0, Eigen::Vector3d& x1, Eigen::Vector3d& x2, Eigen::Vector3d& x3,
							Eigen::Vector3d& v0, Eigen::Vector3d& v1, Eigen::Vector3d& v2, Eigen::Vector3d& v3, double h, double *dt);
        
        static int areEEProximate(Eigen::Vector3d& x0, Eigen::Vector3d& x1, Eigen::Vector3d& x2, Eigen::Vector3d& x3,
							Eigen::Vector3d& v0, Eigen::Vector3d& v1, Eigen::Vector3d& v2, Eigen::Vector3d& v3, double h, double *dt);

	
	static double vertexFaceDist(Eigen::Vector3d& p, Eigen::Vector3d& a, Eigen::Vector3d& b, Eigen::Vector3d& c, double &t1, double &t2, double &t3);
	
	static double edgeEdgeDist(Vector3d &p1, Vector3d &q1, Vector3d &p2, Vector3d &q2, double &s, double &t);
        static int areLineLineInt(Vector3d& x0, Vector3d& x1, Vector3d& x2, Vector3d& x3,
            Vector3d& v0, Vector3d& v1, Vector3d& v2, Vector3d& v3, double h, double *dt);
        static int arePointPointInt(Vector3d& x0, Vector3d& x1,
            Vector3d& v0, Vector3d& v1, double h, double *dt);
        static int arePointLineInt(Vector3d& x0, Vector3d& x1, Vector3d& x2,
            Vector3d& v0, Vector3d& v1, Vector3d& v2, double h, double *dt);
        static int areEdgeEdgeInt(Vector3d& x0, Vector3d& x1, Vector3d& x2, Vector3d& x3,
            Vector3d& v0, Vector3d& v1, Vector3d& v2, Vector3d& v3, double h, double *dt);
	static bool TestSegmentSubdivideSpace(Vector3d &e1v1s, Vector3d &e1v1e, Vector3d &e1v2s, Vector3d &e1v2e,
                                              Vector3d &e2v1s, Vector3d &e2v1e, Vector3d &e2v2s, Vector3d &e2v2e, double h);
        static bool TestSegmentSubdivideTime(Vector3d &e1v1s, Vector3d &e1v1e, Vector3d &e1v2s, Vector3d &e1v2e,
                                             Vector3d &e2v1s, Vector3d &e2v1e, Vector3d &e2v2s, Vector3d &e2v2e, double h);
        static bool TestSegmentSegment(Vector3d &e1v1s, Vector3d &e1v1e, Vector3d &e1v2s, Vector3d &e1v2e,
                                       Vector3d &e2v1s, Vector3d &e2v1e, Vector3d &e2v2s, Vector3d &e2v2e, double h);
        static void SegmentAABB(Vector3d &v1s, Vector3d &v1e, Vector3d &v2s, Vector3d &v2e,
                                       Vector3d &emin, Vector3d &emax, double h);
};

#endif

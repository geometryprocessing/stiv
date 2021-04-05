#include "Octree.h"

#include "Common.h"
#include "ColSurface.h"

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/intersections.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <algorithm>
#include <vector>

using std::vector;
using namespace std;

namespace Intersection
{
	void test1()
	{
//		typedef Intersection::AABB<float,Vector_3> AABB;
//		typedef Intersection::SphereContainer<float,Vector_3,float*> Sphere;
//		Intersection::Octree<AABB,Sphere> octree(AABB(Vector_3(0.5,0.5,0.5),Vector_3(1,1,1)));
//		float temp = 3;
//		Sphere* s1 = new Sphere(Vector_3(0.25,0.25,0.25),0.2,&temp);
//		octree.add(*s1);
//		
//		list<Sphere*> l = octree.windowQueryAABB(AABB(Vector_3(0.25,0.25,0.25),Vector_3(0.25,0.25,0.25)));
//		assert(l.size() == 1);
//		delete s1;
	}

	void test2()
	{
//		typedef Intersection::AABB<float,Vector_3> AABB;
//		typedef Intersection::SphereContainer<float,Vector_3,float*> Sphere;
//		Intersection::Octree<AABB,Sphere> octree(AABB(Vector_3(0.5,0.5,0.5),Vector_3(1,1,1)));
//		float temp = 3;
//		Sphere* s1 = new Sphere(Vector_3(0.49,0.49,0.49),0.2,&temp);
//		octree.add(*s1);
//		Sphere* s2 = new Sphere(Vector_3(0.25,0.25,0.25),0.2,&temp);
//		octree.add(*s2);
//		
//		list<Sphere*> l = octree.windowQueryAABB(AABB(Vector_3(0.75,0.75,0.75),Vector_3(0.5,0.5,0.5)));
//		assert(l.size() == 1);
//
//		delete s1;
//		delete s2;
//
	}
	
//	bool OctreeTestSuite::runTestSuite()
//	{
//		test1();
//		test2();
//		
//	}
//	
	
	
}

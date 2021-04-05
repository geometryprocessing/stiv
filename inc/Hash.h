#ifndef HASH_H_
#define HASH_H_

#include "Common.h"
#include <CGAL/Polyhedron_3.h>
#include <iostream>
#include <algorithm>
#include <parallel/algorithm>
#include <vector>
#include <list>
#include <iomanip>
#include "Collision.h"
//#include "PrismVolume.h"

#include <Eigen/Core>

#include <boost/foreach.hpp>
#include <ext/hash_map>
#include "Octree.h"
typedef Intersection::AABB<double, Eigen::Vector3d> AABB;
using Eigen::Vector3d;

/**
 * This class implements a Hash-based intersaction system
 */

namespace Intersection
{
    using namespace std;
	
	class HashItem
	{
	public:
		HashItem(long k, int i)
		{
			key = k;
			id = i;
		}

		HashItem(){}

		long key;
		int id;
		
		bool operator<(const HashItem& hi) const
		{
			return key < hi.key;
		}
	};
	
	class Hash
	{
	public:
		
		//Hash(Vector3d min, Vector3d max, double cellSize, double h, double periodicLength = 0);
		Hash();
		
		virtual ~Hash();
		
		
		void resetHash(Vector3d min, Vector3d max, double cellSize, double h, int nv_resident, double periodicLength = 0);
		
                void fill(Polyhedron* p, std::vector<std::set<unsigned int> > &subMeshes);
		AABBi makeAABBi(double min[3], double max[3]);
		
		void addTriangle(Polyhedron* p, Facet* fit, int id, vector<HashItem>& hs);
		void addVertex(Polyhedron* p, Vertex* vit, int id, vector<HashItem>& hs);
		void addEdge(Polyhedron* p, Halfedge* eit, int id, vector<HashItem>& hs);
		
		// CPU SORT
		typedef std::pair<int, int> intPair;
		void sort();
		
		//typedef std::pair<Facet_iterator, Vertex_iterator> FV;
		// THIS MAY USE THE GPU IF NEEDED
		void computeVolume(vector<Collision>& collisions, vector<Collision>& collisionsEE);
		void computeVolumeParallel(vector<Collision>& collisions, vector<Collision>& collisionsEE);
		void add(AABBi aabbi, int id, vector<HashItem>& hs);
		void add(Point3i p, int id, vector<HashItem>& hs);

		long hash(Point3i p);		
		void clear();
		
		double m_cellSize;
		double m_thickness;
                int m_nv_resident;
                double m_periodicLength;

                int m_periodicSize;
		int m_gridSize;
		int m_numberTriangle;
		Vector3d m_domainMin;
		Vector3d m_domainMax;
		
		Polyhedron* m_p;

		HashItem& get(int i);
		HashItem& getEdge(int i);
		
		
	private:
		vector<HashItem> m_hash;
		vector<HashItem> m_hash_edge;

                Facet* facet[1];
                Vertex* vs[1];
                Halfedge* es[1];
	};
}

typedef Intersection::Hash MyHash;

#endif /* HASH_H_ */

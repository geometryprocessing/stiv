#ifndef OCTREE_H_
#define OCTREE_H_

#include "Common.h"

#include <CGAL/Polyhedron_3.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <list>
#include <iomanip>


/**
 * This class implements a Loose Multiple Storage Octree
 * Items are not replicated, the memory management of the items is handled by the octree
 */

namespace Intersection
{
    using namespace std;
	
    typedef enum
    {
		NT_INTERNAL,
		NT_LEAF
    } NodeType;
	
    typedef enum
    {
		OP_OK,
		OP_DOWN,
		OP_UP
    } OctreePosition;
	

    template<class MyScalarType, class MyVectorType>
    class AABB
    {
	public:
		
		typedef MyScalarType ScalarType;
		typedef MyVectorType VectorType;
		typedef AABB<ScalarType,VectorType> AABBType;
		
        // for stl containers only
        AABB()
        {
            m_center = VectorType(0,0,0);
            m_size   = VectorType(0,0,0);
			computeMinMax();
        }
		
		AABB(VectorType Center, VectorType dim)
		{
			m_center = Center;
			m_size   = dim;
			computeMinMax();
		}

		AABB(vector<VectorType>& vv)
		{
			assert(vv.size() > 0);
			double min[3] = {vv[0][0],vv[0][1],vv[0][2]};
			double max[3] = {vv[0][0],vv[0][1],vv[0][2]};
			for(unsigned j=0; j<3; ++j)
			{
				for(unsigned i=1;i<vv.size();++i)
				{
					min[j] = std::min(min[j],vv[i][j]);
					max[j] = std::max(max[j],vv[i][j]);
				}
			}
			initMinMax(VectorType(min[0],min[1],min[2]),VectorType(max[0],max[1],max[2]));
		}

		~AABB()
		{
		}

		void initMinMax(VectorType Min, VectorType Max)
		{
			assert(Min[0] <= Max[0]);
			assert(Min[1] <= Max[1]);
			assert(Min[2] <= Max[2]);
			m_center = (Min+Max)/2.0;
			m_size   = Max - Min;
			m_min = Min;
			m_max = Max;
		}
		
		void computeMinMax()
		{
			m_min = getMinSlow();
			m_max = getMaxSlow();
			m_mine = getMinESlow();
			m_maxe = getMaxESlow();
		}
		
		AABBType getAABB()
		{
			return *this;
		}
		
		VectorType getMinSlow()
		{
			return m_center - (m_size/2.0);
		}
		
		VectorType getMaxSlow()
		{
			return m_center + (m_size/2.0);
		}

		VectorType& getMin()
		{
			return m_min;
		}
		
		VectorType& getMax()
		{
			return m_max;
		}

		VectorType getMinESlow()
		{
			return m_center - ((m_size * 1.5)/2.0);
		}
		
		VectorType getMaxESlow()
		{
			return m_center + ((m_size * 1.5)/2.0);
		}
		
		VectorType& getMinE()
		{
			return m_mine;
		}
		
		VectorType& getMaxE()
		{
			return m_maxe;
		}
		
		
		bool intersect(AABB& c)
		{
			VectorType& cmin = c.getMin();
			VectorType& cmax = c.getMax();
			
			VectorType& min = getMin();
			VectorType& max = getMax();
			
			bool nointersect = false;
			nointersect = nointersect || (cmax[0] < min[0]) || (cmin[0] > max[0]);
			nointersect = nointersect || (cmax[1] < min[1]) || (cmin[1] > max[1]);
			nointersect = nointersect || (cmax[2] < min[2]) || (cmin[2] > max[2]);
			return !nointersect;
		}

		bool extendedIntersect(AABB& c)
		{
			VectorType& cmin = c.getMin();
			VectorType& cmax = c.getMax();
			
			VectorType& min = getMinE();
			VectorType& max = getMaxE();
			
			bool nointersect = false;
			nointersect = nointersect || (cmax[0] < min[0]) || (cmin[0] > max[0]);
			nointersect = nointersect || (cmax[1] < min[1]) || (cmin[1] > max[1]);
			nointersect = nointersect || (cmax[2] < min[2]) || (cmin[2] > max[2]);
			return !nointersect;
		}
		
		bool contains(VectorType& p)
		{
			VectorType& min = getMin();
			VectorType& max = getMax();
			
			bool res = true;
			res = res && min[0] <= p[0];
			res = res && min[1] <= p[1];
			res = res && min[2] <= p[2];
			
			res = res && max[0] > p[0];
			res = res && max[1] > p[1];
			res = res && max[2] > p[2];
			
			return res;
		}
		
		bool contains(AABB<ScalarType,VectorType>& c)
		{
			VectorType& mint = c.getMin();
			VectorType& maxt = c.getMax();
            
			return contains(mint) && contains(maxt);
		}
		
		bool split(vector<AABBType >& v)
		{
			v.resize(8);
			VectorType halfSize = m_size/2.0;
			VectorType shift = m_size/4.0;
			
			
			v[0] = AABBType(m_center + VectorType(-shift[0],-shift[1],-shift[2]),VectorType(halfSize[0],halfSize[1],halfSize[2]));
			v[1] = AABBType(m_center + VectorType( shift[0],-shift[1],-shift[2]),VectorType(halfSize[0],halfSize[1],halfSize[2]));
			v[2] = AABBType(m_center + VectorType(-shift[0], shift[1],-shift[2]),VectorType(halfSize[0],halfSize[1],halfSize[2]));
			v[3] = AABBType(m_center + VectorType( shift[0], shift[1],-shift[2]),VectorType(halfSize[0],halfSize[1],halfSize[2]));
			v[4] = AABBType(m_center + VectorType(-shift[0],-shift[1], shift[2]),VectorType(halfSize[0],halfSize[1],halfSize[2]));
			v[5] = AABBType(m_center + VectorType( shift[0],-shift[1], shift[2]),VectorType(halfSize[0],halfSize[1],halfSize[2]));
			v[6] = AABBType(m_center + VectorType(-shift[0], shift[1], shift[2]),VectorType(halfSize[0],halfSize[1],halfSize[2]));
			v[7] = AABBType(m_center + VectorType( shift[0], shift[1], shift[2]),VectorType(halfSize[0],halfSize[1],halfSize[2]));
			
			return true;
		}
		
		VectorType& getCenter()
		{
			return m_center;
		}
		
		void setCenter(VectorType& v)
		{
			m_center = v;
		}
		
		VectorType& getSize()
		{
			return m_size;
		}
		
		void setSize(VectorType& v)
		{
			m_size = v;
		}
		
		friend std::ostream& operator<<(std::ostream& os, AABBType& p)
		{
			cerr << "AABB: " << &p;
			os << std::setprecision(10) << "center: " << p.m_center << " ";
			os << std::setprecision(10) << "size: " << p.m_size << " ";
			
			return os;
		}
		
		VectorType m_center;
		VectorType m_size; 
		VectorType m_min; 
		VectorType m_max; 
		VectorType m_mine; 
		VectorType m_maxe; 
    };
	
//	
//	template<class AABBType, class ElementType>
//	class Node;
//	
//	template<class AABBType, class ElementType>
//	class LeafNode;
//	
//	template<class AABBType, class ElementType>
//	class InternalNode;
//	
//	template<class MyScalarType, class MyVectorType, class MyContainedType>
//	class SphereContainer
//	{
//	public:
//		typedef MyScalarType ScalarType;
//		typedef MyVectorType VectorType;
//		typedef MyContainedType ContainedType;
//		typedef AABB<ScalarType,VectorType> AABBType;
//		typedef SphereContainer<ScalarType,VectorType, ContainedType> SphereType;
//		typedef LeafNode<AABBType,SphereType> LeafType;
//		
//		SphereContainer(VectorType center, ScalarType radius, ContainedType p)
//		{
//			m_center = center;
//			m_radius = radius;
//			m_obj    = p;
//		}
//
//		SphereContainer()
//		{
//		}
//		
//		virtual ~SphereContainer()
//		{
//		}
//		
//		void clear()
//		{
//			octreeNode = 0;
//			m_obj->octreeElem = 0;
//		}
//		
//		// Computes the square distance between a point p and an AABB b
//		float SqDistPointAABB(Vector_3 p, AABBType& b)
//		{
//			float sqDist = 0.0f;
//			Vector_3 min = b.getMin();
//			Vector_3 max = b.getMax();
//			for (int i = 0; i < 3; i++) 
//			{
//				// For each axis count any excess distance outside box extents
//				float v = p[i]; 
//				if (v < min[i]) sqDist += (min[i] - v) * (min[i] - v);
//				if (v > max[i]) sqDist += (v - max[i]) * (v - max[i]);
//			}
//			return sqDist;		
//		}
//		
//		bool intersect(AABBType& aabb)
//		{
//			float sqDist = SqDistPointAABB(m_center, aabb);
//			// Sphere and AABB intersect if the (squared) distance 
//			// between them is less than the (squared) sphere radius 
//			return sqDist <= m_radius * m_radius;		
//		}
//		
//		bool isContainedIn(AABBType& aabb)
//		{
//			Vector_3 c = aabb.getCenter();
//			Vector_3 s = aabb.getSize();
//			float dR = 2.0*m_radius;
//			s = s - Vector_3(dR,dR,dR);
//			
//			return AABBType(c,s).contains(m_center);
//		}
//		
//		OctreePosition checkPosition(AABBType& aabb)
//		{
//			float objsize = m_radius*2.0;
//			float boxsize = aabb.getSize()[0];
//			assert(aabb.getSize()[0] == aabb.getSize()[1]); // this code works only for boxes!
//			assert(aabb.getSize()[1] == aabb.getSize()[2]); // this code works only for boxes!
//			
//			if (objsize <= 0.5*boxsize)
//			{
//				return OP_DOWN;
//			}
//			else if ((0.5*boxsize < objsize) && (objsize < boxsize))
//			{
//				return OP_OK;
//			}
//			else 
//			{
//				return OP_UP;
//			}
//		}
//		
//		VectorType& getCenter()
//		{
//			return m_center;
//		}
//		
//		void setOctreeNode(LeafType* ln)
//		{
//			octreeNode = ln;
//			m_obj->octreeElem = (void*)this; // warning: do not relocate this object!
//		}
//		
//		ContainedType m_obj;
//		ScalarType m_radius;
//		VectorType m_center;
//		LeafType* octreeNode;
//		
//	};
//
//	template<class MyScalarType, class MyVectorType, class MyContainedType>
//	class SegmentContainer
//	{
//	public:
//		typedef MyScalarType ScalarType;
//		typedef MyVectorType VectorType;
//		typedef MyContainedType ContainedType;
//		typedef AABB<ScalarType,VectorType> AABBType;
//		typedef SegmentContainer<ScalarType,VectorType, ContainedType> SegmentType;
//		typedef LeafNode<AABBType,SegmentType> LeafType;
//		
//		SegmentContainer(VectorType a, VectorType b, ContainedType p)
//		{
//			m_a = a;
//			m_b = b;
//			m_obj    = p;
//		}
//		
//		SegmentContainer()
//		{
//		}
//		
//		virtual ~SegmentContainer()
//		{
//		}
//		
//		void clear()
//		{
//			octreeNode = 0;
//			m_obj->octreeElem = 0;
//		}
//
//		
//		// Test if segment specified by points p0 and p1 intersects AABB b
//		static bool TestSegmentAABB(VectorType p0, VectorType p1, AABBType& b) 
//		{
//			VectorType min = b.getMin();
//			VectorType max = b.getMax();
//			VectorType c = (min + max) * 0.5f; // Box center-point
//			VectorType e = max - c;             // Box halflength extents
//			VectorType m = (p0 + p1) * 0.5f;       // Segment midpoint
//			VectorType d = p1 - m;                // Segment halflength vector
//			m = m - c;                        // Translate box and segment to origin
//			
//			// Try world coordinate axes as separating axes
//			double adx = fabs(d[0]);
//			if (fabs(m[0]) > e[0] + adx) return false;
//			double ady = fabs(d[1]);
//			if (fabs(m[1]) > e[1] + ady) return false;
//			double adz = fabs(d[2]);
//			if (fabs(m[2]) > e[2] + adz) return false;
//			
//			// Add in an epsilon term to counteract arithmetic errors when segment is 
//			// (near) parallel to a coordinate axis (see text for detail) 
//			adx += EPSILON; ady += EPSILON; adz += EPSILON;
//			
//			// Try cross products of segment direction vector with coordinate axes
//			if (fabs(m[1] * d[2] - m[2] * d[1]) > e[1] * adz + e[2] * ady) return false;
//			if (fabs(m[2] * d[0] - m[0] * d[2]) > e[0] * adz + e[2] * adx) return false;
//			if (fabs(m[0] * d[1] - m[1] * d[0]) > e[0] * ady + e[1] * adx) return false;
//			
//			// No separating axis found; segment must be overlapping AABB 
//			return true;
//		}
//		
//		
//		bool intersect(AABBType& aabb)
//		{
//			return TestSegmentAABB(m_a,m_b,aabb);
//			
//		}
//		
//		bool isContainedIn(AABBType& aabb)
//		{
//			return aabb.contains(m_a) && aabb.contains(m_b);
//		}
//		
//		OctreePosition checkPosition(AABBType& aabb)
//		{
//			float objsize = sqrt((m_a-m_b).squared_length());
//			float boxsize = aabb.getSize()[0];
//			assert(aabb.getSize()[0] == aabb.getSize()[1]); // this code works only for boxes!
//			assert(aabb.getSize()[1] == aabb.getSize()[2]); // this code works only for boxes!
//			
//			if (objsize <= 0.5*boxsize)
//			{
//				return OP_DOWN;
//			}
//			else if ((0.5*boxsize < objsize) && (objsize < boxsize))
//			{
//				return OP_OK;
//			}
//			else 
//			{
//				return OP_UP;
//			}
//		}
//		
//		VectorType getCenter()
//		{
//			return (m_a + m_b)/2.0;
//		}
//		
//		void setOctreeNode(LeafType* ln)
//		{
//			octreeNode = ln;
//			m_obj->octreeElem = (void*)this; // warning: do not relocate this object!
//		}
//		
//		ContainedType m_obj;
//		VectorType m_a;
//		VectorType m_b;
//		LeafType* octreeNode;
//		
//	};
//	
//	template<class AABBType, class ElementType>
//	class Octree
//	{
//	public:
//		
//		typedef typename AABBType::ScalarType ScalarType;
//		typedef typename AABBType::VectorType VectorType;
//		typedef LeafNode<AABBType,ElementType> LeafType;
//		typedef InternalNode<AABBType,ElementType> InternalType;
//		typedef Node<AABBType,ElementType> NodeType;
//
//		Octree(AABBType c)
//		{
//			root = new LeafNode<AABBType,ElementType>(c);
//			maxlevel = 0;
//			levelLimit = 10;
//		}
//		
//		virtual ~Octree(){};
//		
//		void add(ElementType l)
//		{
////			assert(l.isContainedIn(root->getAABB()));
//			ElementType* lp = new ElementType();
//			*lp = l;
//			
//			add_aux(&root,lp,0);
//		}
//		
//		void remove(ElementType* l)
//		{
//            remove_aux(&root,l);
//		}
//		
//		list<ElementType*> windowQueryAABB(AABBType c)
//		{
//            list<ElementType*> labels;
//            windowQuery_auxAABB(root,c, labels);
//            return labels;
//		}
//		
//		AABBType getDomain()
//		{
//            return root->getAABB();
//		}
//		
//	private:
//		unsigned int maxlevel;
//		
//		void add_aux(NodeType** n, ElementType* l, unsigned int level)
//        {
//			if (maxlevel < level) // just for debug
//				maxlevel = level;
//			
//			OctreePosition op = l->checkPosition((*n)->getAABB());
//			assert(op != OP_UP);
//			if (op == OP_OK || level >= levelLimit)
//			{
//				Vector_3 lcenter = l->getCenter();
//				if ((*n)->getAABB().contains(lcenter)) // if the center is not contained in the AABB that it belongs to a cell in the neighbourhood at the same level
//				{
//					NodeType* ptr = *n;
//					
//					if ((*n)->getType() == NT_LEAF)
//					{
//						LeafType* in = dynamic_cast<LeafType*>(ptr);
//						assert(in);
//						in->addElement(l);
//					}
//					else 
//					{
//						InternalType* in = dynamic_cast<InternalType*>(ptr);
//						assert(in);
//						in->addElement(l);
//					}
//				}
//			}
//			else 
//			{
//				assert(op == OP_DOWN);
//				if((*n)->getType() == NT_LEAF) // if it is a leaf, we replace the leaf with an internal node
//				{
//					if (level < levelLimit)
//					{
//						splitLeaf(n, level);
//					}
//				}
//
//				assert((*n)->getType() == NT_INTERNAL);
//				InternalType* in = dynamic_cast<InternalType*>(*n);
//				for(int i=0; i < 8; ++i)
//				{
//					AABBType c = in->children[i]->getAABB();
//					if(l->intersect(c))
//					{
//						add_aux(&(in->children[i]), l, level+1);
//					}
//				}
//				
//			}
//        }
//		
//        void splitLeaf(Node<AABBType,ElementType>** n, unsigned int level)
//        {
//            assert(level <= levelLimit);
//			
//            if (maxlevel < level) // just for debug
//                maxlevel = level;
//			
//            assert((*n)->getType() == NT_LEAF);
//            LeafType* ln = dynamic_cast<LeafType*>(*n);
//			
//            AABBType backupCuboid = ln->getAABB();
//			
//            *n = new InternalType(backupCuboid,level,ln);
//            delete ln;
//        }
//		
//		
//        void remove_aux(Node<AABBType,ElementType>** n, ElementType* l)
//        {
//			LeafType* ln = dynamic_cast<LeafType*>(*n);
//			if (l->intersect (*n)->getAABB())
//				ln->removeElement(l);
//
//			
//			if ((*n)->getType() == NT_INTERNAL)
//			{
//				InternalType* in = dynamic_cast<InternalType*>(*n);
//				for(int i=0; i < 8; ++i)
//				{
//					AABBType c = in->children[i]->getAABB();
//					if(l->intersect(c))
//					{
//						remove_aux(&(in->children[i]), l);
//					}
//				}
//			}
//        }
//		
//        void windowQuery_auxAABB(NodeType* n, AABBType& c, list<ElementType*>& ll)
//        {
//			typename list<ElementType*>::iterator it, end;
//			
//			if(n->getType() == NT_LEAF)
//			{
//				LeafType* ln = dynamic_cast<LeafType*>(n);
//				assert(ln);
//				it = ln->objs.begin();
//				end = ln->objs.end(); 
//			}
//			else 
//			{
//				InternalType* ln = dynamic_cast<InternalType*>(n);
//				assert(ln);
//				it = ln->objs.begin();
//				end = ln->objs.end(); 
//			}
//
//			for(; it != end; ++it)
//			{
//				ElementType* l = *it;
//				if(l->intersect(c))
//					ll.push_back(l);
//			}
//			
//			if (n->getType() == NT_INTERNAL)
//			{
//				InternalType* in = dynamic_cast<InternalType*>(n);
//				for(int i=0; i < 8; ++i)
//				{
//					AABBType& cc = in->children[i]->getAABB();
//					if (cc.extendedIntersect(c))
//					{
//						windowQuery_auxAABB(in->children[i],c,ll);
//					}
//				}
//			}
//        }
//				
//        void stats_aux(Node<AABBType,ElementType>* n, unsigned int currMarker, int& elemCount, int& rifCount, int& nodeIntCount, int& nodeLeafCount)
//        {
//            if(n->getType() == NT_LEAF)
//            {
//                ++nodeLeafCount;
//                LeafNode<AABBType,ElementType>* ln = dynamic_cast<LeafNode<AABBType,ElementType>*>(n);
//				
//                typename list<ElementType*>::iterator it;
//                for(it = ln->objs.begin(); it != ln->objs.end(); ++it)
//                {
//                    ElementType* l = *it;
//                    if (l->mark != currMarker)
//                    {
//                        elemCount++;
//                        l->mark = currMarker;
//                    }
//                    rifCount++;
//                }
//            }
//            else
//            {
//                ++nodeIntCount;
//                assert(n->getType() == NT_INTERNAL);
//                InternalType* in = dynamic_cast<InternalType*>(n);
//                for(int i=0; i < 8; ++i)
//                {
//                    stats_aux(in->children[i],currMarker, elemCount, rifCount, nodeIntCount, nodeLeafCount);
//                }
//				
//            }
//        }
//		
//		
//        Node<AABBType,ElementType>* root;
//		int threshold;
//		
//		unsigned short marker;
//		unsigned int levelLimit;
//	};
//	
//	template<class AABBType, class ElementType>
//	class Node
//	{
//	public:
//		Node(NodeType nt, AABBType c) : type(nt), cuboid(c) {};
//		virtual ~Node(){};
//		NodeType getType() { return type;}
//        AABBType& getAABB() { return cuboid; };
//		
//		friend std::ostream& operator<<(std::ostream& os, const Node& p)
//		{
//			os << p.getAABB();
//			return os;
//		}
//		
//		OctreePosition checkSize(ElementType* e)
//		{
//			return e->checkSize(this);
//		}
//		
//	private:
//		NodeType type;
//		AABBType cuboid;
//	};
//	
//	template<class AABBType, class ElementType>
//	class LeafNode : public Node<AABBType,ElementType>
//	{
//	public:
//		LeafNode(AABBType c, NodeType ntype = NT_LEAF) : Node<AABBType,ElementType>(ntype,c) {};
//		virtual ~LeafNode() {};
//		
//		void addElement(ElementType* l)
//		{
//			objs.push_back(l);
//			l->setOctreeNode(this);
//		};
//		
//		void removeElement(ElementType* l) // TODO: not use this, it's slow!
//		{
//			assert(std::find(objs.begin(),objs.end(),l) != objs.end());
//			objs.remove(l);
//			l->clear();
//			delete l;
//		};
//		
//		void removeElement(typename list<ElementType>::iterator l)
//		{
//			objs.erase(l);
//			l->setOctreeNode(0);
//			delete l;
//		};
//		
//		size_t getSize() // TODO: it is a list, not use this!
//		{
//			return objs.size();
//		};
//		
//		list<ElementType*> objs;
//	};
//	
//	template<class AABBType, class ElementType>
//	class InternalNode : public LeafNode<AABBType,ElementType>
//	{
//	public:
//		typedef typename AABBType::ScalarType ScalarType;
//		typedef typename AABBType::VectorType VectorType;
//		typedef LeafNode<AABBType,ElementType> LeafType;
//		typedef InternalNode<AABBType,ElementType> InternalType;
//		InternalNode(AABBType c, int Level, LeafType* leaf) : LeafType(c,NT_INTERNAL), level(Level)
//		{
//			vector<AABBType> splittedAABB;
//			c.split(splittedAABB);
//			
//			for (unsigned i=0; i<8; ++i)
//				children[i] = new LeafType(splittedAABB[i]);
//			
//			importNodeFromLeaf(leaf);
//		}
//		
//		void importNodeFromLeaf(LeafType* leaf)
//		{
//			typename list<ElementType*>::iterator it = leaf->objs.begin();
//			while(it != leaf->objs.end())
//			{
//				LeafType::addElement(*it);
//				++it;
//			}
//		}
//		
//		Node<AABBType,ElementType>* children[8];
//		
//		int level;
//	};
//	
//	
//	template<class AABBType, class ElementType>
//	class FakeOctree
//	{
//		typedef typename AABBType::ScalarType ScalarType;
//		typedef typename AABBType::CoordType CoordType;
//		
//	public:
//		FakeOctree(AABBType c) : domain(c) {};
//		virtual ~FakeOctree() {};
//		
//        void add(ElementType* l)
//        {
//			AABBType lcub = l->getAABB();
//			CoordType min = lcub.min;
//			CoordType max = lcub.max;
//			assert(domain.contains(min));
//			assert(domain.contains(max));
//			
//			labels.push_back(l);
//			
//        }
//		
//        void remove(ElementType* l)
//        {
//			labels.remove(l);
//        }
//		
//        AABBType domain;
//        list<ElementType*> windowQuery(AABBType c)
//        {
//			list<ElementType*> ll;
//			
//			typename list<ElementType*>::iterator it;
//			for(it=labels.begin(); it != labels.end(); ++it)
//			{
//				ElementType* l = *it;
//				if(l->intersect(c))
//					ll.push_back(l);
//			}
//			return ll;
//        }
//		
//	private:
//		
//        list<ElementType> labelsStorage;
//        list<ElementType*> labels;
//	};
//	
//	class OctreeTestSuite
//	{
//	public:
//		static bool runTestSuite();
//	};
}

#endif /* OCTREE_H_ */

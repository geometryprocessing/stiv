#include "Hash.h"

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
#include <set>
#include <map>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include "PrismVolume.h"

using std::vector;
using namespace std;
using namespace __gnu_cxx;

namespace Intersection
{
	using namespace std;
	
	typedef Eigen::Matrix<int,3,1> Point3i;
	
	
	//Hash::Hash(Vector3d min, Vector3d max, double cellSize, double h, double periodicLength)
	Hash::Hash()
	{
            /*
		//m_cellSize = cellSize;
                // TODO: set cell size to be something that periodicLength
                // will be multiple of cell size
		m_cellSize = 0.2;
		m_thickness = h;
                m_periodicLength = periodicLength;
               
                // for periodic case assumption is made that periodicLength is multiple of m_cellSize
                m_periodicSize = periodicLength/m_cellSize;
	
                m_domainMin = min;
		m_domainMax = max;

                if(periodicLength>0)
                    m_gridSize = m_periodicSize;
                else
                {
                    m_gridSize = (int)((max[0] - min[0])/m_cellSize);
                    m_gridSize = std::max( m_gridSize, (int)((max[1] - min[1])/m_cellSize) );
                    m_gridSize = std::max( m_gridSize, (int)((max[2] - min[2])/m_cellSize) );
                }
            */
	}
	
	Hash::~Hash(){};
	
        void Hash::resetHash(Vector3d min, Vector3d max, double cellSize, double h, int nv_resident, double periodicLength)
	{
		//m_cellSize = cellSize;
                // TODO: set cell size to be something that periodicLength
                // will be multiple of cell size
		clear();
		m_cellSize = 0.2;
		m_thickness = h;
                m_nv_resident = nv_resident;
                m_periodicLength = periodicLength;
               
                // for periodic case assumption is made that periodicLength is multiple of m_cellSize
                m_periodicSize = round(periodicLength/m_cellSize);
	
                m_domainMin = min;
		m_domainMax = max;

                if(periodicLength>0)
                    m_gridSize = m_periodicSize;
                else
                {
                    m_gridSize = (int)((max[0] - min[0])/m_cellSize);
                    m_gridSize = std::max( m_gridSize, (int)((max[1] - min[1])/m_cellSize) );
                    m_gridSize = std::max( m_gridSize, (int)((max[2] - min[2])/m_cellSize) );
                }

                if(periodicLength>0)
                {
                    //add(Point3i(m_gridSize+1, m_gridSize+1, m_gridSize+1), -1, m_hash);
                    //add(Point3i(m_gridSize+1, m_gridSize+1, m_gridSize+1), -1, m_hash_edge);
                }
                else
                {
                    //add(Point3i(max[0]/m_cellSize+1, max[1]/m_cellSize+1, max[2]/m_cellSize+1), -1, m_hash);
                    //add(Point3i(max[0]/m_cellSize+1, max[1]/m_cellSize+1, max[2]/m_cellSize+1), -1, m_hash_edge);
                }
                // assert (m_gridSize+2)^3 < (int_max_limit)
                assert( (m_gridSize+2) < cbrt(std::numeric_limits<long>::max()) );
	}

	void Hash::fill(Polyhedron* p, vector<set<unsigned int> > &subMeshes)
	{
		
		m_p = p;

		// precompute AABB
		m_p->initAABB(m_thickness);
		
//#undef USEOPENMP
//#ifdef USEOPENMP
////		// parallel version...
//		int faceSize = p->size_of_facets();
//		int vertSize = p->size_of_vertices();
//		
//		#define N 8		
//		int countersFace[N+1];
//		int countersVert[N+1];
//
//		for(int i=0; i<N;++i)
//		{
//			countersFace[i] = (faceSize/N)*i;
//			countersVert[i] = (vertSize/N)*i;
//		}
//		countersFace[N] = faceSize;
//		countersVert[N] = vertSize;
//		
//		vector<vector<HashItem> > tempVs(N);
//		
//		int i;
//		
//        #pragma omp parallel for
//		for(i=0; i<N; i++)
//		{
//			for(int j=countersFace[i]; j < countersFace[i+1]; ++j)
//			{
//				if ((p->f_array[j]->marked) && (!p->f_array[j]->ignore))
//					addTriangle(p,p->f_array[j],j,tempVs[i]);
//			}
//		}
//
//		m_numberTriangle = p->size_of_facets();
//
//        #pragma omp parallel for
//		for(i=0; i<N; i++)
//		{
//			for(int j=countersVert[i]; j < countersVert[i+1]; ++j)
//			{
//				if (p->v_array[j]->marked)
//					addVertex(p, p->v_array[j], faceSize + j, tempVs[i]);
//			}
//		}
//		
//		int sizeTotal = 0;
//		
//		BOOST_FOREACH(vector<HashItem>& v, tempVs) sizeTotal += v.size();
//			
//		m_hash.resize(sizeTotal);
//		
//		for(unsigned i=0; i<tempVs.size(); ++i)
//			for(unsigned j=0; j<tempVs[i].size(); ++j)
//				m_hash.push_back(HashItem(tempVs[i][j]));
//		
//#else
                // begin of new omp fill
                size_t omp_p = omp_get_max_threads();
                vector< vector<HashItem> > hash_vf(omp_p);
                vector< vector<HashItem> > hash_edge(omp_p);
                m_numberTriangle = p->size_of_facets();
                size_t size_v = p->size_of_vertices();
                size_t size_f = p->size_of_facets();
                size_t size_e = p->e_array.size();
                #pragma omp parallel for
                for(size_t tid=0; tid<omp_p; tid++)
                {
                    size_t av = (size_v*(tid+0))/omp_p;
                    size_t bv = (size_v*(tid+1))/omp_p;
                    size_t af = (size_f*(tid+0))/omp_p;
                    size_t bf = (size_f*(tid+1))/omp_p;
                    size_t ae = (size_e*(tid+0))/omp_p;
                    size_t be = (size_e*(tid+1))/omp_p;
                    for(size_t i=av; i<bv; i++)
                    {
                        addVertex(p, p->v_array[i], m_numberTriangle+p->v_array[i]->getIndex(), hash_vf[tid]);
                    }
                    for(size_t i=af; i<bf; i++)
                    {
                        addTriangle(p, p->f_array[i], p->f_array[i]->getIndex(), hash_vf[tid]);
                    }
                    for(size_t i=ae; i<be; i++)
                    {
			addEdge(p, p->e_array[i], p->e_array[i]->getIndex(), hash_edge[tid]);
                    }
                }
                vector<size_t> hash_cnt(omp_p, 0);
                vector<size_t> hash_dsp(omp_p, 0);
                size_t hashvf_size=0;
                for(size_t tid=0; tid<omp_p; tid++)
                {
                    if(tid)
                        hash_dsp[tid] = hash_vf[tid-1].size() + hash_dsp[tid-1];
                    hash_cnt[tid] = hash_vf[tid].size();
                    hashvf_size += hash_vf[tid].size();
                }
                m_hash.resize(hashvf_size);

                #pragma omp parallel for
                for(size_t tid=0; tid<omp_p; tid++)
                {
                    std::memcpy(&m_hash[hash_dsp[tid]], &hash_vf[tid][0], sizeof(HashItem)*hash_cnt[tid]);
                }

                hash_cnt.resize(omp_p, 0);
                hash_dsp.resize(omp_p, 0);
                size_t hashe_size=0;
                for(size_t tid=0; tid<omp_p; tid++)
                {
                    if(tid)
                        hash_dsp[tid] = hash_edge[tid-1].size() + hash_dsp[tid-1];
                    hash_cnt[tid] = hash_edge[tid].size();
                    hashe_size += hash_edge[tid].size();
                }
                m_hash_edge.resize(hashe_size);

                #pragma omp parallel for
                for(size_t tid=0; tid<omp_p; tid++)
                {
                    std::memcpy(&m_hash_edge[hash_dsp[tid]], &hash_edge[tid][0], sizeof(HashItem)*hash_cnt[tid]);
                }

                // end of new omp fill
		
                /*
		//int count = 0;
		// first insert all triangles
		// indexes 0 ...  (#Tri-1)
		Facet_iterator fit = p->facets_begin();
		while (fit != p->facets_end())
		{
			//if (fit->marked && (!fit->ignore))
				addTriangle(p,&*fit,fit->getIndex(),m_hash);
			++fit;
			//++count;
		}
		
		m_numberTriangle = p->size_of_facets();
		
		// then insert all points
		// indexes #Tri ... (#Tri + #Vert)
		Vertex_iterator vit = p->vertices_begin();
		while (vit != p->vertices_end())
		{
			//if (vit->marked)
				addVertex(p, &*vit, m_numberTriangle+vit->getIndex(), m_hash);
			++vit;
			//++count;
		}
#ifdef DETECTEE
		// then insert all edges (in a separate grid)
		// indexes 0 ... #Edges
		for(int i=0; i<p->e_array.size();++i)
		{
			if (p->e_array[i]->getIndex() != i)
				cerr << "This should not happen! Edge data structure corrupted" << endl;
//			cerr << "Inserting Edge in hash:" << p->e_array[i]->getIndex() << " " << p->e_array[i] << endl;

			addEdge(p, p->e_array[i], p->e_array[i]->getIndex(), m_hash_edge);
			//++count;
		}
#endif
                */
                std::cout<<"m_hash size: "<<m_hash.size()<<std::endl;
                std::cout<<"m_edge size: "<<m_hash_edge.size()<<std::endl;
	}
	
	AABBi Hash::makeAABBi(double min[3], double max[3])
	{
		AABBi aabbi;
                if(m_periodicLength > 0)
                {
                    aabbi.m_min = Point3i(floor(min[0]/m_cellSize),floor(min[1]/m_cellSize),floor(min[2]/m_cellSize));
                    aabbi.m_max = Point3i(ceil(max[0]/m_cellSize),ceil(max[1]/m_cellSize),ceil(max[2]/m_cellSize));
                }
                else
                {
                    aabbi.m_min = Point3i(min[0]/m_cellSize,min[1]/m_cellSize,min[2]/m_cellSize);
                    aabbi.m_max = Point3i(max[0]/m_cellSize,max[1]/m_cellSize,max[2]/m_cellSize);
                }
                return aabbi;
	}
	
	
	void Hash::addTriangle(Polyhedron* p, Facet* fit, int id, vector<HashItem>& hs)
	{		
                AABBi aabbi = makeAABBi(fit->min,fit->max);
                add(aabbi,id, hs);
	}
	
	void Hash::addVertex(Polyhedron* p, Vertex* vit, int id, vector<HashItem>& hs)
	{
		AABBi aabbi = makeAABBi(vit->min,vit->max);
		add(aabbi,id,hs);
	}

	void Hash::addEdge(Polyhedron* p, Halfedge* eit, int id, vector<HashItem>& hs)
	{
		AABBi aabbi = makeAABBi(eit->min,eit->max);
		add(aabbi,id,hs);
	}
	
	//typedef std::pair<int, int> intPair;
	void Hash::sort()
	{
		//std::sort(m_hash.begin(),m_hash.end());
		//std::sort(m_hash_edge.begin(),m_hash_edge.end());
            __gnu_parallel::sort(m_hash.begin(),m_hash.end());
            __gnu_parallel::sort(m_hash_edge.begin(),m_hash_edge.end());
	}
	
	struct eqstr
	{
		bool operator()(Vertex* s1, Vertex* s2) const
		{
			return s1 == s2;
		}
	};
	
	class CollisionOpenMpContainer
	{
	public:
		Vertex* v;
		Facet* f;
		Collision c;
		bool full;
	};

        class CollisionOpenMpEEContainer
        {
            public:
                Halfedge* e1;
                Halfedge* e2;
                Collision c;
                bool full;
        };
	
	void Hash::computeVolume(vector<Collision>& collisions, vector<Collision>& collisionsEE)
	{
                double t_vf_start = omp_get_wtime();
	        ////cout<<"begin face-vertex detection"<<endl;	
		// Compute Intersections for face-point pairs
		unsigned size = m_hash.size();
	        ////cout<<"v-f hash map size: "<<size<<endl;	
		assert(size != 0);
		
		long prevH = get(0).key;
		//static Facet* facet[10000];
		//static Vertex* vs[10000];

		int fiSize = 0;
		int viSize = 0;
		
//		vector<Facet_iterator> faces;
//		vector<Vertex_iterator> vs;

		typedef boost::unordered_set<Facet*> setType;
		typedef boost::unordered_map<Vertex*, setType > mapType;
		mapType psCopy;
                
                int icount=0;
                double t_vf_broad_start = omp_get_wtime();
		for(unsigned i=0; i < size; ++i)
		{
			int& currId = get(i).id;
			long& currH = get(i).key;
			
			if ( (prevH != currH) || (i == (size-1)) ) // new bucket
			{
                                //std::cout<<"fSize: "<<fiSize<<" ,vSize: "<<viSize<<std::endl;
				for(int vi=0; vi<viSize; ++vi)
				{
					Vertex* vit = vs[vi]; 
					for(int fi=0; fi<fiSize; ++fi)
					{
						Facet* fit = facet[fi];
                                                if( (m_p->getVertex(fit->vi[0])->meshTag>=m_nv_resident && vit->meshTag>=m_nv_resident) || (m_p->getVertex(fit->vi[0])->meshTag == vit->meshTag) )
                                                  continue;
						if (PrismVolume::validateCollisionFastFast(fit,vit))
						{
							mapType::iterator itr = psCopy.find(&*vit);
							bool found = false;
							if (itr != psCopy.end())
								found = itr->second.find(&*fit) != itr->second.end();

							if(!found && PrismVolume::validateCollisionFast(fit,vit,m_p,m_periodicLength)){
								psCopy[&*vit].insert(&*fit);
                                                                ++icount;
                                                        }
						}
					}
				}
				
				// intersections!
				fiSize = 0;
				viSize    = 0;
			}

			prevH = currH;
			if ((currId < m_numberTriangle) && (currId != -1))
			{
				facet[fiSize] = &*(m_p->facets_begin()+currId);
				++fiSize;
				assert(fiSize < 10000);
			}
			else 
			{
				vs[viSize] = &*(m_p->vertices_begin()+(currId-m_numberTriangle));
				++viSize;
				assert(viSize < 10000);
			}
                        //std::cout<<"fisize: "<<fiSize<<", visize: "<<viSize<<std::endl;
		}
                double t_vf_broad_end = omp_get_wtime();
                cout<<"total time for f-v collision detection broad phase: "<<(t_vf_broad_end - t_vf_broad_start)<<endl;
//#ifdef DEBUGVG
                ////cout<<"total inserted for f-v: "<<icount<<endl;
//#endif	
                double t_vf_narrow_start = omp_get_wtime();
		vector<CollisionOpenMpContainer> openmpVector;
		openmpVector.reserve(psCopy.size() * 5);
		
		BOOST_FOREACH(mapType::value_type& pvf, psCopy)
		{
			BOOST_FOREACH(Facet* f, pvf.second)
			{
				CollisionOpenMpContainer temp;
				temp.v = pvf.first;
				temp.f = f;
				temp.full = false;
				openmpVector.push_back(temp);
			}
		}
		
                std::cout<<"openmpvector size: "<<openmpVector.size()<<std::endl;
		
                #pragma omp parallel for
		for(int i=0; i<openmpVector.size(); ++i)
		{
			CollisionOpenMpContainer& temp = openmpVector[i];
			temp.full = PrismVolume::validateCollisionAndReplace(temp.f,temp.v, m_thickness, temp.c, m_p, m_periodicLength);
		}
		
		for(int i=0; i<openmpVector.size(); ++i)
		{
			CollisionOpenMpContainer& temp = openmpVector[i];
			if (temp.full)
				collisions.push_back(temp.c);
		}
                double t_vf_narrow_end = omp_get_wtime();
                cout<<"total time for f-v collision detection narrow phase: "<<(t_vf_narrow_end - t_vf_narrow_start)<<endl;
		
                std::cout<<"total validate f-v: "<<collisions.size()<<std::endl;
                double t_vf_end = omp_get_wtime();
                ////cout<<"total time for f-v collision detection: "<<(t_vf_end - t_vf_start)<<endl;
//		BOOST_FOREACH(mapType::value_type& pvf, psCopy)
//		{
//			BOOST_FOREACH(Facet* f, pvf.second)
//			{
//				PrismVolume::validateCollisionAndAdd(f,pvf.first, m_thickness, collisions);
//			}
//		}
		
#ifdef DETECTEE
		// Compute Intersections for edge-edge pairs;
                double t_ee_start = omp_get_wtime();
	        ////cout<<"begin edge-egde detection"<<endl;	
		size = m_hash_edge.size();
	        ////cout<<"e-e hash map size: "<<size<<endl;	
		assert(size != 0);
		prevH = getEdge(0).key;
		
		//static Halfedge* es[10000];
		int eSize = 0;
		
		typedef boost::unordered_set<Halfedge*> setTypeE;
		typedef boost::unordered_map<Halfedge*, setTypeE > mapTypeE;
		mapTypeE eCopy;

                icount=0;
                double t_ee_broad_start = omp_get_wtime();
		for(unsigned i=0; i < size; ++i)
		{
			int& currId = getEdge(i).id;
			long& currH = getEdge(i).key;
			
			if ( (prevH != currH) || (i == (size-1)) ) // new bucket
			{
                                //std::cout<<"eSize: "<<eSize<<std::endl;
				for(int ei=0; ei<eSize; ++ei)
				{
					Halfedge* eip1 = es[ei];
					for(int ei2=(ei+1); ei2<eSize; ++ei2)
					{
						Halfedge* eip2 = es[ei2];
                                                if( (eip1->vertex()->meshTag>=m_nv_resident && eip2->vertex()->meshTag>=m_nv_resident) || (eip1->vertex()->meshTag == eip2->vertex()->meshTag) )
                                                  continue;
						//if (eip1 != eip2 && PrismVolume::validateCollisionEEFastFast(eip1,eip2))
                                                if(1)
						{
						        if(1){
					                Halfedge* eip1tmp = eip1;
					                Halfedge* eip2tmp = eip2;
                                                        if (eip1tmp > eip2tmp)
								std::swap(eip1tmp,eip2tmp);
                                                        mapTypeE::iterator itr = eCopy.find(&*eip1tmp);
							bool found = false;
							if (itr != eCopy.end())
								found = (itr->second.find(&*eip2tmp) != itr->second.end());
							
                                                        if(!found && PrismVolume::validateCollisionEEFast(eip1tmp,eip2tmp,m_p,m_periodicLength,m_thickness)){
								eCopy[&*eip1tmp].insert(&*eip2tmp);
                                                                ++icount;
                                                        }
                                                        } 

                                                        if(0){
                                                        if (eip1 > eip2)
								std::swap(eip1,eip2);
                                                        mapTypeE::iterator itr = eCopy.find(&*eip1);
							bool found = false;
							if (itr != eCopy.end())
								found = itr->second.find(&*eip2) != itr->second.end();
							
                                                        if(!found && PrismVolume::validateCollisionEEFast(eip1,eip2,m_p))
								eCopy[&*eip1].insert(&*eip2);
                                                        }
						}
					}
				}
				// intersections!
				eSize = 0;
			}
			
			prevH = currH;
                        if(currId != -1)
                        {
                            es[eSize] = &*(m_p->e_array[currId]);
                            ++eSize;
                            assert(eSize < 10000);
                            //std::cout<<"esize: "<<eSize<<std::endl;
                        }
		}
                double t_ee_broad_end = omp_get_wtime();
                cout<<"total time for e-e collision broad detection: "<<(t_ee_broad_end - t_ee_broad_start)<<endl;
//#ifdef DEBUGVG
                ////cout<<"total inserted for e-e: "<<icount<<endl;
//#endif	
                double t_ee_narrow_start = omp_get_wtime();
                vector<CollisionOpenMpEEContainer> openmpEEVector;
                openmpEEVector.reserve(eCopy.size() * 10); // 10?

                BOOST_FOREACH(mapTypeE::value_type& pee, eCopy)
                {
                    //cout<<"edge has: "<<pee.second.size()<<" candidate edges"<<endl;
                    BOOST_FOREACH(Halfedge* e, pee.second)
                    {
                        CollisionOpenMpEEContainer temp;
                        temp.e1 = pee.first;
                        temp.e2 = e;
                        temp.full = false;
                        openmpEEVector.push_back(temp);
                    }
                }

                std::cout<<"openmpEEvector size: "<<openmpEEVector.size()<<std::endl;

                #pragma omp parallel for
                for(int i=0; i<openmpEEVector.size(); ++i)
                {
                    CollisionOpenMpEEContainer& temp = openmpEEVector[i];
                    temp.full = PrismVolume::validateCollisionEEAndReplace(temp.e1, temp.e2, m_thickness, temp.c, m_p, m_periodicLength);
                }

                for(int i=0; i<openmpEEVector.size(); ++i)
                {
                    CollisionOpenMpEEContainer& temp = openmpEEVector[i];
                    if(temp.full)
                        collisionsEE.push_back(temp.c);
                }
                /*
		BOOST_FOREACH(mapTypeE::value_type& pvf, eCopy)
		{
			Halfedge* e1 = pvf.first;
			
			BOOST_FOREACH(Halfedge* e2, pvf.second)
			{ 
                                Collision ce;
				if (PrismVolume::validateCollisionEEAndReplace(e1, e2, m_thickness, ce, m_periodicLength))
				{
					collisionsEE.push_back(ce);
				}
			}
		}
                */
                double t_ee_narrow_end = omp_get_wtime();
                cout<<"total time for e-e collision narrow detection: "<<(t_ee_narrow_end - t_ee_narrow_start)<<endl;

                std::cout<<"total validate e-e: "<<collisionsEE.size()<<std::endl;
                double t_ee_end = omp_get_wtime();
                ////cout<<"total time for e-e collision detection: "<<(t_ee_end - t_ee_start)<<endl;
		//cerr << "EE det: " << collisionsEE.size() << endl;
#endif
	}
	
        void Hash::computeVolumeParallel(vector<Collision>& collisions, vector<Collision>& collisionsEE)
	{
                double t_vf_start = omp_get_wtime();
		unsigned size = m_hash.size();
		assert(size != 0);

                double t_vf_broad_start = omp_get_wtime();
                // new omp parallel broad phase
                size_t omp_p = omp_get_max_threads();
                std::vector< std::pair<Vertex*, Facet*> > v_f_pairs[omp_p];
                #pragma omp parallel for
                for(size_t tid=0; tid<omp_p; tid++)
                {
                    size_t a = (size*(tid+0))/omp_p;
                    size_t b = (size*(tid+1))/omp_p;
                    if(a>0    ) a = std::lower_bound(&m_hash[0], &m_hash[0]+size, m_hash[a]) - &m_hash[0];
                    if(b<size ) b = std::lower_bound(&m_hash[0], &m_hash[0]+size, m_hash[b]) - &m_hash[0];
                    for(size_t i=a; i<b;)
                    {
                        size_t j=i; while(j<b && get(j).key==get(i).key) j++;

                        for(size_t id_i=i; id_i<j; id_i++)
                        {
                            int hashid_i = m_hash[id_i].id;
                            for(size_t id_j=id_i+1; id_j<j; id_j++)
                            {
                                int hashid_j = m_hash[id_j].id;
                                if( ( (hashid_i<m_numberTriangle) && (hashid_j<m_numberTriangle) ) ||
                                    ( (hashid_i>=m_numberTriangle) && (hashid_j>=m_numberTriangle) )
                                  )
                                    continue;

                                Vertex* vit;
                                Facet* fit;

                                if(hashid_i>=m_numberTriangle)
                                {
                                    vit = &*(m_p->vertices_begin()+(hashid_i - m_numberTriangle));
                                    fit = &*(m_p->facets_begin()+hashid_j);
                                }
                                else
                                {
                                    vit = &*(m_p->vertices_begin()+(hashid_j - m_numberTriangle));
                                    fit = &*(m_p->facets_begin()+hashid_i);
                                }


                                if( (m_p->getVertex(fit->vi[0])->meshTag>=m_nv_resident && vit->meshTag>=m_nv_resident) || (m_p->getVertex(fit->vi[0])->meshTag == vit->meshTag) )
                                    continue;
                                if (PrismVolume::validateCollisionFastFast(fit,vit))
                                {
                                    if(PrismVolume::validateCollisionFast(fit,vit,m_p,m_periodicLength))
                                        v_f_pairs[tid].push_back(pair<Vertex*, Facet*>(vit, fit));
                                }
                            }
                        }
                        i = j;
                    }
                    std::sort(v_f_pairs[tid].begin(), v_f_pairs[tid].end());
                    v_f_pairs[tid].erase( std::unique(v_f_pairs[tid].begin(), v_f_pairs[tid].end()), v_f_pairs[tid].end());
                }

                size_t vf_size=0;
                vector<size_t> pair_cnt(omp_p, 0);
                vector<size_t> pair_dsp(omp_p, 0);
                for(size_t tid=0; tid<omp_p; tid++)
                {
                    if(tid)
                        pair_dsp[tid] = v_f_pairs[tid-1].size() + pair_dsp[tid-1];
                    pair_cnt[tid] = v_f_pairs[tid].size();
                    vf_size += v_f_pairs[tid].size();
                }
                std::vector< std::pair<Vertex*, Facet*> > total_v_f_pairs(vf_size), total_v_f_pairs_unique;
                #pragma omp parallel for
                for(size_t tid=0; tid<omp_p; tid++)
                {
                    size_t dsp = pair_dsp[tid];
                    size_t cnt = pair_cnt[tid];
                    std::vector<std::pair<Vertex*, Facet*> >& pair_tmp = v_f_pairs[tid];
                    for(size_t i=0; i<cnt; i++)
                    {
                        total_v_f_pairs[(dsp+i)].first = pair_tmp[i].first;
                        total_v_f_pairs[(dsp+i)].second = pair_tmp[i].second;
                    }
                }
                __gnu_parallel::sort(total_v_f_pairs.begin(), total_v_f_pairs.end());
                __gnu_parallel::unique_copy(total_v_f_pairs.begin(), total_v_f_pairs.end(), back_inserter(total_v_f_pairs_unique));
                // end of new omp parallel broad phase
                double t_vf_broad_end = omp_get_wtime();
                cout<<"total time for f-v collision detection broad phase: "<<(t_vf_broad_end - t_vf_broad_start)<<endl;
                
                double t_vf_narrow_start = omp_get_wtime();
		vector<CollisionOpenMpContainer> openmpVector(total_v_f_pairs_unique.size());
                #pragma omp parallel for
                for(size_t i=0; i<openmpVector.size(); i++)
                {
                    openmpVector[i].v = total_v_f_pairs_unique[i].first;
                    openmpVector[i].f = total_v_f_pairs_unique[i].second;
                    openmpVector[i].full = false;
                }
                std::cout<<"openmpvector size: "<<openmpVector.size()<<std::endl;
		
                #pragma omp parallel for
		for(int i=0; i<openmpVector.size(); ++i)
		{
			CollisionOpenMpContainer& temp = openmpVector[i];
			temp.full = PrismVolume::validateCollisionAndReplace(temp.f,temp.v, m_thickness, temp.c, m_p, m_periodicLength);
		}
		
		for(int i=0; i<openmpVector.size(); ++i)
		{
			CollisionOpenMpContainer& temp = openmpVector[i];
			if (temp.full)
				collisions.push_back(temp.c);
		}
                double t_vf_narrow_end = omp_get_wtime();
                cout<<"total time for f-v collision detection narrow phase: "<<(t_vf_narrow_end - t_vf_narrow_start)<<endl;
		
                std::cout<<"total validate f-v: "<<collisions.size()<<std::endl;
                double t_vf_end = omp_get_wtime();
		
#ifdef DETECTEE
                double t_ee_start = omp_get_wtime();
		size = m_hash_edge.size();
		assert(size != 0);
                
                double t_ee_broad_start = omp_get_wtime();
                // new omp parallel broad phase
                std::vector< std::pair<Halfedge*, Halfedge*> > e_e_pairs[omp_p];
                #pragma omp parallel for
                for(size_t tid=0; tid<omp_p; tid++)
                {
                    size_t a = (size*(tid+0))/omp_p;
                    size_t b = (size*(tid+1))/omp_p;
                    if(a>0    ) a = std::lower_bound(&m_hash_edge[0], &m_hash_edge[0]+size, m_hash_edge[a]) - &m_hash_edge[0];
                    if(b<size ) b = std::lower_bound(&m_hash_edge[0], &m_hash_edge[0]+size, m_hash_edge[b]) - &m_hash_edge[0];
                    for(size_t i=a; i<b;)
                    {
                        size_t j=i; while(j<b && getEdge(j).key==getEdge(i).key) j++;

                        for(size_t id_i=i; id_i<j; id_i++)
                        {
                            int hashid_i = m_hash_edge[id_i].id;
                            for(size_t id_j=id_i+1; id_j<j; id_j++)
                            {
                                int hashid_j = m_hash_edge[id_j].id;

                                Halfedge* eip1; 
                                Halfedge* eip2;

                                if(hashid_i>hashid_j)
                                {
                                    eip1 = &*(m_p->e_array[hashid_i]);
                                    eip2 = &*(m_p->e_array[hashid_j]);
                                }
                                else
                                {
                                    eip2 = &*(m_p->e_array[hashid_i]);
                                    eip1 = &*(m_p->e_array[hashid_j]);
                                }

                                if( (eip1->vertex()->meshTag>=m_nv_resident && eip2->vertex()->meshTag>=m_nv_resident) || (eip1->vertex()->meshTag == eip2->vertex()->meshTag) )
                                    continue;

                                if(PrismVolume::validateCollisionEEFast(eip1,eip2,m_p,m_periodicLength,m_thickness))
                                {
                                    e_e_pairs[tid].push_back(pair<Halfedge*, Halfedge*>(eip1, eip2));
                                }

                            }
                        }
                        i = j;
                    }
                    std::sort(e_e_pairs[tid].begin(), e_e_pairs[tid].end());
                    e_e_pairs[tid].erase( std::unique(e_e_pairs[tid].begin(), e_e_pairs[tid].end()), e_e_pairs[tid].end());
                }

                size_t ee_size=0;
                pair_cnt.resize(omp_p, 0);
                pair_dsp.resize(omp_p, 0);
                for(size_t tid=0; tid<omp_p; tid++)
                {
                    if(tid)
                        pair_dsp[tid] = e_e_pairs[tid-1].size() + pair_dsp[tid-1];
                    pair_cnt[tid] = e_e_pairs[tid].size();
                    ee_size += e_e_pairs[tid].size();
                }
                std::vector< std::pair<Halfedge*, Halfedge*> > total_e_e_pairs(ee_size), total_e_e_pairs_unique;
                #pragma omp parallel for
                for(size_t tid=0; tid<omp_p; tid++)
                {
                    size_t dsp = pair_dsp[tid];
                    size_t cnt = pair_cnt[tid];
                    std::vector<std::pair<Halfedge*, Halfedge*> >& pair_tmp = e_e_pairs[tid];
                    for(size_t i=0; i<cnt; i++)
                    {
                        total_e_e_pairs[(dsp+i)].first = pair_tmp[i].first;
                        total_e_e_pairs[(dsp+i)].second = pair_tmp[i].second;
                    }
                }
                __gnu_parallel::sort(total_e_e_pairs.begin(), total_e_e_pairs.end());
                __gnu_parallel::unique_copy(total_e_e_pairs.begin(), total_e_e_pairs.end(), back_inserter(total_e_e_pairs_unique));
                // end of new omp parallel broad phase
                double t_ee_broad_end = omp_get_wtime();
                cout<<"total time for e-e collision broad detection: "<<(t_ee_broad_end - t_ee_broad_start)<<endl;

                double t_ee_narrow_start = omp_get_wtime();
                vector<CollisionOpenMpEEContainer> openmpEEVector(total_e_e_pairs_unique.size());
                #pragma omp parallel for
                for(size_t i=0; i<openmpEEVector.size(); i++)
                {
                    openmpEEVector[i].e1 = total_e_e_pairs_unique[i].first;
                    openmpEEVector[i].e2 = total_e_e_pairs_unique[i].second;
                    openmpEEVector[i].full = false;
                }
                std::cout<<"openmpEEvector size: "<<openmpEEVector.size()<<std::endl;

                #pragma omp parallel for
                for(int i=0; i<openmpEEVector.size(); ++i)
                {
                    CollisionOpenMpEEContainer& temp = openmpEEVector[i];
                    temp.full = PrismVolume::validateCollisionEEAndReplace(temp.e1, temp.e2, m_thickness, temp.c, m_p, m_periodicLength);
                }

                for(int i=0; i<openmpEEVector.size(); ++i)
                {
                    CollisionOpenMpEEContainer& temp = openmpEEVector[i];
                    if(temp.full)
                        collisionsEE.push_back(temp.c);
                }
                double t_ee_narrow_end = omp_get_wtime();
                cout<<"total time for e-e collision narrow detection: "<<(t_ee_narrow_end - t_ee_narrow_start)<<endl;

                std::cout<<"total validate e-e: "<<collisionsEE.size()<<std::endl;
                double t_ee_end = omp_get_wtime();
#endif
	}
	
	void Hash::add(AABBi aabbi, int id, vector<HashItem>& hs)
	{
		for(int x=aabbi.getMin().x(); x<=aabbi.getMax().x(); ++x)
			for(int y=aabbi.getMin().y(); y<=aabbi.getMax().y(); ++y)
				for(int z=aabbi.getMin().z(); z<=aabbi.getMax().z(); ++z)
				{
                                        if(m_periodicLength > 0)
                                        {
                                            /*
                                            if(x >= m_periodicSize/2)
                                                x -= m_periodicSize;
                                            if(x < -m_periodicSize/2)
                                                x += m_periodicSize;

                                            if(y >= m_periodicSize/2)
                                                y -= m_periodicSize;
                                            if(y < -m_periodicSize/2)
                                                y += m_periodicSize;

                                            if(z >= m_periodicSize/2)
                                                z -= m_periodicSize;
                                            if(z < -m_periodicSize/2)
                                                z += m_periodicSize;

                                        
                                            add(Point3i(x,y,z),id, hs);
                                            */

                                            add(Point3i((x+m_periodicSize)%m_periodicSize,
                                                        (y+m_periodicSize)%m_periodicSize,
                                                        (z+m_periodicSize)%m_periodicSize),id, hs);
                                        }
                                        else
                                        {
                                            add(Point3i(x,y,z),id, hs);
                                        }
				}
	}
	
	void Hash::add(Point3i p, int id, vector<HashItem>& hs)
	{
		hs.push_back(HashItem(hash(p),id));
		//cerr << p.transpose() << " (" << hash(p) << ") " << "---->" << id << endl;

	}
	
	long Hash::hash(Point3i p)
	{
		return ((long)p.z()) * ((long)m_gridSize) * ((long)m_gridSize) + ((long)p.y())*((long)m_gridSize) + ((long)p.x());
	}
	
	void Hash::clear()
	{
//		m_hashId.clear();
//		m_hashCode.clear();
		m_hash.clear();
		m_hash_edge.clear();
	}
	
	HashItem& Hash::get(int i)
	{
		return m_hash[i];
	}

	HashItem& Hash::getEdge(int i)
	{
		return m_hash_edge[i];
	}
	
//	void Hash::insert(HashItem hi)
//	{
//		m_hash.push_back(hi);
//	}
}

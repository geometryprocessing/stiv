// Surface.cpp
//

#include <iostream>
#include "ColSurface.h"
#include "PrismVolume.h"
#include "lcp.h"
#include "ode_common.h"
#include <list>
#include <set>
#include <queue>
#include <iomanip>
#include <algorithm>
#include <boost/foreach.hpp>
#include "timer.h"
#include <boost/unordered_set.hpp>

using namespace std;

ColSurface::ColSurface(Polyhedron *p)
 : m_controlMesh(p), m_slaveMesh(0), m_aabb(0), m_rigid(false), m_moveAll(false)
{
        m_thickness = 0.00;
        ////cout<<"thickness is: "<<m_thickness<<endl;
}

ColSurface::~ColSurface()
{
        if(m_controlMesh == m_slaveMesh)
        {
            if (m_controlMesh)
		delete m_controlMesh;
        }
        else
        {
            if (m_controlMesh)
		delete m_controlMesh;

            if (m_slaveMesh)
		delete m_slaveMesh;
        }

	if (m_aabb)
		delete m_aabb;
}

void ColSurface::initialize(bool precompRigigGradients)
{
	// No Slave mesh? Just make a copy of the control mesh
	//
	if (!getSlaveMesh())
	{
		m_slaveMesh = new Polyhedron(*getControlMesh());
		m_slaveMesh->reset();

		Vertex_iterator vit = m_slaveMesh->vertices_begin();
		while (vit != m_slaveMesh->vertices_end())
		{
			vit->Pori() = vit->P();
			assert(vit->Pori() == vit->P());
			vit++;
		}
	}
        
	if (precompRigigGradients)
	{
            precomputeRigidGradients();
	}
	
	if (!m_aabb)
	{
		// INIT BOUNDING BOX
		vector<Vector3d> temp;
		temp.reserve(m_slaveMesh->size_of_vertices());
		Vertex_iterator vit = m_slaveMesh->vertices_begin();
		while (vit != m_slaveMesh->vertices_end())
		{
			temp.push_back(vit->P());
			++vit;
		}

		AABB aabb(temp);
		vector<double> sizes;
		sizes.push_back(aabb.getSize()[0]);
		sizes.push_back(aabb.getSize()[1]);
		sizes.push_back(aabb.getSize()[2]);
		double m = *max_element(sizes.begin(), sizes.end());
		m *= 1.1;
                
                m = max(m, m_periodicLength);
		
		m_aabb = new AABB(aabb.getCenter(),Vector3d(m,m,m));
		// END BOUNDING BOX

		// AVERAGE EDGE LENGTH
		
		double avgl = 0;
		double minl = 1e50;
		double maxl = 0.0;
		Halfedge_iterator heit = m_controlMesh->halfedges_begin();
		while (heit != m_controlMesh->halfedges_end())
		{
			double l = (heit->vertex()->P() - heit->opposite()->vertex()->P()).norm();
			if (l < minl)
				minl = l;
			if (l > maxl)
				maxl = l;
			
                        avgl += (heit->vertex()->P() - heit->opposite()->vertex()->P()).norm();

                        // record original length
                        heit->length = l;
                        heit->opposite()->length = l;
                        //cout << "this half edge is: "<<heit->length<<endl;

			++heit;
		}
		avgl = avgl / m_controlMesh->size_of_halfedges();
		m_edgeLength = avgl;
                cout << "NUM of E: " << m_controlMesh->size_of_halfedges() <<  endl;
		cout << "AVG: " << avgl << endl;
		cout << "MIN: " << minl << endl;
		cout << "MAX: " << maxl << endl;

		// END AVERAGE EDGE LENGTH
	}
}

void ColSurface::precomputeRigidGradients()
{
        m_subMeshVertices.clear(); m_subMeshVertices.resize(m_nv+m_np);
        m_subMeshes.clear(); m_subMeshes.resize(m_nv+m_np);

        std::cout<<"vertex size: "<<m_controlMesh->size_of_vertices()<<"\n";
        std::cout<<"face size: "<<m_controlMesh->size_of_facets()<<"\n";

        int p_per_v = (m_p+1)*m_p*2 + 2;
        int f_per_v = 2*2*m_p*m_p + 2*m_p*2;
        int p_per_p = m_pp*m_pp;
        int f_per_p = (m_pp-1)*(m_pp-1)*2;
        
        std::cout<<p_per_v*m_nv+p_per_p*m_np<<", "<<f_per_v*m_nv+f_per_p*m_np<<"\n";

        double t_vtag_start = omp_get_wtime();
        #pragma omp parallel for
        for(size_t i=0; i<m_nv; i++)
        {
            unsigned int basev = i*p_per_v;
            unsigned int basef = i*f_per_v;
            set<unsigned int> &currMeshv = m_subMeshVertices[i];
            set<unsigned int> &currMeshf = m_subMeshes[i];
            //fill vesicle i vertices
            for(size_t j=0; j<p_per_v; j++)
            {
                unsigned int v_ind = basev + j;
                m_slaveMesh->getVertex(v_ind)->meshTag = i;
                currMeshv.insert(v_ind);
            }
            //fill vesicle i faces
            for(size_t j=0; j<f_per_v; j++)
            {
                currMeshf.insert(basef+j);
            }
            //m_subMeshVertices.push_back(currMeshv);
            //m_subMeshes.push_back(currMeshf);
        }
        double t_vtag_end = omp_get_wtime();
        std::cout<<"vtag build time: "<<t_vtag_end - t_vtag_start<<std::endl;

        double t_ftag_start = omp_get_wtime();
        #pragma omp parallel for
        for(size_t i=0; i<m_np; i++)
        {
            unsigned int basev = m_nv*p_per_v + i*p_per_p;
            unsigned int basef = m_nv*f_per_v + i*f_per_p;
            set<unsigned int> &currMeshv = m_subMeshVertices[m_nv+i];
            set<unsigned int> &currMeshf = m_subMeshes[m_nv+i] ;

            for(size_t j=0; j<p_per_p; j++)
            {
                unsigned int v_ind = basev + j;
                m_slaveMesh->getVertex(v_ind)->meshTag = m_nv;
                currMeshv.insert(v_ind);
            }
            for(size_t j=0; j<f_per_p; j++)
            {
                currMeshf.insert(basef+j);
            }
            //m_subMeshVertices.push_back(currMeshv);
            //m_subMeshes.push_back(currMeshf);
        }
        double t_ftag_end = omp_get_wtime();
        std::cout<<"ftag build time: "<<t_ftag_end - t_ftag_start<<std::endl;

        std::cout<<"number of vesicles: "<<m_nv<<std::endl;
        cout << m_subMeshVertices.size() << " submeshes. " << endl;

        /*
	// Group mesh vertices into disjoint sub-meshes
	//
	//vector<set<unsigned int> > subMeshes;
	set<unsigned int> vertices;
	for (size_t i=0; i<m_controlMesh->size_of_vertices(); ++i)
		vertices.insert(i);

	m_subMeshVertices.clear();
	m_subMeshes.clear();
	
        int meshTag = 0;
        std::cout<<"number of vesicles: "<<m_nv<<std::endl;
        double t_vtag_start = omp_get_wtime();
	while (!vertices.empty())
	{
                std::cout<<"meshTag: "<<meshTag<<std::endl;
		unsigned int seed = *vertices.begin();
		set<unsigned int> currMesh;

		queue<unsigned int> q;
		q.push(seed);

		currMesh.insert(seed);
                m_slaveMesh->getVertex(seed)->meshTag = meshTag;

		while (!q.empty())
		{
			unsigned int curr = q.front();
			q.pop();

			Vertex_handle vh = m_controlMesh->vertices_begin() + curr;
			Halfedge_Vertex_circulator hvcit = vh->vertex_begin();
			do 
			{
				Vertex_handle neighbor_v = hvcit->opposite()->vertex();
				// Have we already visited this vertex?
				//
				if (!currMesh.count(neighbor_v->getIndex()))
				{
					q.push(neighbor_v->getIndex());
					currMesh.insert(neighbor_v->getIndex());
                                        m_slaveMesh->getVertex(neighbor_v->getIndex())->meshTag = meshTag;
				}

				hvcit++;
			} while (hvcit != vh->vertex_begin());
		}

		m_subMeshVertices.push_back(currMesh);

		vector<unsigned int> newVertices(vertices.size() - currMesh.size());
		std::set_difference(vertices.begin(), vertices.end(), currMesh.begin(), currMesh.end(), newVertices.begin());
		vertices.clear();
		vertices.insert(newVertices.begin(), newVertices.end());
                
                meshTag++;
                if(meshTag >= m_nv)
                    meshTag = m_nv;
	}
        double t_vtag_end = omp_get_wtime();
        std::cout<<"vtag build time: "<<t_vtag_end - t_vtag_start<<std::endl;


        double t_ftag_start = omp_get_wtime();
	// There's probably a better way to do this, but just do another flood-fill on the faces of the slave mesh
	//
	set<unsigned int> faces;
	for (size_t i=0; i<m_slaveMesh->size_of_facets(); ++i)
		faces.insert(i);
	
	while (!faces.empty())
	{
		unsigned int seed = *faces.begin();
		set<unsigned int> currMesh;

		queue<unsigned int> q;
		q.push(seed);

		currMesh.insert(seed);

		while (!q.empty())
		{
			unsigned int curr = q.front();
			q.pop();

			Halfedge_Facet_circulator hfcit = m_slaveMesh->getFacet(curr)->facet_begin();
			do 
			{
				if (!hfcit->opposite()->is_border())
				{
					Facet_handle neighbor_f = hfcit->opposite()->facet();

					// Have we already visited this vertex?
					//
					if (!currMesh.count(neighbor_f->getIndex()))
					{
						q.push(neighbor_f->getIndex());
						currMesh.insert(neighbor_f->getIndex());
					}
				}

				hfcit++;
			} while (hfcit != m_slaveMesh->getFacet(curr)->facet_begin());
		}

		m_subMeshes.push_back(currMesh);

		vector<unsigned int> newFaces(faces.size() - currMesh.size());
		std::set_difference(faces.begin(), faces.end(), currMesh.begin(), currMesh.end(), newFaces.begin());
		faces.clear();
		faces.insert(newFaces.begin(), newFaces.end());
	}
        double t_ftag_end = omp_get_wtime();
        std::cout<<"ftag build time: "<<t_ftag_end - t_ftag_start<<std::endl;
	
        cout << m_subMeshVertices.size() << " submeshes. " << endl;
        */
}

void ColSurface::controlMeshUpdateIteration(set<unsigned int> &selectedVertices, map<unsigned int, unsigned int> &colliding, 
    bool handleCollisions, bool handleElasticity, Vector_3 &translate, bool rigid, bool moveAll, int limiter, 
    Eigen::SparseVector<double> &gradient, std::vector<double> &IV, Eigen::SparseVector<int> &gradient_index, int &IV_size, 
    std::vector<double> &IS_BDV)
{
	m_rigid = rigid;
	
	m_moveAll = moveAll;
	
        computeNormals(m_slaveMesh); 	

	if (handleCollisions)
	{
		int iteration = 0;
		vector<Collisions> collisions;
		
                if(getContinuousIntersections(colliding, collisions))
		{       
                        double t_ivs_start = omp_get_wtime();
                        cout<<"computing volumes\n";
			IntersectionVolumes ivs;
			for (vector<Collisions>::iterator cItr=collisions.begin(); cItr!=collisions.end(); ++cItr)
				ivs.push_back(IntersectionVolume(*getSlaveMesh(), m_thickness, *cItr, m_periodicLength, m_nv));
                        cout<<"after computing volumes\n";
                        double t_ivs_end = omp_get_wtime();
                        cout<<"total time for ivs: "<<(t_ivs_end - t_ivs_start)<<endl;
                        
                        IV.resize(ivs.size(), 0); 
                        IS_BDV.resize(ivs.size(), 1.0); 
                        IV_size = ivs.size();
                        int ivs_count = 0;
                        double t_sumgrad_start = omp_get_wtime();
                        for (IntersectionVolumesIterator ivItr=ivs.begin(); ivItr!=ivs.end(); ++ivItr)
                        {
                          IV[ivs_count] = ivItr->getVolume();
                          if(ivItr->m_bdv)
                              IS_BDV[ivs_count]=0;
                          ivs_count++;
                                                    
                          gradient = gradient + ivItr->getGradients();
                          for (Eigen::SparseVector<double>::InnerIterator grad_it(ivItr->getGradients()); grad_it; ++grad_it)
                          {
                            gradient_index.insert(grad_it.index()) = ivs_count;
                          }
                        }	
                        double t_sumgrad_end = omp_get_wtime();
                        cout<<"total time for sumgrad: "<<(t_sumgrad_end - t_sumgrad_start)<<endl;
		}
	}
        else
        {
            assert(false);
        }
}

void ColSurface::controlMeshUpdate(std::set<unsigned int> &selectedVertices, map<unsigned int, unsigned int> &colliding, 
    bool handleCollisions, bool handleElasticity, Vector_3 &translate, bool rigid, bool moveAll, int limiter, 
    Eigen::SparseVector<double> &gradient, std::vector<double> &IV, Eigen::SparseVector<int> &gradient_index, int &IV_size, 
    std::vector<double> &IS_BDV)
{
	controlMeshUpdateIteration(selectedVertices, colliding, 
            handleCollisions, handleElasticity, translate, rigid, moveAll, limiter,
            gradient, IV, gradient_index, IV_size, IS_BDV);
}

bool ColSurface::intersection(const Polyhedron &P, Triangles_IS &triangles_is)
{
    std::vector<Box> boxes;
    boxes.reserve( P.size_of_facets());
    for ( Facet_const_iterator i = P.facets_begin(); i != P.facets_end(); ++i){
        boxes.push_back(
            Box( i->halfedge()->vertex()->point().bbox()
               + i->halfedge()->next()->vertex()->point().bbox()
               + i->halfedge()->next()->next()->vertex()->point().bbox(),
                 i));
    }
    std::vector<const Box*> box_ptr;
    box_ptr.reserve( P.size_of_facets());
    for ( std::vector<Box>::iterator j = boxes.begin(); j != boxes.end(); ++j){
        box_ptr.push_back( &*j);
    }
    CGAL::box_self_intersection_d( box_ptr.begin(), box_ptr.end(),
                                   Intersect_facets(triangles_is), std::ptrdiff_t(2000));
    if(triangles_is.size()>0)
      return true;
    else
      return false;
      
}

bool ColSurface::getContinuousIntersections(map<unsigned int, unsigned int> &colliding, vector<Collisions> &collisions)
{
        // INIT BOUNDING BOX
	vector<Vector3d> temp;
	temp.reserve(m_slaveMesh->size_of_vertices());
	Vertex_iterator vit = m_slaveMesh->vertices_begin();
	while (vit != m_slaveMesh->vertices_end())
	{
		temp.push_back(vit->P());
		++vit;
	}

	AABB aabb(temp);
	vector<double> sizes;
	sizes.push_back(aabb.getSize()[0]);
	sizes.push_back(aabb.getSize()[1]);
	sizes.push_back(aabb.getSize()[2]);
	double m = *max_element(sizes.begin(), sizes.end());
	m *= 1.1;
                
        m = max(m, m_periodicLength);

        if(m_aabb)
            delete m_aabb;
	m_aabb = new AABB(aabb.getCenter(),Vector3d(m,m,m));
	// END BOUNDING BOX

	Collisions colls;
	vector<Collision> collsEE;
	
        double t_cvwh_start = omp_get_wtime();
        cout<<"Begin collision detection "<<endl;
	PrismVolume::computeVolumeWithHash(*m_slaveMesh, *m_aabb, m_thickness, m_edgeLength, m_subMeshes, colls, collsEE, 
                hash_grid, m_nv_resident, m_periodicLength);
        cout<<"End collision detection "<<endl;
        double t_cvwh_end = omp_get_wtime();
        cout<<"total time for cvwh: "<<(t_cvwh_end - t_cvwh_start)<<endl;
#ifdef DEBUGVG
        ////cout<<"ee collision size "<<collsEE.size()<<endl;
        ////cout<<"fv collision size "<<colls.size()<<endl;
#endif
	
        if (colls.empty() && collsEE.empty())
	{
		collisions.clear();
                cout<<"No group, col size:"<<collisions.size()<<endl;
		return false;
	}
        
        // merge fv and ee collisions
        vector<Collision> colls_total;
        colls_total.reserve(colls.size()+collsEE.size());
        colls_total.insert(colls_total.end(),colls.begin(),colls.end());
        colls_total.insert(colls_total.end(),collsEE.begin(),collsEE.end());
        cout<<"Begin group "<<endl;
	groupCollisions(colls_total, colliding, collisions);
        cout<<"End group, col size:"<<collisions.size()<<endl;

	return (collisions.size() > 0);
}

// Groups a list of collisions into disjoint sets
//
void ColSurface::groupCollisions(Collisions &inputColls, map<unsigned int, unsigned int> &colliding, vector<Collisions> &collisions)
{
	collisions.clear();

	// First put every collision in its own list
	//
	std::list<Collisions> collisionsLists;
	std::list<boost::unordered_set<unsigned int> > vertexLists;
	for (CollisionsIterator cItr=inputColls.begin(); cItr!=inputColls.end(); ++cItr)
	{
		Collisions c;
		c.push_back(*cItr);
		collisionsLists.push_back(c);

		colliding[cItr->getVertex()] = 60;
		colliding[cItr->getTriangleVertex1()] = 60;
		colliding[cItr->getTriangleVertex2()] = 60;
		colliding[cItr->getTriangleVertex3()] = 60;

		boost::unordered_set<unsigned int> s;
		s.insert(cItr->getVertex());
		s.insert(cItr->getTriangleVertex1());
		s.insert(cItr->getTriangleVertex2());
		s.insert(cItr->getTriangleVertex3());

		unsigned int verts[4];
		verts[0] = cItr->getVertex();
		verts[1] = cItr->getTriangleVertex1();
		verts[2] = cItr->getTriangleVertex2();
		verts[3] = cItr->getTriangleVertex3();

		// Insert one-ring of each vertex into the lists
		// to look for overlap
		//
		for (size_t i=0; i<4; ++i)
		{
			Halfedge_Vertex_const_circulator hvcir = m_slaveMesh->getVertex(verts[i])->vertex_begin();
			do
			{
				s.insert(hvcir->opposite()->vertex()->getIndex());
				hvcir++;
			} while (hvcir != m_slaveMesh->getVertex(verts[i])->vertex_begin());
		}
		vertexLists.push_back(s);
	}
	

	// Merge collisions with common vertices
	//
	while (!collisionsLists.empty())
	{
		// Start with first collision "zone"
		//
		Collisions c = collisionsLists.front();
		collisionsLists.pop_front();
		boost::unordered_set<unsigned int> si = vertexLists.front();
		vertexLists.pop_front();

		int s = collisionsLists.size();
		bool hits = false;

		// Loop through other zones
		//
		for (int i=0; i<s; ++i)
		{
			Collisions n = collisionsLists.front();
			collisionsLists.pop_front();
			boost::unordered_set<unsigned int> next = vertexLists.front();
			vertexLists.pop_front();

			// Search for common elements
			//
			bool common = false;
			for (boost::unordered_set<unsigned int>::iterator sItr=next.begin(); sItr!=next.end(); ++sItr)
			{
				if (si.find(*sItr) != si.end())
				{
					// They share a vertex
					//
					common = true;
					break;
				}
			}

			// If there are elements in common, then take the union of the two sets
			//
			if (common)
			{
				// Merge
				//
				si.insert(next.begin(), next.end());
				c.insert(c.end(), n.begin(), n.end());

				// We need to keep track of new zones for later, so signify that there's a new one
				//
				hits = true;
			}
			else
			{
				// Otherwise put it back on the queue
				//
				vertexLists.push_back(next);
				collisionsLists.push_back(n);
			}
		}

		// If the search returned no hits, then put this set on the final list
		//
		if (!hits)
			collisions.push_back(c);
		else
		{
			collisionsLists.push_back(c);
			vertexLists.push_back(si);
		}
	}
}

void ColSurface::computeNormals(Polyhedron* p)
{
	Polyhedron& P = *p;
	
	for (Vertex_iterator vit = P.vertices_begin(); vit != P.vertices_end(); vit++)
		vit->normal = Vector_3(0,0,0);
    
        for (Facet_iterator fit = P.facets_begin(); fit != P.facets_end(); fit++) 
        {
		Vertex* v0 = P.v_array[fit->vi[0]];
		Vertex* v1 = P.v_array[fit->vi[1]];
		Vertex* v2 = P.v_array[fit->vi[2]];
        
		fit->normal = Vector_3(0,0,0);
                fit->normal = CGAL::cross_product((v1->Pcgal()-v0->Pcgal()),(v2->Pcgal()-v1->Pcgal()));
		fit->normal = fit->normal / sqrt(fit->normal.squared_length());
                
                v0->normal = v0->normal + fit->normal;
		v1->normal = v1->normal + fit->normal;
		v2->normal = v2->normal + fit->normal;
        }
	
	for (Vertex_iterator vit = P.vertices_begin(); vit != P.vertices_end(); vit++)
		vit->normal = vit->normal / sqrt(vit->normal.squared_length());
}

double ColSurface::minDistance()
{
  PQP_REAL R1[3][3], R2[3][3], T1[3], T2[3];
  
  R1[0][0] = R1[1][1] = R1[2][2] = 1.0;
  R1[0][1] = R1[1][0] = R1[2][0] = 0.0;
  R1[0][2] = R1[1][2] = R1[2][1] = 0.0;

  R2[0][0] = R2[1][1] = R2[2][2] = 1.0;
  R2[0][1] = R2[1][0] = R2[2][0] = 0.0;
  R2[0][2] = R2[1][2] = R2[2][1] = 0.0;
  
  T1[0] = 0.0;  T1[1] = 0.0; T1[2] = 0.0;
  T2[0] = 0.0;  T2[1] = 0.0; T2[2] = 0.0;

  double minD = 1e50;
  int count = m_subMeshes.size();
  //#pragma omp parallel for collapse(2)
  for(int i=0; i<count; ++i)
  {
    std::set<unsigned int> subMesh1 = m_subMeshes[i];
    PQP_Model *b1 = new PQP_Model;
    loadModel(subMesh1,b1);
    //#pragma omp parallel for
    for(int j=i+1; j<count; ++j)
    {
      std::set<unsigned int> subMesh2 = m_subMeshes[j];
      PQP_Model *b2 = new PQP_Model;
      loadModel(subMesh2,b2);

      PQP_DistanceResult dres;
      PQP_Distance(&dres, R1, T1, b1, R2, T2, b2, 0.0, 0.0);

      //#pragma omp critical(dataupdate1)
      {
        if(dres.Distance()<minD)
          minD = dres.Distance();
      }


      delete b2;
    }
    delete b1;
  }

  return minD;
}

double ColSurface::minDistanceBD(bool flag_orig)
{
  PQP_REAL R1[3][3], R2[3][3], T1[3], T2[3];
  
  R1[0][0] = R1[1][1] = R1[2][2] = 1.0;
  R1[0][1] = R1[1][0] = R1[2][0] = 0.0;
  R1[0][2] = R1[1][2] = R1[2][1] = 0.0;

  R2[0][0] = R2[1][1] = R2[2][2] = 1.0;
  R2[0][1] = R2[1][0] = R2[2][0] = 0.0;
  R2[0][2] = R2[1][2] = R2[2][1] = 0.0;
  
  T1[0] = 0.0;  T1[1] = 0.0; T1[2] = 0.0;
  T2[0] = 0.0;  T2[1] = 0.0; T2[2] = 0.0;

  double minD = 1e50;
  int count = m_subMeshes.size();
  //std::cout<<"count: "<<count<<std::endl;
  //std::cout<<"nv: "<<m_nv<<std::endl;
  //#pragma omp parallel for collapse(2)
  #pragma omp parallel for
  for(int i=0; i<count; ++i)
  {
    std::set<unsigned int> subMesh1 = m_subMeshes[i];
    PQP_Model *b1 = new PQP_Model;
    loadModel(subMesh1,b1,flag_orig);
    //#pragma omp parallel for
    for(int j=(i>=m_nv)?0:i+1; j<m_nv; ++j)
    {
      //std::cout<<"pair: ("<<i<<", "<<j<<")"<<std::endl;
      if(!subMeshCheck(i,j,flag_orig))
          continue;
      //std::cout<<"checking pair: ("<<i<<", "<<j<<")"<<std::endl;
      std::set<unsigned int> subMesh2 = m_subMeshes[j];
      PQP_Model *b2 = new PQP_Model;
      loadModel(subMesh2,b2,flag_orig);

      PQP_DistanceResult dres;
      PQP_Distance(&dres, R1, T1, b1, R2, T2, b2, 0.0, 0.0);

      #pragma omp critical(dataupdate1)
      {
        if(dres.Distance()<minD)
          minD = dres.Distance();
      }


      delete b2;
    }
    delete b1;
  }

  return minD;
}

void ColSurface::loadModel(std::set<unsigned int> subMesh, PQP_Model *b, bool flag_orig)
{ 
  int count = 0;
  b->BeginModel();
  //#pragma omp parallel for
  for (set<unsigned int>::iterator sItr=subMesh.begin(); sItr!=subMesh.end(); ++sItr)
  {
    Facet* f = m_slaveMesh->f_array[*sItr];
    Vector3d v0,v1,v2;
    if(flag_orig)
    {
        v0 = m_slaveMesh->getVertex(f->vi[0])->Pori();
        v1 = m_slaveMesh->getVertex(f->vi[1])->Pori();
        v2 = m_slaveMesh->getVertex(f->vi[2])->Pori();
    }
    else
    {
        v0 = m_slaveMesh->getVertex(f->vi[0])->P();
        v1 = m_slaveMesh->getVertex(f->vi[1])->P();
        v2 = m_slaveMesh->getVertex(f->vi[2])->P();
    }

    PQP_REAL p0[3], p1[3], p2[3];
    p0[0] = v0[0];
    p0[1] = v0[1];
    p0[2] = v0[2];
    
    p1[0] = v1[0];
    p1[1] = v1[1];
    p1[2] = v1[2];

    p2[0] = v2[0];
    p2[1] = v2[1];
    p2[2] = v2[2];

    //#pragma omp critical(dataupdate2)
    //{
      b->AddTri(p0,p1,p2,count);
      count += 1;
    //}
  }
  b->EndModel();
}

bool ColSurface::subMeshCheck(int i, int j, bool flag_orig)
{
    Vector3d v[3];
    Vector3d mini, maxi, minj, maxj;

    if(flag_orig)
    {
        mini = m_slaveMesh->getVertex( m_slaveMesh->f_array[*m_subMeshes[i].begin()]->vi[0] )->Pori();
        maxi = m_slaveMesh->getVertex( m_slaveMesh->f_array[*m_subMeshes[i].begin()]->vi[0] )->Pori();
        
        minj = m_slaveMesh->getVertex( m_slaveMesh->f_array[*m_subMeshes[j].begin()]->vi[0] )->Pori();
        maxj = m_slaveMesh->getVertex( m_slaveMesh->f_array[*m_subMeshes[j].begin()]->vi[0] )->Pori();
    }
    else
    {
        mini = m_slaveMesh->getVertex( m_slaveMesh->f_array[*m_subMeshes[i].begin()]->vi[0] )->P();
        maxi = m_slaveMesh->getVertex( m_slaveMesh->f_array[*m_subMeshes[i].begin()]->vi[0] )->P();
        
        minj = m_slaveMesh->getVertex( m_slaveMesh->f_array[*m_subMeshes[j].begin()]->vi[0] )->P();
        maxj = m_slaveMesh->getVertex( m_slaveMesh->f_array[*m_subMeshes[j].begin()]->vi[0] )->P();
    }

    for (set<unsigned int>::iterator sItr=m_subMeshes[i].begin(); sItr!=m_subMeshes[i].end(); ++sItr)
    {
        Facet* f = m_slaveMesh->f_array[*sItr];
        if(flag_orig)
        {
            for(int ii=0; ii<3; ii++)
                v[ii] = m_slaveMesh->getVertex(f->vi[ii])->Pori();
        }
        else
        {
            for(int ii=0; ii<3; ii++)
                v[ii] = m_slaveMesh->getVertex(f->vi[ii])->P();
        }

        for(int ii=0; ii<3; ii++)
        {
            for(int jj=0; jj<3; jj++)
            {
                if(mini[jj] > v[ii][jj])
                    mini[jj] = v[ii][jj];
                if(maxi[jj] < v[ii][jj])
                    maxi[jj] = v[ii][jj];
            }
        }
    }
    
    for (set<unsigned int>::iterator sItr=m_subMeshes[j].begin(); sItr!=m_subMeshes[j].end(); ++sItr)
    {
        Facet* f = m_slaveMesh->f_array[*sItr];
        if(flag_orig)
        {
            for(int ii=0; ii<3; ii++)
                v[ii] = m_slaveMesh->getVertex(f->vi[ii])->Pori();
        }
        else
        {
            for(int ii=0; ii<3; ii++)
                v[ii] = m_slaveMesh->getVertex(f->vi[ii])->P();
        }

        for(int ii=0; ii<3; ii++)
        {
            for(int jj=0; jj<3; jj++)
            {
                if(minj[jj] > v[ii][jj])
                    minj[jj] = v[ii][jj];
                if(maxj[jj] < v[ii][jj])
                    maxj[jj] = v[ii][jj];
            }
        }
    }

    for (int jj=0; jj<3; ++jj)
    {
	mini[jj] -= m_thickness;
        maxi[jj] += m_thickness;

        minj[jj] -= m_thickness;
        maxj[jj] += m_thickness;
    }

        
    bool nointersect = false;
    nointersect = nointersect || (maxi[0] < minj[0]) || (mini[0] > maxj[0]);
    nointersect = nointersect || (maxi[1] < minj[1]) || (mini[1] > maxj[1]);
    nointersect = nointersect || (maxi[2] < minj[2]) || (mini[2] > maxj[2]);
    return !nointersect;
}

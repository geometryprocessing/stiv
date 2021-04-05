#include "Structures.h"
#include "Common.h"
#include <boost/foreach.hpp>

ExtendedMesh::ExtendedMesh() :
    visible(true),name(""),filename(""),active(false) {
}

void ExtendedMesh::reset() {
    visible = true;
    Vertex_iterator vit = vertices_begin();
    while (vit != vertices_end()) {
        vit->reset();
        vit++;
    }
	
    Face_iterator fit = facets_begin();
    while (fit != facets_end()) {
        fit->reset();
        fit++;
    }
    
	init();
}

void ExtendedMesh::init() {
    int count = 0;
    v_array.clear();
    for (Vertex_iterator i = this->vertices_begin();
            i != this->vertices_end(); i++) {
        i->setIndex(count);
        v_array.push_back(&*i);
        count++;
    }

    count = 0;
    f_array.clear();
    for (Facet_iterator i = this->facets_begin();
            i != this->facets_end(); i++) {
        i->setIndex(count);
        f_array.push_back(&*i);
        count++;
    }

	// generate edges and tag them
	
	// mark all edges as non-marked
	for (Halfedge_iterator i = this->halfedges_begin();i!= this->halfedges_end(); ++i)
	{
		i->marked = false;
	}
	count = 0;
    e_array.clear();
	for (Halfedge_iterator i = this->halfedges_begin();i!= this->halfedges_end(); ++i)
	{
		if (!i->marked)
		{
			i->marked = true;
			i->opposite()->marked = true;
			i->setIndex(count);
			i->opposite()->setIndex(count);
			
			e_array.push_back(&*i);
			count++;
		}
	}

	size_t i=0;
	for (Edge_iterator eItr=edges_begin(); eItr!=edges_end(); ++eItr)
		count++;

	////std::cout << e_array.size() << endl;
	////std::cout << (size_of_halfedges() / 2) << endl;
	////std::cout << count << endl;
	
    for (Facet_iterator i = this->facets_begin(); i != this->facets_end(); i++) 
	{
		Halfedge_around_facet_const_circulator hfc = i->facet_begin();
		i->vi[0] = hfc->vertex()->getIndex();
		++hfc;
		i->vi[1] = hfc->vertex()->getIndex();
		++hfc;
		i->vi[2] = hfc->vertex()->getIndex();
    }

//	m_pos = Eigen::VectorXd::Zero(this->size_of_vertices()*3);
//	m_oripos = Eigen::VectorXd::Zero(this->size_of_vertices()*3);
	count = 0;
	for (Vertex_iterator i = this->vertices_begin(); i != this->vertices_end(); i++) 
	{
		i->P() = Eigen::Vector3d(i->point()[0],i->point()[1],i->point()[2]);
		i->Pori() = i->P();
//		m_oripos.segment(count*3,3) = Eigen::Vector3d(i->point()[0],i->point()[1],i->point()[2]);
//		i->pos = &m_pos;
//		i->oripos = &m_oripos;
		++count;
    }
	
//	updateWeight(); // TODO
}

void ExtendedMesh::updateAllCotanPori()
{
	// compute weights
	Halfedge_iterator hit = halfedges_begin();
    while (hit != halfedges_end()) 
	{
        hit->weight = cotangentWeightPori(&(*hit));
        hit++;
    }
}

void ExtendedMesh::updateWeight()
{
	// compute weights
	Halfedge_iterator hit = halfedges_begin();
    while (hit != halfedges_end()) 
	{
        hit->weight = computeWeight(&(*hit));
        hit++;
    }
	
	// precompute 1-rings
    for (Vertex_iterator i = this->vertices_begin(); i != this->vertices_end(); i++)  
	{
		i->ns.clear();
		i->nsW.clear();
		
		typedef std::pair<int, Halfedge*> PairID;
		std::vector<PairID > vp = get1RingH(i->getIndex());
		BOOST_FOREACH(PairID p, vp)
		{
			i->ns.push_back(p.first);
			i->nsW.push_back(getWeightVs(p.second));
		}
    }
}

void ExtendedMesh::initAABB(double h)
{
	// AABBs for faces
	int size = this->size_of_facets();

	//double *maxSize = new double[omp_get_num_threads()];
	//for (size_t j=0; j<omp_get_num_threads(); ++j)
	//	maxSize[j] = 0.0;

	#pragma omp parallel for
        for (int i = 0; i < size; i++)
	{
		Facet* fit = this->f_array[i];
		Point3d p[6];
		p[0] = v_array[fit->vi[0]]->Pori();
		p[1] = v_array[fit->vi[1]]->Pori();
		p[2] = v_array[fit->vi[2]]->Pori();
		p[3] = v_array[fit->vi[0]]->P();
		p[4] = v_array[fit->vi[1]]->P();
		p[5] = v_array[fit->vi[2]]->P();
		
		for(unsigned j = 0; j < 3; ++j) 
		{
			fit->min[j] = p[0][j];
			fit->max[j] = p[0][j];
		}

		for(unsigned j = 1; j < 6; ++j)
		{
			for(unsigned k = 0; k < 3; ++k) 
			{
				if (p[j][k] < fit->min[k])
					fit->min[k] = p[j][k];
				if (p[j][k] > fit->max[k])
					fit->max[k] = p[j][k];
			}
		}

		//int threadID = omp_get_thread_num();
		for (unsigned j=0; j<3; ++j)
		{
			fit->min[j] -= h; fit->min[j] -= 1e-8;
			fit->max[j] += h; fit->max[j] += 1e-8;
		
			//if (fabs(fit->max[j] - fit->min[j]) > maxSize[threadID])
			//	maxSize[threadID] = fabs(fit->max[j] - fit->min[j]);
		}
        }

	// AABBs for vertices
	int sizev = this->size_of_vertices();
	
        #pragma omp parallel for
        for (int i = 0; i < sizev; i++)
	{
		Vertex* vit = this->v_array[i];
		Point3d p[2];
		p[0] = vit->Pori();
		p[1] = vit->P();
		
		for(unsigned j = 0; j < 3; ++j) 
		{
			vit->min[j] = p[0][j];
			vit->max[j] = p[0][j];
		}
		
		for(unsigned j = 0; j < 3; ++j) 
		{
			if (p[1][j] < vit->min[j])
				vit->min[j] = p[1][j];
			if (p[1][j] > vit->max[j])
				vit->max[j] = p[1][j];
		}
		
		for (unsigned j=0; j<3; ++j)
		{
			vit->min[j] -= h/2; vit->min[j] -= 1e-8;
			vit->max[j] += h/2; vit->min[j] += 1e-8;
			
		}
        }

#ifdef DETECTEE
	// AABBs for edges
	size = this->e_array.size();
	//cerr << "Size of edge array: " << size << endl;
	
        #pragma omp parallel for
        for (int i = 0; i < size; i++)
	{
		Halfedge* eit = this->e_array[i];
		
		Point3d p[4];
		p[0] = eit->vertex()->Pori();
		p[1] = eit->opposite()->vertex()->Pori();
		
		p[2] = eit->vertex()->P();
		p[3] = eit->opposite()->vertex()->P();

//		cerr << "-----> BOX for edge (" << eit->vertex()->getIndex() << "," << eit->opposite()->vertex()->getIndex() << ")" << endl;
//		
//		cerr << p[0].transpose() << " | \n" << p[1].transpose() << " | \n" << p[2].transpose() << " | \n" << p[3].transpose() << endl;
		
		for(unsigned j = 0; j < 3; ++j) 
		{
			eit->min[j] = p[0][j];
			eit->max[j] = p[0][j];
		}
		
		for(unsigned j = 1; j < 4; ++j)
		{
			for(unsigned k = 0; k < 3; ++k) 
			{
				if (p[j][k] < eit->min[k])
					eit->min[k] = p[j][k];
				if (p[j][k] > eit->max[k])
					eit->max[k] = p[j][k];
			}
		}
		
		for (unsigned j=0; j<3; ++j)
		{
			//eit->min[j] -= h;
			//eit->max[j] += h;
                        eit->min[j] -= h/2; eit->min[j] -= 1e-8; 
			eit->max[j] += h/2; eit->max[j] += 1e-8;
		}
		
//		cerr << "BOX min: " << eit->min[0] << " " << eit->min[1] << " " << eit->min[2] <<endl;
//		cerr << "BOX MAX: " << eit->max[0] << " " << eit->max[1] << " " << eit->max[2] <<endl;
		
		// copy to the opposite he
		for (unsigned j=0; j<3; ++j)
		{
			eit->opposite()->min[j] = eit->min[j];
			eit->opposite()->max[j] = eit->max[j];
		}
        }
#endif
}

double ExtendedMesh::getArea(unsigned int idx)
{
	double A = 0.0;
	Halfedge_Vertex_circulator hvcir1 = getVertex(idx)->vertex_begin();
	Halfedge_Vertex_circulator hvcir2 = getVertex(idx)->vertex_begin(); ++hvcir2;
	do
	{
		A += (hvcir1->opposite()->vertex()->Pori() - getVertex(idx)->Pori()).cross(
			  hvcir2->opposite()->vertex()->Pori() - getVertex(idx)->Pori()).norm();
		++hvcir1;
		++hvcir2;
	} while (hvcir1 != getVertex(idx)->vertex_begin());
	A = A / 6.0;

	return A;
}

std::vector<std::pair<int, double> > ExtendedMesh::get1Ring(int i)
{
	assert(i < v_array.size());
	Vertex* v = v_array[i];
	
	vector<pair<int, double> > ret;
	
	Halfedge_Vertex_circulator hvcit = v->vertex_begin();
	do 
	{
		pair<int, double> pair;
		Vertex_handle neighbor_v = hvcit->opposite()->vertex();
		pair.first  = neighbor_v->getIndex();
		pair.second = hvcit->weight;
		
		ret.push_back(pair);
		hvcit++;
	} while (hvcit != v->vertex_begin());
	return ret;
}

std::vector<std::pair<int, Halfedge*> > ExtendedMesh::get1RingH(int i)
{
	assert(i < v_array.size());
	Vertex* v = v_array[i];
	
	vector<pair<int, Halfedge*> > ret;
	
	Halfedge_Vertex_circulator hvcit = v->vertex_begin();
	do 
	{
		pair<int, Halfedge*> pa;
		Vertex_handle neighbor_v = hvcit->opposite()->vertex();
		pa.first  = neighbor_v->getIndex();
		pa.second = &*hvcit;
		
		ret.push_back(pa);
		hvcit++;
	} while (hvcit != v->vertex_begin());
	return ret;
}

double ExtendedMesh::computeWeight(Halfedge* h) 
{
    return cotangentWeight(h);
	//    return uniformWeight(h);
}

double ExtendedMesh::uniformWeight(Halfedge* h) const 
{
    return 1.0;
}

double ExtendedMesh::cotangentWeight(Halfedge* h) const 
{
    //     p1
    //    /^\
    //  p2 | p3
    //    \|/
    //     p0
    Point3d p0 = h->opposite()->vertex()->P();
    Point3d p1 = h->vertex()->P();
    Point3d p2 = h->next()->vertex()->P();
    Point3d p3 = h->opposite()->next()->vertex()->P();
    double cot_left  = 0.0;
    double cot_right = 0.0;
    if (!h->is_border())
        cot_left =  max(0.0, cotangent(p0, p2, p1));
    if (!h->opposite()->is_border())
        cot_right = max(0.0, cotangent(p0, p3, p1));
    if (cot_left < 0 || cot_right < 0) {
        h->vertex()->marked = true;
        h->opposite()->vertex()->marked = true;
    }
    return 0.5 * (cot_left + cot_right);
}

double ExtendedMesh::cotangentWeightPori(Halfedge* h) const 
{
    //     p1
    //    /^\
    //  p2 | p3
    //    \|/
    //     p0
    Point3d p0 = h->opposite()->vertex()->Pori();
    Point3d p1 = h->vertex()->Pori();
    Point3d p2 = h->next()->vertex()->Pori();
    Point3d p3 = h->opposite()->next()->vertex()->Pori();
    double cot_left  = 0.0;
    double cot_right = 0.0;
    if (!h->is_border())
        cot_left =  max(0.0, cotangent(p0, p2, p1));
    if (!h->opposite()->is_border())
        cot_right = max(0.0, cotangent(p0, p3, p1));
    if (cot_left < 0 || cot_right < 0) {
        h->vertex()->marked = true;
        h->opposite()->vertex()->marked = true;
    }
    return 0.5 * (cot_left + cot_right);
}

/// Return cotangent of (P,Q,R) corner (i.e. cotan of QP,QR angle).
double ExtendedMesh::cotangent(const Point3d& P,
								   const Point3d& Q,
								   const Point3d& R) const
{
    Point3d u = P - Q;
    Point3d v = R - Q;
    // (u . v)/((u x v).len)
    double dot = (u.dot(v));
    Point3d cross_vector = u.cross(v);
    double cross_norm = cross_vector.norm();
    assert(cross_norm != 0.0);
    if(cross_norm > EPSILON)
        return (dot/cross_norm);
    else
        return 0.0; // undefined
}

WeightInfo ExtendedMesh::getWeightVs(Halfedge* h)
{
	WeightInfo wi;
	//     p1
    //    /^\
    //  p2 | p3
    //    \|/
    //     p0
    wi.vis[0] = h->opposite()->vertex()->getIndex();
    wi.vis[1] = h->vertex()->getIndex();
    wi.vis[2] = h->next()->vertex()->getIndex();
    wi.vis[3] = h->opposite()->next()->vertex()->getIndex();
    
	wi.is_border = h->is_border();
	wi.opposite_is_border = h->opposite()->is_border();
	
	wi.w = computeWeight(h);
	
	// DEBUG
	
	vector<Eigen::Vector3d> ps;
	
	ps.push_back(h->opposite()->vertex()->P());
	ps.push_back(h->vertex()->P());
	ps.push_back(h->next()->vertex()->P());
	ps.push_back(h->opposite()->next()->vertex()->P());

	assert(computeWeightVs(ps, wi) == computeWeight(h));
	
	return wi;
}

double ExtendedMesh::computeWeightVs(std::vector<Eigen::Vector3d>& ps, WeightInfo& wi)
{
	//     p1
    //    /^\
    //  p2 | p3
    //    \|/
    //     p0
	assert(ps.size() == 4);
    Point3d p0 = ps[0];
    Point3d p1 = ps[1];
    Point3d p2 = ps[2];
    Point3d p3 = ps[3];
    double cot_left  = 0.0;
    double cot_right = 0.0;
    if (!wi.is_border)
        cot_left =  max(0.0, cotangent(p0, p2, p1));
    if (!wi.opposite_is_border)
        cot_right = max(0.0, cotangent(p0, p3, p1));
//    if (cot_left < 0 || cot_right < 0) {
//        h->vertex()->marked = true;
//        h->opposite()->vertex()->marked = true;
//    }
    return 0.5 * (cot_left + cot_right);
	
}

void ExtendedMesh::copyCgalToEigen()
{
	for (Vertex_iterator i = this->vertices_begin(); i != this->vertices_end(); i++) 
	{
		i->P() = Eigen::Vector3d(i->point()[0],i->point()[1],i->point()[2]);
    }
}

void ExtendedMesh::copyEigenToCgal()
{
	for (Vertex_iterator i = this->vertices_begin(); i != this->vertices_end(); i++) 
	{
		i->point() = Point_3(i->P()[0],i->P()[1],i->P()[2]);
    }	
}


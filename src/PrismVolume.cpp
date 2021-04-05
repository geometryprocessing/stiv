#include "PrismVolume.h"

#include "rpoly.h"
#include <boost/foreach.hpp>
#include "Hash.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "MeshManager.h"

using std::vector;
using namespace Eigen;



void PrismVolume::validateCollisionAndAdd(Facet* fit, Vertex* vit, double h, vector<Collision>& collisions, Polyhedron *p)
{
	//Polyhedron* p = MeshManager::getInstance()->getMesh()->getSlaveMesh();

	Vertex* vh0 = p->v_array[fit->vi[0]];
	Vertex* vh1 = p->v_array[fit->vi[1]];
	Vertex* vh2 = p->v_array[fit->vi[2]];
		
	double tRet;
	if (validateCollision(vit, vh0->Pori(), vh1->Pori(), vh2->Pori(), vh0->P(), vh1->P(), vh2->P(), h, tRet))
	{
		collisions.push_back(Collision(vit->getIndex(),
									   vh0->getIndex(),
									   vh1->getIndex(),
									   vh2->getIndex(),
									   tRet, true));
		collisions.back().setFace(fit->getIndex());
	}
}

bool PrismVolume::validateCollisionAndReplace(Facet* fit, Vertex* vit, double h, Collision& c, Polyhedron *p, double periodicLength)
{
    if(periodicLength > 0)
    {
        //Polyhedron* p = MeshManager::getInstance()->getMesh()->getSlaveMesh();
	
	Vertex* vh0 = p->v_array[fit->vi[0]];
	Vertex* vh1 = p->v_array[fit->vi[1]];
	Vertex* vh2 = p->v_array[fit->vi[2]];

        Vector3d transV(0,0,0);
        getTranslationVector(transV, vh0, vit, periodicLength);
	
	double tRet;
        Vector3d v0pori = vh0->Pori()+transV;
        Vector3d v1pori = vh1->Pori()+transV;
        Vector3d v2pori = vh2->Pori()+transV;
        
        Vector3d v0p = vh0->P()+transV;
        Vector3d v1p = vh1->P()+transV;
        Vector3d v2p = vh2->P()+transV;

	if (validateCollision(vit, v0pori, v1pori, v2pori,
                                   v0p,    v1p,    v2p, h, tRet))
	{
		c = Collision(vit->getIndex(),
			   	      vh0->getIndex(),
					  vh1->getIndex(),
					  vh2->getIndex(),
					  tRet, true);
		c.setFace(fit->getIndex());
		
		return true;
	}
	
	return false;
    }
    else
    {
	//Polyhedron* p = MeshManager::getInstance()->getMesh()->getSlaveMesh();
	
	Vertex* vh0 = p->v_array[fit->vi[0]];
	Vertex* vh1 = p->v_array[fit->vi[1]];
	Vertex* vh2 = p->v_array[fit->vi[2]];
	
	double tRet;
	if (validateCollision(vit, vh0->Pori(), vh1->Pori(), vh2->Pori(), vh0->P(), vh1->P(), vh2->P(), h, tRet))
	{
		c = Collision(vit->getIndex(),
			   	      vh0->getIndex(),
					  vh1->getIndex(),
					  vh2->getIndex(),
					  tRet, true);
		c.setFace(fit->getIndex());
		
		return true;
	}
	
	return false;
    }
}

// Test if segment specified by points p0 and p1 intersects AABB b
bool PrismVolume::TestSegmentAABB(Vector3d &p0, Vector3d &p1, Vector3d &min, Vector3d &max)
{
	Vector3d c = (min + max) * 0.5f; // Box center-point
	Vector3d e = max - c;             // Box halflength extents
	Vector3d m = (p0 + p1) * 0.5f;       // Segment midpoint
	Vector3d d = p1 - m;                // Segment halflength vector
	m = m - c;                        // Translate box and segment to origin
	
	// Try world coordinate axes as separating axes
	double adx = fabs(d[0]);
	if (fabs(m[0]) > e[0] + adx) return false;
	double ady = fabs(d[1]);
	if (fabs(m[1]) > e[1] + ady) return false;
	double adz = fabs(d[2]);
	if (fabs(m[2]) > e[2] + adz) return false;
	
	// Add in an epsilon term to counteract arithmetic errors when segment is 
	// (near) parallel to a coordinate axis (see text for detail) 
	adx += EPSILON; ady += EPSILON; adz += EPSILON;
	
	// Try cross products of segment direction vector with coordinate axes
	if (fabs(m[1] * d[2] - m[2] * d[1]) > e[1] * adz + e[2] * ady) return false;
	if (fabs(m[2] * d[0] - m[0] * d[2]) > e[0] * adz + e[2] * adx) return false;
	if (fabs(m[0] * d[1] - m[1] * d[0]) > e[0] * ady + e[1] * adx) return false;
	
	// No separating axis found; segment must be overlapping AABB 
	return true;
}

bool PrismVolume::validateCollisionFastFast(Facet* fit, Vertex* vit)
{
	int vi = vit->getIndex();
	int f0 = fit->vi[0];
	int f1 = fit->vi[1];
	int f2 = fit->vi[2];
	
	if ((vi == f0) || (vi == f1) || (vi == f2))
		return false;
	
	return true;
}

void PrismVolume::getTranslationVector(Vector3d &transV, Vertex* vit1, Vertex* vit2, double periodicLength)
{
        double dx = vit1->P()[0] - vit2->P()[0];
        double dy = vit1->P()[1] - vit2->P()[1];
        double dz = vit1->P()[2] - vit2->P()[2];

        if(dx >  0.5*periodicLength)
            transV[0] = -periodicLength;
        if(dx < -0.5*periodicLength)
            transV[0] =  periodicLength;

        if(dy >  0.5*periodicLength)
            transV[1] = -periodicLength;
        if(dy < -0.5*periodicLength)
            transV[1] =  periodicLength;

        if(dz >  0.5*periodicLength)
            transV[2] = -periodicLength;
        if(dz < -0.5*periodicLength)
            transV[2] =  periodicLength;
}

bool PrismVolume::validateCollisionFast(Facet* fit, Vertex* vit, Polyhedron *p, double periodicLength)
{
    if(periodicLength>0)
    {
        //Polyhedron* p = MeshManager::getInstance()->getMesh()->getSlaveMesh();
	Vertex* vit2 = p->v_array[fit->vi[0]];
        
        Vector3d transV(0,0,0);
        getTranslationVector(transV, vit, vit2, periodicLength);
                
        for (size_t i=0; i<3; ++i)
	{
		if (vit->min[i]+transV[i] > fit->max[i] || vit->max[i]+transV[i] < fit->min[i])
			return false;
	}

	Vector3d min(fit->min[0], fit->min[1], fit->min[2]);
	Vector3d max(fit->max[0], fit->max[1], fit->max[2]);
	
        Vector3d vitpori = vit->Pori()+transV;
        Vector3d vitp = vit->P()+transV;
	return TestSegmentAABB(vitpori, vitp, min, max);
    }
    else
    {
	for (size_t i=0; i<3; ++i)
	{
		if (vit->min[i] > fit->max[i] || vit->max[i] < fit->min[i])
			return false;
	}

	Vector3d min(fit->min[0], fit->min[1], fit->min[2]);
	Vector3d max(fit->max[0], fit->max[1], fit->max[2]);
	
	return TestSegmentAABB(vit->Pori(), vit->P(), min, max);
    }
}

bool PrismVolume::validateCollisionEEFastFast(Halfedge* eit1, Halfedge* eit2)
{
	int v11 = eit1->vertex()->getIndex();
	int v12 = eit1->opposite()->vertex()->getIndex();
	int v21 = eit2->vertex()->getIndex();
	int v22 = eit2->opposite()->vertex()->getIndex();
	
	//cerr << "Checking: " << v11 << " " << v12 << "---" << v21 << " " << v22 << endl;
	
	if ((v11 == v21) || (v11 == v22) || (v12 == v21) || (v12 == v22))
	{
	//	cerr << "false" << endl;
        	std::cout<<"worng!!!"<<std::endl;
		return false;
	}

	//cerr << "true" << endl;
	return true;
}

bool PrismVolume::validateCollisionEEFast(Halfedge* eit1, Halfedge* eit2, Polyhedron *p, double periodicLength, double h)
{
    if(periodicLength > 0)
    {
        //Polyhedron* p = MeshManager::getInstance()->getMesh()->getSlaveMesh();
	int v1 = eit1->vertex()->getIndex();
	int v2 = eit2->vertex()->getIndex();
	Vertex *vh1 = p->v_array[v1];
	Vertex *vh2 = p->v_array[v2];
        
        Vector3d transV(0,0,0);
        getTranslationVector(transV, vh1, vh2, periodicLength);

        for (size_t i=0; i<3; ++i)
	{
		if (eit1->min[i]+transV[i] > eit2->max[i] || eit1->max[i]+transV[i] < eit2->min[i])
			return false;
	}

	//Vector3d min(eit2->min[0], eit2->min[1], eit2->min[2]);
	//Vector3d max(eit2->max[0], eit2->max[1], eit2->max[2]);

	//return (TestSegmentAABB(eit1->vertex()->Pori(), eit1->vertex()->P(), min, max) || 
        //    TestSegmentAABB(eit1->opposite()->vertex()->Pori(), eit1->opposite()->vertex()->P(), min, max));

        // TODO: Add fast detection code
        /*
        Vector3d e1v1s = eit1->vertex()->Pori()+transV;             Vector3d e1v1e = eit1->vertex()->P()+transV;
        Vector3d e1v2s = eit1->opposite()->vertex()->Pori()+transV; Vector3d e1v2e = eit1->opposite()->vertex()->P()+transV;
        Vector3d e2v1s = eit2->vertex()->Pori();                    Vector3d e2v1e = eit2->vertex()->P();
        Vector3d e2v2s = eit2->opposite()->vertex()->Pori();        Vector3d e2v2e = eit2->opposite()->vertex()->P();
        return TestSegmentSubdivideSpace(e1v1s, e1v1e, e1v2s, e1v2e, e2v1s, e2v1e, e2v2s, e2v2e, h);
        */
        return true;
    }
    else
    {
	for (size_t i=0; i<3; ++i)
	{
		if (eit1->min[i] > eit2->max[i] || eit1->max[i] < eit2->min[i])
			return false;
	}

	//Vector3d min(eit2->min[0], eit2->min[1], eit2->min[2]);
	//Vector3d max(eit2->max[0], eit2->max[1], eit2->max[2]);

	//return (TestSegmentAABB(eit1->vertex()->Pori(), eit1->vertex()->P(), min, max) || 
        //    TestSegmentAABB(eit1->opposite()->vertex()->Pori(), eit1->opposite()->vertex()->P(), min, max));

        // TODO: Add fast detection code
        /*
        Vector3d e1v1s = eit1->vertex()->Pori();             Vector3d e1v1e = eit1->vertex()->P();
        Vector3d e1v2s = eit1->opposite()->vertex()->Pori(); Vector3d e1v2e = eit1->opposite()->vertex()->P();
        Vector3d e2v1s = eit2->vertex()->Pori();             Vector3d e2v1e = eit2->vertex()->P();
        Vector3d e2v2s = eit2->opposite()->vertex()->Pori(); Vector3d e2v2e = eit2->opposite()->vertex()->P();
        return TestSegmentSubdivideSpace(e1v1s, e1v1e, e1v2s, e1v2e, e2v1s, e2v1e, e2v2s, e2v2e, h);
        */
        return true;
    }
}

bool PrismVolume::validateCollisionEEAndReplace(Halfedge* eit1, Halfedge* eit2, double h, Collision& c, Polyhedron *p, double periodicLength)
{
    if(periodicLength > 0)
    {
	//Polyhedron* p = MeshManager::getInstance()->getMesh()->getSlaveMesh();
	
	int v11 = eit1->vertex()->getIndex();
	int v12 = eit1->opposite()->vertex()->getIndex();
	int v21 = eit2->vertex()->getIndex();
	int v22 = eit2->opposite()->vertex()->getIndex();
	
	Vertex *vh11 = p->v_array[v11];
	Vertex *vh12 = p->v_array[v12];
	Vertex *vh21 = p->v_array[v21];
	Vertex *vh22 = p->v_array[v22];
        
        Vector3d transV(0,0,0);
        getTranslationVector(transV, vh11, vh21, periodicLength);
	
	double tRet;
        Vector3d v11pori = vh11->Pori()+transV;
        Vector3d v12pori = vh12->Pori()+transV;
        Vector3d v11p = vh11->P()+transV;
        Vector3d v12p = vh12->P()+transV;
	if (validateCollision(v11pori, v12pori, vh21->Pori(), vh22->Pori(),
			      v11p,    v12p,    vh21->P(),    vh22->P(), h, tRet))
	{
		c = Collision(v11, v12, v21, v22, tRet, false);
		return true;
	}
	
	return false;
    }
    else
    {
        //Polyhedron* p = MeshManager::getInstance()->getMesh()->getSlaveMesh();
	
	int v11 = eit1->vertex()->getIndex();
	int v12 = eit1->opposite()->vertex()->getIndex();
	int v21 = eit2->vertex()->getIndex();
	int v22 = eit2->opposite()->vertex()->getIndex();
	
	Vertex *vh11 = p->v_array[v11];
	Vertex *vh12 = p->v_array[v12];
	Vertex *vh21 = p->v_array[v21];
	Vertex *vh22 = p->v_array[v22];
	
	double tRet;
	if (validateCollision(vh11->Pori(), vh12->Pori(), vh21->Pori(), vh22->Pori(),
						  vh11->P(), vh12->P(), vh21->P(), vh22->P(), h, tRet))
	{
		c = Collision(v11, v12, v21, v22, tRet, false);
		return true;
	}
	
	return false;
    }
}

void PrismVolume::computeVolumeWithHash(Polyhedron& m, AABB& aabb, double thickness, double edgeLength,
vector<set<unsigned int> > &subMeshes, vector<Collision> &collisions, vector<Collision> &collisionsEE, Intersection::Hash &h, int nv_resident, double periodicLength)
{
	h.resetHash(aabb.getMin(), aabb.getMax(), edgeLength, thickness, nv_resident, periodicLength);

        ////cout<<"Begin fill hash grid "<<endl;
        double t_fill_start = omp_get_wtime();
	h.fill(&m, subMeshes);
        double t_fill_end = omp_get_wtime();
        cout<<"fill time: "<<t_fill_end-t_fill_start<<"\n";
        ////cout<<"End fill hash grid "<<endl;

        ////cout<<"Begin sort hash grid "<<endl;
        double t_sort_start = omp_get_wtime();
	h.sort();
        double t_sort_end = omp_get_wtime();
        cout<<"sort time: "<<t_sort_end-t_sort_start<<"\n";
        ////cout<<"End sort hash grid "<<endl;

        ////cout<<"Begin check hash grid collision "<<endl;
        double t_compute_start = omp_get_wtime();
	h.computeVolumeParallel(collisions, collisionsEE);
        double t_compute_end = omp_get_wtime();
        cout<<"compute volume time: "<<t_compute_end-t_compute_start<<"\n";
        ////cout<<"End check hash grid collision "<<endl;
        h.clear();
}

bool PrismVolume::validateCollision(Vector3d &s11, Vector3d &s12, Vector3d &s21, Vector3d &s22,
								    Vector3d &e11, Vector3d &e12, Vector3d &e21, Vector3d &e22, double h, double& tret)
{
	// If thickness is non-zero, we are doing proximity detection, otherwise do normal collision detection
	//
	if (h)
	{
		//assert(0);
                //cout<<"thickness not zero!"<<endl;
                Vector3d v11 = e11 - s11;
                Vector3d v12 = e12 - s12;
                Vector3d v21 = e21 - s21;
                Vector3d v22 = e22 - s22;

                bool found = EdgeEdgeProximity(s11, s12, s21, s22, v11, v12, v21, v22, h, tret);

                return found;

	}

	Vector3d v11 = e11 - s11;
	Vector3d v12 = e12 - s12;
	Vector3d v21 = e21 - s21;
	Vector3d v22 = e22 - s22;

	bool found = EdgeEdgeIntersection(s11, s12, s21, s22, v11, v12, v21, v22, tret);

	return found;
}

bool PrismVolume::validateCollision(Vertex*& vh, Vector3d& before0, Vector3d& before1, Vector3d& before2,
									Vector3d& after0, Vector3d& after1, Vector3d& after2,
									double h, double& tret)
{
	// If thickness is non-zero, we are doing proximity detection, otherwise do normal collision detection
	//
	if (h)
	{
		Vector3d v0 = after0 - before0;
		Vector3d v1 = after1 - before1;
		Vector3d v2 = after2 - before2;
		Vector3d v3 = vh->P() - vh->Pori();

		bool found = VertexFaceProximity(before0, before1, before2, vh->Pori(),
										 v0, v1, v2, v3, h, tret);

		return found;
	}

	Vector3d v0 = after0 - before0;
	Vector3d v1 = after1 - before1;
	Vector3d v2 = after2 - before2;
	Vector3d v3 = vh->P() - vh->Pori();

	bool found = VertexFaceIntersection(before0, before1, before2, vh->Pori(),
							            v0, v1, v2, v3, tret);
	
	return found;
}

bool PrismVolume::VertexFaceIntersection(Vector3d &x0,Vector3d &x1,Vector3d &x2,Vector3d &x3,
									     Vector3d &v0,Vector3d &v1,Vector3d &v2,Vector3d &v3, double &tRet)
{
	double t[4];
	int count = areCoplanar(x0, x1, x2, x3, v0, v1, v2, v3, t);
	
	for (int k=0; k<count; ++k)
	{
		if (t[k] > 0.0 && t[k] <= 1.0)
		{
			Vector3d p0 = x0 + t[k] * v0;
			Vector3d p1 = x1 + t[k] * v1;
			Vector3d p2 = x2 + t[k] * v2;
			Vector3d p3 = x3 + t[k] * v3;
			
			double t1, t2, t3;
			double distance = vertexFaceDist(p3, p0, p1, p2, t1, t2, t3);
			
			Vector3d relVel = v3 - (t1 * v0 + t2 * v1 + t3 * v2);
			Vector3d e1 = p1 - p0;
			Vector3d e2 = p2 - p0;
			Vector3d n  = e1.cross(e2);
					
			double Vn = relVel.dot(n);
			//double epsilon = Vn * Vn * 1e-10;
			double epsilon = 1e-9;
			
			if (distance < epsilon)
			{
				//tRet = max(t[k]-0.01, 0.0);
				//tRet = max(t[k]-0.02, 0.0);
				tRet = t[k];

				return true;
			}
		}
	}
	
	return false;
	
}

bool PrismVolume::EdgeEdgeIntersection(Vector3d &x0,Vector3d &x1,Vector3d &x2,Vector3d &x3,
									   Vector3d &v0,Vector3d &v1,Vector3d &v2,Vector3d &v3, double &tRet)
{
	double t[4];
	int count = areCoplanar(x0, x1, x2, x3, v0, v1, v2, v3, t);
	
	for (int k=0; k<count; ++k)
	{
		if (t[k] > 0.0 && t[k] <= 1.0)
		{
			Vector3d p0 = x0 + t[k] * v0;
			Vector3d p1 = x1 + t[k] * v1;
			Vector3d p2 = x2 + t[k] * v2;
			Vector3d p3 = x3 + t[k] * v3;
			
			double r, s;
			double distance = edgeEdgeDist(p0, p1, p2, p3, r, s);
			
			Vector3d relVel = (v0 * (1.0 - r) + v1 * r) - (v2 * (1.0 - s) + v3 * s);
			Vector3d e1 = p1 - p0;
			Vector3d e2 = p3 - p2;
			Vector3d n  = e1.cross(e2);
					
			double Vn = relVel.dot(n);
			//double epsilon = Vn * Vn * 1e-10;
			double epsilon = 1e-9;
			
			if (distance < epsilon)
			{
				//tRet = max(t[k]-0.01, 0.0);
				//tRet = max(t[k]-0.02, 0.0);
				tRet = t[k];

				return true;
			}
		}
	}
	
	return false;
	
}

void PrismVolume::getBarycentricCoordinates(Vector3d &x0, Vector3d &x1, Vector3d &x2, Vector3d &x3, bool vt, double w[4])
{
      if (vt)
      {
        /*
	Vector3d v0 = x1 - x0;
	Vector3d v1 = x2 - x0;
	Vector3d v2 = x3 - x0;

	double dot00 = v0.dot(v0);
	double dot01 = v0.dot(v1);
	double dot02 = v0.dot(v2);
	double dot11 = v1.dot(v1);
	double dot12 = v1.dot(v2);

	double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);

	w[1] = (dot11 * dot02 - dot01 * dot12) * invDenom;
	w[2] = (dot00 * dot12 - dot01 * dot02) * invDenom;
	
        w[0] = 1.0 - w[1] - w[2];
	w[3] = 1.0;
        
        if (w[0] < 0.0)
          w[0] = 0.0;

        if (w[0] > 1.0)
          w[0] = 1.0;

        if (w[1] < 0.0)
          w[1] = 0.0;
        
        if (w[1] > 1.0)
          w[1] = 1.0;

        if (w[2] < 0.0)
          w[2] = 0.0;

        if (w[2] > 1.0)
          w[2] = 1.0;
        */
        double ab[3], ac[3], ap[3], bp[3];
	
	ab[0] = x1[0] - x0[0];
	ab[1] = x1[1] - x0[1];
	ab[2] = x1[2] - x0[2];
	
	ac[0] = x2[0] - x0[0];
	ac[1] = x2[1] - x0[1];
	ac[2] = x2[2] - x0[2];
	
	ap[0] = x3[0] - x0[0];
	ap[1] = x3[1] - x0[1];
	ap[2] = x3[2] - x0[2];
	
	double d1 = ab[0]*ap[0] + ab[1]*ap[1] + ab[2]*ap[2];
	double d2 = ac[0]*ap[0] + ac[1]*ap[1] + ac[2]*ap[2];
	
	if ((d1 <= 0.0f) && (d2 <= 0.0f))
	{
		w[0] = 1.0f;
		w[1] = 0.0f;
		w[2] = 0.0f;
                w[3] = 1.0f;
		
                return;
	}
	
	bp[0] = x3[0] - x1[0];
	bp[1] = x3[1] - x1[1];
	bp[2] = x3[2] - x1[2];
	
	double d3 = ab[0]*bp[0] + ab[1]*bp[1] + ab[2]*bp[2];
	double d4 = ac[0]*bp[0] + ac[1]*bp[1] + ac[2]*bp[2];
	
	if ((d3 >= 0.0f) && (d4 <= d3))
	{
		w[0] = 0.0f;
		w[1] = 1.0f;
		w[2] = 0.0f;
                w[3] = 1.0f;
		
                return;
	}
	
	double vc = d1*d4 - d3*d2;
	
	if ((vc <= 0.0f) && (d1 >= 0.0f) && (d3 <= 0.0f))
	{
		double v = d1 / (d1 - d3);
		
		w[0] = 1-v;
		w[1] = v;
		w[2] = 0;
                w[3] = 1.0;
		
                return;
	}
	
	double cp[3];
	cp[0] = x3[0] - x2[0];
	cp[1] = x3[1] - x2[1];
	cp[2] = x3[2] - x2[2];
	
	double d5 = ab[0]*cp[0] + ab[1]*cp[1] + ab[2]*cp[2];
	double d6 = ac[0]*cp[0] + ac[1]*cp[1] + ac[2]*cp[2];
	
	if ((d6 >= 0.0f) && (d5 <= d6))
	{
		w[0] = 0;
		w[1] = 0;
		w[2] = 1;
                w[3] = 1;
		
                return;
	}
	
	double vb = d5*d2 - d1*d6;
	
	if ((vb <= 0.0f) && (d2 >= 0.0f) && (d6 <= 0.0f))
	{
		double v = d2 / (d2 - d6);
		
		w[0] = 1-v;
		w[1] = 0;
		w[2] = v;
		w[3] = 1;

                return;
	}
	
	double va = d3*d6 - d5*d4;
	
	if ((va <= 0.0f) && ((d4-d3) >= 0.0f) && ((d5-d6) >= 0.0f))
	{
		double v = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		
		w[0] = 0;
		w[1] = 1-v;
		w[2] = v;
		w[3] = 1.0;

                return;
	}
	
	double denom = 1.0f / (va + vb + vc);
	w[1] = vb * denom;
	w[2] = vc * denom;
	w[0] = 1.0 - w[1] - w[2];
        w[3] = 1.0;
      }
      else
      {
	Vector3d v0 = x1 - x0;
	Vector3d v1 = x3 - x2;
	Vector3d v2 = x0 - x2;

	double dot00 = v0.dot(v0);
	double dot01 = v0.dot(v1);
	double dot02 = v0.dot(v2);
	double dot11 = v1.dot(v1);
	double dot12 = v1.dot(v2);

	//double invDenom = 1.0 / (dot00 * dot11 + dot01 * dot01);
	double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);

        w[1] = (-dot11 * dot02 + dot01 * dot12) * invDenom;
	w[3] = (dot00 * dot12 - dot01 * dot02) * invDenom;
       
        if (w[1] < 0.0)
          w[1] = 0.0;
        if (w[1] > 1.0)
          w[1] = 1.0;
        if (w[3] < 0.0)
          w[3] = 0.0;
        if (w[3] > 1.0)
          w[3] = 1.0;

	w[0] = 1.0 - w[1];
	w[2] = 1.0 - w[3];
	//w[1] = (dot11 * dot02 + dot01 * dot12) * invDenom;
	//w[3] = (dot00 * dot12 + dot01 * dot02) * invDenom;
	//w[0] = 1.0 - w[1];
	//w[2] = 1.0 - w[3];
      }
}

void PrismVolume::getEETimeGradient(Vector3d &x1, Vector3d &x2, Vector3d &x3, Vector3d &x4,
								  Vector3d &v1, Vector3d &v2, Vector3d &v3, Vector3d &v4, double w[4], double t, double h,
								  Vector3d tGrad[4])
{
	Vector3d x10 = x2 - x1;
	Vector3d x20 = x3 - x1;
	Vector3d x30 = x4 - x1;
	Vector3d v10 = v2 - v1;
	Vector3d v20 = v3 - v1;
	Vector3d v30 = v4 - v1;
	
	Vector3d p = x10.cross(x20);
	Vector3d q = v10.cross(x20) + x10.cross(v20);
	Vector3d r = v10.cross(v20);

	double a = v30.dot(r);
	double b = x30.dot(r) + v30.dot(q);
	double c = x30.dot(q) + v30.dot(p);
	double d = x30.dot(p);

	double scale = 1.0 / (3.0 * a * t * t + 2.0 * b * t + c);

	Vector3d z[4];
	z[0] = x1 + v1 * t;
	z[1] = x2 + v2 * t;
	z[2] = x3 + v3 * t;
	z[3] = x4 + v4 * t;

	Vector3d n = t * scale * (z[1] - z[0]).cross(z[3] - z[2]);
	
	tGrad[0] = -w[0] * n;
	tGrad[1] = -w[1] * n;
        tGrad[2] = w[2] * n;
	tGrad[3] = w[3] * n;
}

void PrismVolume::getTimeGradient(Vector3d &x1, Vector3d &x2, Vector3d &x3, Vector3d &x4,
								  Vector3d &v1, Vector3d &v2, Vector3d &v3, Vector3d &v4, double w[4], double t, double h,
								  Vector3d tGrad[4])
{
	Vector3d x10 = x2 - x1;
	Vector3d x20 = x3 - x1;
	Vector3d x30 = x4 - x1;
	Vector3d v10 = v2 - v1;
	Vector3d v20 = v3 - v1;
	Vector3d v30 = v4 - v1;
	
	Vector3d p = x10.cross(x20);
	Vector3d q = v10.cross(x20) + x10.cross(v20);
	Vector3d r = v10.cross(v20);

	double a = v30.dot(r);
	double b = x30.dot(r) + v30.dot(q);
	double c = x30.dot(q) + v30.dot(p);
	double d = x30.dot(p);

	double scale = 1.0 / (3.0 * a * t * t + 2.0 * b * t + c);

	Vector3d z[4];
	z[0] = x1 + v1 * t;
	z[1] = x2 + v2 * t;
	z[2] = x3 + v3 * t;

	Vector3d n = t * scale * (z[1] - z[0]).cross(z[2] - z[0]);
	
	tGrad[0] = w[0] * n;
	tGrad[1] = w[1] * n;
	tGrad[2] = w[2] * n;
	tGrad[3] = -w[3] * n;
}
	
void PrismVolume::getEEBarycentricGradients(Vector3d &x0, Vector3d &x1, Vector3d &x2, Vector3d &x3, Vector3d wGrad[4][4])
{

  double w0dx[12];
  double w1dx[12];
  double w2dx[12];
  double w3dx[12];


  w0dx[0] = ((x0[0]*-2.0+x1[0]+x2[0])*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))+(x2[0]-x3[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x2[0]-x3[0])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x2[0]-x3[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[0]*2.0-x1[0]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w0dx[1] = ((x0[1]*-2.0+x1[1]+x2[1])*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))+(x2[1]-x3[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x2[1]-x3[1])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x2[1]-x3[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[1]*2.0-x1[1]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w0dx[2] = ((x0[2]*-2.0+x1[2]+x2[2])*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))+(x2[2]-x3[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x2[2]-x3[2])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x2[2]-x3[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[2]*2.0-x1[2]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w0dx[3] = ((x0[0]-x2[0])*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-(x2[0]-x3[0])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x2[0]-x3[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[0]*2.0-x1[0]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w0dx[4] = ((x0[1]-x2[1])*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-(x2[1]-x3[1])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x2[1]-x3[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[1]*2.0-x1[1]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w0dx[5] = ((x0[2]-x2[2])*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-(x2[2]-x3[2])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x2[2]-x3[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[2]*2.0-x1[2]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w0dx[6] = (-(x2[0]*2.0-x3[0]*2.0)*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[0]-x2[0]*2.0+x3[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x0[0]-x1[0])*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))+(x0[0]-x1[0])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x0[0]-x1[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[0]*2.0-x3[0]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w0dx[7] = (-(x2[1]*2.0-x3[1]*2.0)*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[1]-x2[1]*2.0+x3[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x0[1]-x1[1])*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))+(x0[1]-x1[1])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x0[1]-x1[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[1]*2.0-x3[1]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w0dx[8] = (-(x2[2]*2.0-x3[2]*2.0)*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[2]-x2[2]*2.0+x3[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x0[2]-x1[2])*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))+(x0[2]-x1[2])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x0[2]-x1[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[2]*2.0-x3[2]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w0dx[9] = -(-(x2[0]*2.0-x3[0]*2.0)*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[0]-x2[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x0[0]-x1[0])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x0[0]-x1[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[0]*2.0-x3[0]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w0dx[10] = -(-(x2[1]*2.0-x3[1]*2.0)*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[1]-x2[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x0[1]-x1[1])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x0[1]-x1[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[1]*2.0-x3[1]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w0dx[11] = -(-(x2[2]*2.0-x3[2]*2.0)*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[2]-x2[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x0[2]-x1[2])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x0[2]-x1[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[2]*2.0-x3[2]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w1dx[0] = -((x0[0]*-2.0+x1[0]+x2[0])*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))+(x2[0]-x3[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x2[0]-x3[0])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x2[0]-x3[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[0]*2.0-x1[0]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w1dx[1] = -((x0[1]*-2.0+x1[1]+x2[1])*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))+(x2[1]-x3[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x2[1]-x3[1])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x2[1]-x3[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[1]*2.0-x1[1]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w1dx[2] = -((x0[2]*-2.0+x1[2]+x2[2])*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))+(x2[2]-x3[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x2[2]-x3[2])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x2[2]-x3[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[2]*2.0-x1[2]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w1dx[3] = -((x0[0]-x2[0])*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-(x2[0]-x3[0])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x2[0]-x3[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[0]*2.0-x1[0]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w1dx[4] = -((x0[1]-x2[1])*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-(x2[1]-x3[1])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x2[1]-x3[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[1]*2.0-x1[1]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w1dx[5] = -((x0[2]-x2[2])*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-(x2[2]-x3[2])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x2[2]-x3[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[2]*2.0-x1[2]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w1dx[6] = -(-(x2[0]*2.0-x3[0]*2.0)*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[0]-x2[0]*2.0+x3[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x0[0]-x1[0])*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))+(x0[0]-x1[0])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x0[0]-x1[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[0]*2.0-x3[0]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w1dx[7] = -(-(x2[1]*2.0-x3[1]*2.0)*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[1]-x2[1]*2.0+x3[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x0[1]-x1[1])*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))+(x0[1]-x1[1])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x0[1]-x1[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[1]*2.0-x3[1]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w1dx[8] = -(-(x2[2]*2.0-x3[2]*2.0)*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[2]-x2[2]*2.0+x3[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x0[2]-x1[2])*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))+(x0[2]-x1[2])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x0[2]-x1[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[2]*2.0-x3[2]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w1dx[9] = (-(x2[0]*2.0-x3[0]*2.0)*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[0]-x2[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x0[0]-x1[0])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x0[0]-x1[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[0]*2.0-x3[0]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w1dx[10] = (-(x2[1]*2.0-x3[1]*2.0)*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[1]-x2[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x0[1]-x1[1])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x0[1]-x1[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[1]*2.0-x3[1]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w1dx[11] = (-(x2[2]*2.0-x3[2]*2.0)*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[2]-x2[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x0[2]-x1[2])*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))-((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2])))*((x0[2]-x1[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[2]*2.0-x3[2]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w2dx[0] = ((x0[0]*2.0-x1[0]*2.0)*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))+(x0[0]*-2.0+x1[0]+x2[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x2[0]-x3[0])*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))-(x2[0]-x3[0])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x2[0]-x3[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[0]*2.0-x1[0]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w2dx[1] = ((x0[1]*2.0-x1[1]*2.0)*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))+(x0[1]*-2.0+x1[1]+x2[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x2[1]-x3[1])*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))-(x2[1]-x3[1])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x2[1]-x3[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[1]*2.0-x1[1]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w2dx[2] = ((x0[2]*2.0-x1[2]*2.0)*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))+(x0[2]*-2.0+x1[2]+x2[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x2[2]-x3[2])*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))-(x2[2]-x3[2])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x2[2]-x3[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[2]*2.0-x1[2]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w2dx[3] = (-(x0[0]*2.0-x1[0]*2.0)*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))+(x2[0]-x3[0])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[0]-x2[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x2[0]-x3[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[0]*2.0-x1[0]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w2dx[4] = (-(x0[1]*2.0-x1[1]*2.0)*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))+(x2[1]-x3[1])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[1]-x2[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x2[1]-x3[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[1]*2.0-x1[1]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w2dx[5] = (-(x0[2]*2.0-x1[2]*2.0)*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))+(x2[2]-x3[2])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[2]-x2[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x2[2]-x3[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[2]*2.0-x1[2]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w2dx[6] = ((x0[0]-x2[0]*2.0+x3[0])*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))-(x0[0]-x1[0])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[0]-x1[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x0[0]-x1[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[0]*2.0-x3[0]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w2dx[7] = ((x0[1]-x2[1]*2.0+x3[1])*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))-(x0[1]-x1[1])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[1]-x1[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x0[1]-x1[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[1]*2.0-x3[1]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w2dx[8] = ((x0[2]-x2[2]*2.0+x3[2])*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))-(x0[2]-x1[2])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[2]-x1[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x0[2]-x1[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[2]*2.0-x3[2]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w2dx[9] = -((x0[0]-x2[0])*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))-(x0[0]-x1[0])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x0[0]-x1[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[0]*2.0-x3[0]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w2dx[10] = -((x0[1]-x2[1])*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))-(x0[1]-x1[1])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x0[1]-x1[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[1]*2.0-x3[1]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w2dx[11] = -((x0[2]-x2[2])*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))-(x0[2]-x1[2])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x0[2]-x1[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[2]*2.0-x3[2]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w3dx[0] = -((x0[0]*2.0-x1[0]*2.0)*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))+(x0[0]*-2.0+x1[0]+x2[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x2[0]-x3[0])*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))-(x2[0]-x3[0])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x2[0]-x3[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[0]*2.0-x1[0]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w3dx[1] = -((x0[1]*2.0-x1[1]*2.0)*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))+(x0[1]*-2.0+x1[1]+x2[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x2[1]-x3[1])*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))-(x2[1]-x3[1])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x2[1]-x3[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[1]*2.0-x1[1]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w3dx[2] = -((x0[2]*2.0-x1[2]*2.0)*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))+(x0[2]*-2.0+x1[2]+x2[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))+(x2[2]-x3[2])*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))-(x2[2]-x3[2])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x2[2]-x3[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[2]*2.0-x1[2]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w3dx[3] = -(-(x0[0]*2.0-x1[0]*2.0)*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))+(x2[0]-x3[0])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[0]-x2[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x2[0]-x3[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[0]*2.0-x1[0]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w3dx[4] = -(-(x0[1]*2.0-x1[1]*2.0)*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))+(x2[1]-x3[1])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[1]-x2[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x2[1]-x3[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[1]*2.0-x1[1]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w3dx[5] = -(-(x0[2]*2.0-x1[2]*2.0)*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))+(x2[2]-x3[2])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[2]-x2[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x2[2]-x3[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x0[2]*2.0-x1[2]*2.0)*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0)));
  w3dx[6] = -((x0[0]-x2[0]*2.0+x3[0])*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))-(x0[0]-x1[0])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[0]-x1[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x0[0]-x1[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[0]*2.0-x3[0]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w3dx[7] = -((x0[1]-x2[1]*2.0+x3[1])*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))-(x0[1]-x1[1])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[1]-x1[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x0[1]-x1[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[1]*2.0-x3[1]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w3dx[8] = -((x0[2]-x2[2]*2.0+x3[2])*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))-(x0[2]-x1[2])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))+(x0[2]-x1[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))-1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x0[2]-x1[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[2]*2.0-x3[2]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w3dx[9] = ((x0[0]-x2[0])*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))-(x0[0]-x1[0])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x0[0]-x1[0])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[0]*2.0-x3[0]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w3dx[10] = ((x0[1]-x2[1])*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))-(x0[1]-x1[1])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x0[1]-x1[1])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[1]*2.0-x3[1]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
  w3dx[11] = ((x0[2]-x2[2])*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))-(x0[2]-x1[2])*((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2])))/((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0))+1.0/pow((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*(pow(x2[0]-x3[0],2.0)+pow(x2[1]-x3[1],2.0)+pow(x2[2]-x3[2],2.0))-pow((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]),2.0),2.0)*((pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0))*((x0[0]-x2[0])*(x2[0]-x3[0])+(x0[1]-x2[1])*(x2[1]-x3[1])+(x0[2]-x2[2])*(x2[2]-x3[2]))-((x0[0]-x1[0])*(x0[0]-x2[0])+(x0[1]-x1[1])*(x0[1]-x2[1])+(x0[2]-x1[2])*(x0[2]-x2[2]))*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2])))*((x0[2]-x1[2])*((x0[0]-x1[0])*(x2[0]-x3[0])+(x0[1]-x1[1])*(x2[1]-x3[1])+(x0[2]-x1[2])*(x2[2]-x3[2]))*2.0-(x2[2]*2.0-x3[2]*2.0)*(pow(x0[0]-x1[0],2.0)+pow(x0[1]-x1[1],2.0)+pow(x0[2]-x1[2],2.0)));
 
  wGrad[0][0] = Vector3d(w0dx[0], w0dx[1], w0dx[2]);
  wGrad[0][1] = Vector3d(w0dx[3], w0dx[4], w0dx[5]);
  wGrad[0][2] = Vector3d(w0dx[6], w0dx[7], w0dx[8]);
  wGrad[0][3] = Vector3d(w0dx[9], w0dx[10], w0dx[11]);

  wGrad[1][0] = Vector3d(w1dx[0], w1dx[1], w1dx[2]);
  wGrad[1][1] = Vector3d(w1dx[3], w1dx[4], w1dx[5]);
  wGrad[1][2] = Vector3d(w1dx[6], w1dx[7], w1dx[8]);
  wGrad[1][3] = Vector3d(w1dx[9], w1dx[10], w1dx[11]);

  wGrad[2][0] = Vector3d(w2dx[0], w2dx[1], w2dx[2]);
  wGrad[2][1] = Vector3d(w2dx[3], w2dx[4], w2dx[5]);
  wGrad[2][2] = Vector3d(w2dx[6], w2dx[7], w2dx[8]);
  wGrad[2][3] = Vector3d(w2dx[9], w2dx[10], w2dx[11]);

  wGrad[3][0] = Vector3d(w3dx[0], w3dx[1], w3dx[2]);
  wGrad[3][1] = Vector3d(w3dx[3], w3dx[4], w3dx[5]);
  wGrad[3][2] = Vector3d(w3dx[6], w3dx[7], w3dx[8]);
  wGrad[3][3] = Vector3d(w3dx[9], w3dx[10], w3dx[11]);

}

void PrismVolume::getBarycentricGradients(Vector3d &x0, Vector3d &x1, Vector3d &x2, Vector3d &x3, Vector3d wGrad[4][4])
{
	// Gradients of the barycentric coordinates
	//

double u[12] =
{(-2*x0[0] + x1[0] + x3[0])/
     (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))\
     - (2*(-x0[0] + x1[0])*(-((-x0[0] + x1[0])*(-x0[0] + x3[0])) - 
         (-x0[1] + x1[1])*(-x0[1] + x3[1]) - (-x0[2] + x1[2])*(-x0[2] + x3[2])))/
     pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
       pow(-x0[2] + x1[2],2),2) + 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(2*x0[0] - x1[0] - x3[0]) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(2*x0[0] - x2[0] - x3[0]) + 
         (2*x0[0] - x1[0] - x2[0])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         2*(-x0[0] + x1[0])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) - 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       (2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[0] + x2[0]) - 
         2*(-2*x0[0] + x1[0] + x2[0])*
          ((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2])) + 
         2*(-x0[0] + x1[0])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)),2)) + 
    ((-2*x0[0] + x1[0] + x2[0])*((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) - 
    (2*(-x0[0] + x1[0])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
         (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     (pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
         pow(-x0[2] + x1[2],2),2)*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   (-2*x0[1] + x1[1] + x3[1])/
     (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))\
     - (2*(-x0[1] + x1[1])*(-((-x0[0] + x1[0])*(-x0[0] + x3[0])) - 
         (-x0[1] + x1[1])*(-x0[1] + x3[1]) - (-x0[2] + x1[2])*(-x0[2] + x3[2])))/
     pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
       pow(-x0[2] + x1[2],2),2) + 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(2*x0[1] - x1[1] - x3[1]) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(2*x0[1] - x2[1] - x3[1]) + 
         (2*x0[1] - x1[1] - x2[1])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         2*(-x0[1] + x1[1])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) - 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       (2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[1] + x2[1]) - 
         2*(-2*x0[1] + x1[1] + x2[1])*
          ((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2])) + 
         2*(-x0[1] + x1[1])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)),2)) + 
    ((-2*x0[1] + x1[1] + x2[1])*((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) - 
    (2*(-x0[1] + x1[1])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
         (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     (pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
         pow(-x0[2] + x1[2],2),2)*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   (-2*x0[2] + x1[2] + x3[2])/
     (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))\
     - (2*(-x0[2] + x1[2])*(-((-x0[0] + x1[0])*(-x0[0] + x3[0])) - 
         (-x0[1] + x1[1])*(-x0[1] + x3[1]) - (-x0[2] + x1[2])*(-x0[2] + x3[2])))/
     pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
       pow(-x0[2] + x1[2],2),2) + 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(2*x0[2] - x1[2] - x3[2]) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(2*x0[2] - x2[2] - x3[2]) + 
         (2*x0[2] - x1[2] - x2[2])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         2*(-x0[2] + x1[2])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) - 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       (2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[2] + x2[2]) - 
         2*(-2*x0[2] + x1[2] + x2[2])*
          ((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2])) + 
         2*(-x0[2] + x1[2])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)),2)) + 
    ((-2*x0[2] + x1[2] + x2[2])*((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) - 
    (2*(-x0[2] + x1[2])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
         (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     (pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
         pow(-x0[2] + x1[2],2),2)*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   -((-x0[0] + x3[0])/
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))
       ) + (2*(-x0[0] + x1[0])*(-((-x0[0] + x1[0])*(-x0[0] + x3[0])) - 
         (-x0[1] + x1[1])*(-x0[1] + x3[1]) - (-x0[2] + x1[2])*(-x0[2] + x3[2])))/
     pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
       pow(-x0[2] + x1[2],2),2) + 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(-x0[0] + x3[0]) + 
         (-x0[0] + x2[0])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) - 
         2*(-x0[0] + x1[0])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) - 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       (-2*(x0[0] - x2[0])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])) - 
         2*(-x0[0] + x1[0])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)),2)) + 
    ((x0[0] - x2[0])*((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) + 
    (2*(-x0[0] + x1[0])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
         (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     (pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
         pow(-x0[2] + x1[2],2),2)*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   -((-x0[1] + x3[1])/
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))
       ) + (2*(-x0[1] + x1[1])*(-((-x0[0] + x1[0])*(-x0[0] + x3[0])) - 
         (-x0[1] + x1[1])*(-x0[1] + x3[1]) - (-x0[2] + x1[2])*(-x0[2] + x3[2])))/
     pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
       pow(-x0[2] + x1[2],2),2) + 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(-x0[1] + x3[1]) + 
         (-x0[1] + x2[1])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) - 
         2*(-x0[1] + x1[1])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) - 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       (2*(-x0[1] + x2[1])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])) - 
         2*(-x0[1] + x1[1])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)),2)) - 
    ((-x0[1] + x2[1])*((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) + 
    (2*(-x0[1] + x1[1])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
         (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     (pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
         pow(-x0[2] + x1[2],2),2)*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   -((-x0[2] + x3[2])/
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))
       ) + (2*(-x0[2] + x1[2])*(-((-x0[0] + x1[0])*(-x0[0] + x3[0])) - 
         (-x0[1] + x1[1])*(-x0[1] + x3[1]) - (-x0[2] + x1[2])*(-x0[2] + x3[2])))/
     pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
       pow(-x0[2] + x1[2],2),2) + 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(-x0[2] + x3[2]) + 
         (-x0[2] + x2[2])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) - 
         2*(-x0[2] + x1[2])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) - 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       (2*(-x0[2] + x2[2])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])) - 
         2*(-x0[2] + x1[2])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)),2)) - 
    ((-x0[2] + x2[2])*((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) + 
    (2*(-x0[2] + x1[2])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
         (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     (pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
         pow(-x0[2] + x1[2],2),2)*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[0] + x3[0]) - 
         (x0[0] - x1[0])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) - 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       (-2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[0] + x2[0]) - 
         2*(x0[0] - x1[0])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)),2)) + 
    ((x0[0] - x1[0])*((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[1] + x3[1]) + 
         (-x0[1] + x1[1])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) - 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       (-2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[1] + x2[1]) + 
         2*(-x0[1] + x1[1])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)),2)) - 
    ((-x0[1] + x1[1])*((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[2] + x3[2]) + 
         (-x0[2] + x1[2])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) - 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       (-2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[2] + x2[2]) + 
         2*(-x0[2] + x1[2])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)),2)) - 
    ((-x0[2] + x1[2])*((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   -((-x0[0] + x1[0])/
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))
       ) + (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[0] + x2[0]) + 
         (-x0[0] + x1[0])*(-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   -((-x0[1] + x1[1])/
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))
       ) + (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[1] + x2[1]) + 
         (-x0[1] + x1[1])*(-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   -((-x0[2] + x1[2])/
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))
       ) + (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[2] + x2[2]) + 
         (-x0[2] + x1[2])*(-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) };


double v[12] = 
{(-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(2*x0[0] - x1[0] - x3[0])) - 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*(2*x0[0] - x2[0] - x3[0]) - 
       (2*x0[0] - x1[0] - x2[0])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) - 
       2*(-x0[0] + x1[0])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2])))/
     (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2))) - 
    ((2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[0] + x2[0]) - 
         2*(-2*x0[0] + x1[0] + x2[0])*
          ((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2])) + 
         2*(-x0[0] + x1[0])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
            ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x3[2]))) - 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2)),2),
   (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(2*x0[1] - x1[1] - x3[1])) - 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*(2*x0[1] - x2[1] - x3[1]) - 
       (2*x0[1] - x1[1] - x2[1])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) - 
       2*(-x0[1] + x1[1])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2])))/
     (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2))) - 
    ((2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[1] + x2[1]) - 
         2*(-2*x0[1] + x1[1] + x2[1])*
          ((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2])) + 
         2*(-x0[1] + x1[1])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
            ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x3[2]))) - 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2)),2),
   (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(2*x0[2] - x1[2] - x3[2])) - 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*(2*x0[2] - x2[2] - x3[2]) - 
       (2*x0[2] - x1[2] - x2[2])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) - 
       2*(-x0[2] + x1[2])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2])))/
     (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2))) - 
    ((2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[2] + x2[2]) - 
         2*(-2*x0[2] + x1[2] + x2[2])*
          ((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2])) + 
         2*(-x0[2] + x1[2])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
            ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x3[2]))) - 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2)),2),
   (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(-x0[0] + x3[0])) - 
       (-x0[0] + x2[0])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
       2*(-x0[0] + x1[0])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2])))/
     (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2))) - 
    ((-2*(x0[0] - x2[0])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])) - 
         2*(-x0[0] + x1[0])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
            ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x3[2]))) - 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2)),2),
   (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(-x0[1] + x3[1])) - 
       (-x0[1] + x2[1])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
       2*(-x0[1] + x1[1])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2])))/
     (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2))) - 
    ((2*(-x0[1] + x2[1])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])) - 
         2*(-x0[1] + x1[1])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
            ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x3[2]))) - 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2)),2),
   (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(-x0[2] + x3[2])) - 
       (-x0[2] + x2[2])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
       2*(-x0[2] + x1[2])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2])))/
     (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2))) - 
    ((2*(-x0[2] + x2[2])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])) - 
         2*(-x0[2] + x1[2])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
            ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x3[2]))) - 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2)),2),
   (-((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
          (-x0[0] + x3[0])) + (x0[0] - x1[0])*
        ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
          (-x0[2] + x1[2])*(-x0[2] + x3[2])))/
     (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2))) - 
    ((-2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[0] + x2[0]) - 
         2*(x0[0] - x1[0])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])))*
       (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
            ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x3[2]))) - 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2)),2),
   (-((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
          (-x0[1] + x3[1])) - (-x0[1] + x1[1])*
        ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
          (-x0[2] + x1[2])*(-x0[2] + x3[2])))/
     (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2))) - 
    ((-2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[1] + x2[1]) + 
         2*(-x0[1] + x1[1])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])))*
       (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
            ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x3[2]))) - 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2)),2),
   (-((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
          (-x0[2] + x3[2])) - (-x0[2] + x1[2])*
        ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
          (-x0[2] + x1[2])*(-x0[2] + x3[2])))/
     (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2))) - 
    ((-2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[2] + x2[2]) + 
         2*(-x0[2] + x1[2])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])))*
       (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
            ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x3[2]))) - 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2)),2),
   (-((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
         (-x0[0] + x2[0])) - (-x0[0] + x1[0])*
       (-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
         (-x0[2] + x1[2])*(-x0[2] + x2[2])))/
    (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
      (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - pow(-x0[2] + x2[2],2))
      ),(-((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
           pow(-x0[2] + x1[2],2))*(-x0[1] + x2[1])) - 
      (-x0[1] + x1[1])*(-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
         (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2])))/
    (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
      (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - pow(-x0[2] + x2[2],2))
      ),(-((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
           pow(-x0[2] + x1[2],2))*(-x0[2] + x2[2])) - 
      (-x0[2] + x1[2])*(-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
         (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2])))/
    (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
      (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - pow(-x0[2] + x2[2],2))
      )};
 

//double w[12];
//for (size_t i=0; i<12; ++i)
//	w[i] = -u[i] - v[i];

double w[12] = 
{-((-2*x0[0] + x1[0] + x3[0])/
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))
       ) + (2*(-x0[0] + x1[0])*(-((-x0[0] + x1[0])*(-x0[0] + x3[0])) - 
         (-x0[1] + x1[1])*(-x0[1] + x3[1]) - (-x0[2] + x1[2])*(-x0[2] + x3[2])))/
     pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
       pow(-x0[2] + x1[2],2),2) - 
    (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(2*x0[0] - x1[0] - x3[0])) - 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*(2*x0[0] - x2[0] - x3[0]) - 
       (2*x0[0] - x1[0] - x2[0])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) - 
       2*(-x0[0] + x1[0])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2])))/
     (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2))) - 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(2*x0[0] - x1[0] - x3[0]) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(2*x0[0] - x2[0] - x3[0]) + 
         (2*x0[0] - x1[0] - x2[0])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         2*(-x0[0] + x1[0])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) + 
    ((2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[0] + x2[0]) - 
         2*(-2*x0[0] + x1[0] + x2[0])*
          ((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2])) + 
         2*(-x0[0] + x1[0])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
            ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x3[2]))) - 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2)),2) + 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       (2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[0] + x2[0]) - 
         2*(-2*x0[0] + x1[0] + x2[0])*
          ((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2])) + 
         2*(-x0[0] + x1[0])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)),2)) - 
    ((-2*x0[0] + x1[0] + x2[0])*((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) + 
    (2*(-x0[0] + x1[0])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
         (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     (pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
         pow(-x0[2] + x1[2],2),2)*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   -((-2*x0[1] + x1[1] + x3[1])/
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))
       ) + (2*(-x0[1] + x1[1])*(-((-x0[0] + x1[0])*(-x0[0] + x3[0])) - 
         (-x0[1] + x1[1])*(-x0[1] + x3[1]) - (-x0[2] + x1[2])*(-x0[2] + x3[2])))/
     pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
       pow(-x0[2] + x1[2],2),2) - 
    (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(2*x0[1] - x1[1] - x3[1])) - 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*(2*x0[1] - x2[1] - x3[1]) - 
       (2*x0[1] - x1[1] - x2[1])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) - 
       2*(-x0[1] + x1[1])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2])))/
     (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2))) - 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(2*x0[1] - x1[1] - x3[1]) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(2*x0[1] - x2[1] - x3[1]) + 
         (2*x0[1] - x1[1] - x2[1])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         2*(-x0[1] + x1[1])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) + 
    ((2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[1] + x2[1]) - 
         2*(-2*x0[1] + x1[1] + x2[1])*
          ((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2])) + 
         2*(-x0[1] + x1[1])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
            ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x3[2]))) - 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2)),2) + 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       (2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[1] + x2[1]) - 
         2*(-2*x0[1] + x1[1] + x2[1])*
          ((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2])) + 
         2*(-x0[1] + x1[1])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)),2)) - 
    ((-2*x0[1] + x1[1] + x2[1])*((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) + 
    (2*(-x0[1] + x1[1])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
         (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     (pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
         pow(-x0[2] + x1[2],2),2)*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   -((-2*x0[2] + x1[2] + x3[2])/
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))
       ) + (2*(-x0[2] + x1[2])*(-((-x0[0] + x1[0])*(-x0[0] + x3[0])) - 
         (-x0[1] + x1[1])*(-x0[1] + x3[1]) - (-x0[2] + x1[2])*(-x0[2] + x3[2])))/
     pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
       pow(-x0[2] + x1[2],2),2) - 
    (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(2*x0[2] - x1[2] - x3[2])) - 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*(2*x0[2] - x2[2] - x3[2]) - 
       (2*x0[2] - x1[2] - x2[2])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) - 
       2*(-x0[2] + x1[2])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2])))/
     (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2))) - 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(2*x0[2] - x1[2] - x3[2]) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(2*x0[2] - x2[2] - x3[2]) + 
         (2*x0[2] - x1[2] - x2[2])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         2*(-x0[2] + x1[2])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) + 
    ((2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[2] + x2[2]) - 
         2*(-2*x0[2] + x1[2] + x2[2])*
          ((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2])) + 
         2*(-x0[2] + x1[2])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
            ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x3[2]))) - 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2)),2) + 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       (2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[2] + x2[2]) - 
         2*(-2*x0[2] + x1[2] + x2[2])*
          ((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2])) + 
         2*(-x0[2] + x1[2])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)),2)) - 
    ((-2*x0[2] + x1[2] + x2[2])*((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) + 
    (2*(-x0[2] + x1[2])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
         (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     (pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
         pow(-x0[2] + x1[2],2),2)*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   (-x0[0] + x3[0])/
     (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))\
     - (2*(-x0[0] + x1[0])*(-((-x0[0] + x1[0])*(-x0[0] + x3[0])) - 
         (-x0[1] + x1[1])*(-x0[1] + x3[1]) - (-x0[2] + x1[2])*(-x0[2] + x3[2])))/
     pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
       pow(-x0[2] + x1[2],2),2) - 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(-x0[0] + x3[0]) + 
         (-x0[0] + x2[0])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) - 
         2*(-x0[0] + x1[0])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) - 
    (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(-x0[0] + x3[0])) - 
       (-x0[0] + x2[0])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
       2*(-x0[0] + x1[0])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2])))/
     (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2))) + 
    ((-2*(x0[0] - x2[0])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])) - 
         2*(-x0[0] + x1[0])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
            ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x3[2]))) - 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2)),2) + 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       (-2*(x0[0] - x2[0])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])) - 
         2*(-x0[0] + x1[0])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)),2)) - 
    ((x0[0] - x2[0])*((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) - 
    (2*(-x0[0] + x1[0])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
         (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     (pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
         pow(-x0[2] + x1[2],2),2)*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   (-x0[1] + x3[1])/
     (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))\
     - (2*(-x0[1] + x1[1])*(-((-x0[0] + x1[0])*(-x0[0] + x3[0])) - 
         (-x0[1] + x1[1])*(-x0[1] + x3[1]) - (-x0[2] + x1[2])*(-x0[2] + x3[2])))/
     pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
       pow(-x0[2] + x1[2],2),2) - 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(-x0[1] + x3[1]) + 
         (-x0[1] + x2[1])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) - 
         2*(-x0[1] + x1[1])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) - 
    (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(-x0[1] + x3[1])) - 
       (-x0[1] + x2[1])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
       2*(-x0[1] + x1[1])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2])))/
     (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2))) + 
    ((2*(-x0[1] + x2[1])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])) - 
         2*(-x0[1] + x1[1])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
            ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x3[2]))) - 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2)),2) + 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       (2*(-x0[1] + x2[1])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])) - 
         2*(-x0[1] + x1[1])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)),2)) + 
    ((-x0[1] + x2[1])*((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) - 
    (2*(-x0[1] + x1[1])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
         (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     (pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
         pow(-x0[2] + x1[2],2),2)*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   (-x0[2] + x3[2])/
     (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))\
     - (2*(-x0[2] + x1[2])*(-((-x0[0] + x1[0])*(-x0[0] + x3[0])) - 
         (-x0[1] + x1[1])*(-x0[1] + x3[1]) - (-x0[2] + x1[2])*(-x0[2] + x3[2])))/
     pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
       pow(-x0[2] + x1[2],2),2) - 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(-x0[2] + x3[2]) + 
         (-x0[2] + x2[2])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) - 
         2*(-x0[2] + x1[2])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) - 
    (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*(-x0[2] + x3[2])) - 
       (-x0[2] + x2[2])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
       2*(-x0[2] + x1[2])*((-x0[0] + x2[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x2[1])*(-x0[1] + x3[1]) + (-x0[2] + x2[2])*(-x0[2] + x3[2])))/
     (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2))) + 
    ((2*(-x0[2] + x2[2])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])) - 
         2*(-x0[2] + x1[2])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
            ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x3[2]))) - 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2)),2) + 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       (2*(-x0[2] + x2[2])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])) - 
         2*(-x0[2] + x1[2])*(-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)),2)) + 
    ((-x0[2] + x2[2])*((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) - 
    (2*(-x0[2] + x1[2])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
         (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     (pow(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
         pow(-x0[2] + x1[2],2),2)*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   -((((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
           (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
         ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
              pow(-x0[2] + x1[2],2))*(-x0[0] + x3[0]) - 
           (x0[0] - x1[0])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
              (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2]))))/
       ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
           pow(-x0[2] + x1[2],2))*
         (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
              (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
           (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
              pow(-x0[2] + x1[2],2))*
            (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
              pow(-x0[2] + x2[2],2))))) - 
    (-((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[0] + x3[0])) + 
       (x0[0] - x1[0])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
          (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])))/
     (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2))) + 
    ((-2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[0] + x2[0]) - 
         2*(x0[0] - x1[0])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])))*
       (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
            ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x3[2]))) - 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2)),2) + 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       (-2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[0] + x2[0]) - 
         2*(x0[0] - x1[0])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)),2)) - 
    ((x0[0] - x1[0])*((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   -((-((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
              pow(-x0[2] + x1[2],2))*(-x0[1] + x3[1])) - 
         (-x0[1] + x1[1])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])))/
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) - 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[1] + x3[1]) + 
         (-x0[1] + x1[1])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) + 
    ((-2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[1] + x2[1]) + 
         2*(-x0[1] + x1[1])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])))*
       (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
            ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x3[2]))) - 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2)),2) + 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       (-2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[1] + x2[1]) + 
         2*(-x0[1] + x1[1])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)),2)) + 
    ((-x0[1] + x1[1])*((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   -((-((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
              pow(-x0[2] + x1[2],2))*(-x0[2] + x3[2])) - 
         (-x0[2] + x1[2])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2])))/
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) - 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[2] + x3[2]) + 
         (-x0[2] + x1[2])*((-x0[0] + x1[0])*(-x0[0] + x3[0]) + 
            (-x0[1] + x1[1])*(-x0[1] + x3[1]) + (-x0[2] + x1[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))) + 
    ((-2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[2] + x2[2]) + 
         2*(-x0[2] + x1[2])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])))*
       (-((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
            ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
              (-x0[2] + x1[2])*(-x0[2] + x3[2]))) - 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2)),2) + 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       (-2*(-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[2] + x2[2]) + 
         2*(-x0[2] + x1[2])*((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2])))*
       ((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + (-x0[1] + x1[1])*(-x0[1] + x2[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       pow(-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) - (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)),2)) + 
    ((-x0[2] + x1[2])*((-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
          ((-x0[0] + x1[0])*(-x0[0] + x3[0]) + (-x0[1] + x1[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x1[2])*(-x0[2] + x3[2])) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          ((-x0[0] + x2[0])*(-x0[0] + x3[0]) + (-x0[1] + x2[1])*(-x0[1] + x3[1]) + 
            (-x0[2] + x2[2])*(-x0[2] + x3[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   (-x0[0] + x1[0])/
     (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))\
     - (-((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[0] + x2[0])) - 
       (-x0[0] + x1[0])*(-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
          (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2])))/
     (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2))) - 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[0] + x2[0]) + 
         (-x0[0] + x1[0])*(-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   (-x0[1] + x1[1])/
     (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))\
     - (-((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[1] + x2[1])) - 
       (-x0[1] + x1[1])*(-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
          (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2])))/
     (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2))) - 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[1] + x2[1]) + 
         (-x0[1] + x1[1])*(-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2)))),
   (-x0[2] + x1[2])/
     (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))\
     - (-((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[2] + x2[2])) - 
       (-x0[2] + x1[2])*(-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
          (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2])))/
     (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
          (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
       (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
          pow(-x0[2] + x1[2],2))*
        (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
          pow(-x0[2] + x2[2],2))) - 
    (((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
         (-x0[2] + x1[2])*(-x0[2] + x2[2]))*
       ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*(-x0[2] + x2[2]) + 
         (-x0[2] + x1[2])*(-((x0[0] - x1[0])*(-x0[0] + x2[0])) + 
            (-x0[1] + x1[1])*(-x0[1] + x2[1]) + (-x0[2] + x1[2])*(-x0[2] + x2[2]))))/
     ((-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - pow(-x0[2] + x1[2],2))*
       (-pow((x0[0] - x1[0])*(-x0[0] + x2[0]) - (-x0[1] + x1[1])*(-x0[1] + x2[1]) - 
            (-x0[2] + x1[2])*(-x0[2] + x2[2]),2) + 
         (-pow(-x0[0] + x1[0],2) - pow(-x0[1] + x1[1],2) - 
            pow(-x0[2] + x1[2],2))*
          (-pow(-x0[0] + x2[0],2) - pow(-x0[1] + x2[1],2) - 
            pow(-x0[2] + x2[2],2))))};

	wGrad[0][0] = Vector3d(w[0], w[1], w[2]);
	wGrad[0][1] = Vector3d(w[3], w[4], w[5]);
	wGrad[0][2] = Vector3d(w[6], w[7], w[8]);
	wGrad[0][3] = Vector3d(w[9], w[10], w[11]);
	
	wGrad[1][0] = Vector3d(u[0], u[1], u[2]);
	wGrad[1][1] = Vector3d(u[3], u[4], u[5]);
	wGrad[1][2] = Vector3d(u[6], u[7], u[8]);
	wGrad[1][3] = Vector3d(u[9], u[10], u[11]);
	
	wGrad[2][0] = Vector3d(v[0], v[1], v[2]);
	wGrad[2][1] = Vector3d(v[3], v[4], v[5]);
	wGrad[2][2] = Vector3d(v[6], v[7], v[8]);
	wGrad[2][3] = Vector3d(v[9], v[10], v[11]);

	wGrad[3][0] = Vector3d::Zero();
	wGrad[3][1] = Vector3d::Zero();
	wGrad[3][2] = Vector3d::Zero();
	wGrad[3][3] = Vector3d::Zero();
}

int PrismVolume::areCoplanar(Vector3d& x0, Vector3d& x1, Vector3d& x2, Vector3d& x3,
							 Vector3d& x0v, Vector3d& x1v, Vector3d& x2v, Vector3d& x3v, double *dt)
{
	Vector3d x10 = x1  - x0;
	Vector3d x20 = x2  - x0;
	Vector3d x30 = x3  - x0;
	Vector3d v10 = x1v - x0v;
	Vector3d v20 = x2v - x0v;
	Vector3d v30 = x3v - x0v;
	
	Vector3d p = x10.cross(x20);
	Vector3d q = v10.cross(x20) + x10.cross(v20);
	Vector3d r = v10.cross(v20);

	double a = v30.dot(r);
	double b = x30.dot(r) + v30.dot(q);
	double c = x30.dot(q) + v30.dot(p);
	double d = x30.dot(p);

	// The polynomial coefficients
	//
	//double a = -v21[2]*v31[1]*v41[0] + v21[1]*v31[2]*v41[0] + v21[2]*v31[0]*v41[1] -
	//            v21[0]*v31[2]*v41[1] - v21[1]*v31[0]*v41[2] + v21[0]*v31[1]*v41[2];
	//
	//double b = -v31[2]*v41[1]*x21[0] + v31[1]*v41[2]*x21[0] + v31[2]*v41[0]*x21[1] -
	//            v31[0]*v41[2]*x21[1] - v31[1]*v41[0]*x21[2] + v31[0]*v41[1]*x21[2] +
	//            v21[2]*v41[1]*x31[0] - v21[1]*v41[2]*x31[0] - v21[2]*v41[0]*x31[1] +
	//            v21[0]*v41[2]*x31[1] + v21[1]*v41[0]*x31[2] - v21[0]*v41[1]*x31[2] -
	//            v21[2]*v31[1]*x41[0] + v21[1]*v31[2]*x41[0] + v21[2]*v31[0]*x41[1] -
	//            v21[0]*v31[2]*x41[1] - v21[1]*v31[0]*x41[2] + v21[0]*v31[1]*x41[2];
	//
	//double c = -v41[2]*x21[1]*x31[0] + v41[1]*x21[2]*x31[0] + v41[2]*x21[0]*x31[1] -
	//            v41[0]*x21[2]*x31[1] - v41[1]*x21[0]*x31[2] + v41[0]*x21[1]*x31[2] +
	//            v31[2]*x21[1]*x41[0] - v31[1]*x21[2]*x41[0] - v31[2]*x21[0]*x41[1] +
	//            v31[0]*x21[2]*x41[1] + v21[2]*x31[0]*x41[1] - v21[0]*x31[2]*x41[1] +
	//            v31[1]*x21[0]*x41[2] - v31[0]*x21[1]*x41[2] - v21[1]*x31[0]*x41[2] -
	//            v21[2]*x31[1]*x41[0] + v21[1]*x31[2]*x41[0] + v21[0]*x31[1]*x41[2];
	//
	//double d = -x21[2]*x31[1]*x41[0] + x21[1]*x31[2]*x41[0] + x21[2]*x31[0]*x41[1] -
	//            x21[0]*x31[2]*x41[1] - x21[1]*x31[0]*x41[2] + x21[0]*x31[1]*x41[2];
	
	int leadCoeff = 3;
	double coeffs[4];
	if (a == 0.0)
	{
		leadCoeff = 2;
		if (b == 0.0)
		{
			leadCoeff = 1;
			if (c == 0.0)
			{
				// Degenerate polynomial
				//
				return 0;
			}
			else
			{
				coeffs[0] = c;
				coeffs[1] = d;
			}
		}
		else
		{
			coeffs[0] = b;
			coeffs[1] = c;
			coeffs[2] = d;
		}
	}
	else
	{
		coeffs[0] = a;
		coeffs[1] = b;
		coeffs[2] = c;
		coeffs[3] = d;
	}

	// Use of Sturm's theorem, not really needed for cubics...
	//
	//double seq1[4];
	//double seq2[4];

	//seq1[0] = coeffs[3];
	//seq2[0] = coeffs[0] + coeffs[1] + coeffs[2] + coeffs[3];

	//seq1[1] = coeffs[2];
	//seq2[1] = 3.0 * coeffs[0] + 2.0 * coeffs[1] + coeffs[2];
	//for (size_t i=2; i<4; ++i)
	//{
	//	double nbr1 = int(seq1[i-2] / seq1[i-1]);
	//	double nbr2 = int(seq2[i-2] / seq2[i-1]);

	//	seq1[i] = -(seq1[i-2] - (seq1[i-1] * nbr1));
	//	seq2[i] = -(seq2[i-2] - (seq2[i-1] * nbr2));
	//}

	//int nbrChanges1 = 0;
	//int nbrChanges2 = 0;
	//for (size_t i=1; i<4; ++i)
	//{
	//	if ((seq1[i-1] * seq1[i]) <= 0.0)
	//		nbrChanges1++;
	//	
	//	if ((seq2[i-1] * seq2[i]) <= 0.0)
	//		nbrChanges2++;
	//}

	//if ((nbrChanges1 - nbrChanges2) == 0)
		//return 0;

	//cout << coeffs[0] << " " << coeffs[1] << " " << coeffs[2] << " " << coeffs[3] << " " << leadCoeff << endl;
	
	// Find the roots
	//
	//int count = cubicRoots(a, b, c, d, timeStep, dt);
	RootFinder rf;
	double rl[3];
	double im[3];
	int cnt = rf.rpoly(coeffs, leadCoeff, rl, im);
	int count = 0;
	if (cnt != -1)
	{
		for (size_t i=0; i<cnt; ++i)
		{
			if (im[i] == 0.0)
				dt[count++] = rl[i];
		}
	}

	return count;
}

bool PrismVolume::VertexFaceProximity(Vector3d& x0, Vector3d& x1, Vector3d& x2, Vector3d& x3,
									  Vector3d& v0, Vector3d& v1, Vector3d& v2, Vector3d& v3,
									  double h, double& tRet)
{
	//double t[7];
	//int count = areProximate(x0, x1, x2, x3, v0, v1, v2, v3, h, t);


        double t_total[30];
        int count_total = 0;
        double t_tmp[7];
        int count_tmp = 0;
        
	count_tmp = areProximate(x0, x1, x2, x3, v0, v1, v2, v3, h, t_tmp);
        for (size_t i=0; i<count_tmp; ++i)
        {
          if(t_tmp[i]>=0 && t_tmp[i]<=1)
          {
            t_total[count_total] = t_tmp[i];
            count_total++;
          }
        }

        
        count_tmp = arePointPointInt(x0,x3,v0,v3,h,t_tmp);
        for (size_t i=0; i<count_tmp; ++i)
        {
          if(t_tmp[i]>=0 && t_tmp[i]<=1)
          {
            t_total[count_total] = t_tmp[i];
            count_total++;
          }
        }
        count_tmp = arePointPointInt(x1,x3,v1,v3,h,t_tmp);
        for (size_t i=0; i<count_tmp; ++i)
        {
          if(t_tmp[i]>=0 && t_tmp[i]<=1)
          {
            t_total[count_total] = t_tmp[i];
            count_total++;
          }
        }
        count_tmp = arePointPointInt(x2,x3,v2,v3,h,t_tmp);
        for (size_t i=0; i<count_tmp; ++i)
        {
          if(t_tmp[i]>=0 && t_tmp[i]<=1)
          {
            t_total[count_total] = t_tmp[i];
            count_total++;
          }
        }

	count_tmp = arePointLineInt(x0, x1, x3, v0, v1, v3, h, t_tmp);
        for (size_t i=0; i<count_tmp; ++i)
        {
          if(t_tmp[i]>=0 && t_tmp[i]<=1)
          {
            t_total[count_total] = t_tmp[i];
            count_total++;
          }
        }
        count_tmp = arePointLineInt(x1, x2, x3, v1, v2, v3, h, t_tmp);
        for (size_t i=0; i<count_tmp; ++i)
        {
          if(t_tmp[i]>=0 && t_tmp[i]<=1)
          {
            t_total[count_total] = t_tmp[i];
            count_total++;
          }
        }
        count_tmp = arePointLineInt(x2, x0, x3, v2, v0, v3, h, t_tmp);
        for (size_t i=0; i<count_tmp; ++i)
        {
          if(t_tmp[i]>=0 && t_tmp[i]<=1)
          {
            t_total[count_total] = t_tmp[i];
            count_total++;
          }
        }
        
        assert(count_total <= 30); 
        if(count_total > 30)
            std::cout<<"v-f num of roots: "<<count_total<<std::endl;

        bool flag = false;
        tRet = 2;
	for (int k=0; k<count_total; ++k)
	{
               if(t_total[k] == 0)
                   std::cout<<"v-f root is zero????\n";

		if (t_total[k] > 0.0 && t_total[k] <= 1.0)
		{
			Vector3d p0 = x0 + t_total[k] * v0;
			Vector3d p1 = x1 + t_total[k] * v1;
			Vector3d p2 = x2 + t_total[k] * v2;
			Vector3d p3 = x3 + t_total[k] * v3;
			
			double t1, t2, t3;
			double distance = vertexFaceDist(p3, p0, p1, p2, t1, t2, t3);

			Vector3d relVel = v3 - (t1 * v0 + t2 * v1 + t3 * v2);
			Vector3d e1 = p1 - p0;
			Vector3d e2 = p2 - p0;
			Vector3d n = e1.cross(e2);
					
			double Vn = relVel.dot(n);
			//double thiseps = Vn * Vn * 1e-8;
			//double thiseps = h * 1e-8;
			double thiseps = 1e-9;

			if (std::fabs(distance - h*h) < thiseps)
			{
				//tRet = max(t[k]-0.01, 0.0);
				//tRet = max(t[k]-0.02, 0.0);
                                if(t_total[k] < tRet)
				  tRet = t_total[k];

                                flag = true;
				//return true;
			}
		}
	}
	
	//return false;
        return flag;
}

int PrismVolume::areProximate(Vector3d& x0, Vector3d& x1, Vector3d& x2, Vector3d& x3,
							  Vector3d& v0, Vector3d& v1, Vector3d& v2, Vector3d& v3, double h, double *dt)
{
	Vector3d x10 = x1 - x0;
	Vector3d x20 = x2 - x0;
	Vector3d x30 = x3 - x0;
	Vector3d v10 = v1 - v0;
	Vector3d v20 = v2 - v0;
	Vector3d v30 = v3 - v0;
	
	Vector3d p = x10.cross(x20);
	Vector3d q = v10.cross(x20) + x10.cross(v20);
	Vector3d r = v10.cross(v20);

	double a = x30.dot(p);
	double b = x30.dot(q) + v30.dot(p);
	double c = x30.dot(r) + v30.dot(q);
	double d = v30.dot(r);
	
	// Coefficients are in order of decreasing powers (coeffs[0] is t^6 term)
	//
	double coeffs[7];
	coeffs[6] = ((a * a - p.dot(p) * h * h)); //
	coeffs[5] = ((2.0 * a * b - (2.0 * p.dot(q)) * h * h)); //
	coeffs[4] = ((b * b + 2.0 * a * c - (q.dot(q) + 2.0 * p.dot(r)) * h * h)); //
	coeffs[3] = ((2.0 * a * d + 2.0 * b * c - (2.0 * q.dot(r)) * h * h)); //
	coeffs[2] = ((c * c + 2.0 * b * d - r.dot(r) * h * h)); //
	coeffs[1] = ((2.0 * c * d)); //
	coeffs[0] = (d * d); //

	int leadCoeff = 6;
	double *co = &coeffs[0];
	while (coeffs[6-leadCoeff] == 0.0)
	{
		leadCoeff--;
		co++;
	}

	//double poly = Options::getRealValue("polycheck", 1);
        double poly = 1;

	// Fast polynomial checks
	if (poly > 0.001)
	{
//		cerr << "poly! " << poly << endl;
		double fa[7];
		double a[7];
		int signsum  = 0;
		int signzero = 0;

		for(int i=0; i<7; ++i)
		{
			a[i] = coeffs[6-i];
			if (a[i] == 0)
			{
				signzero++;
			}
			else
			{
				if (a[i]>=0)
					signsum++;
			    else
					signsum--;
			}
			
			fa[i] = fabs(a[i]);
		}
//#define CHECKPLOT 
		
#ifdef CHECKPLOT
		cerr << "signsum = " << signsum << endl;
		for(int i=0; i<7; ++i)
		{
			cerr << a[i] << endl;
		}
		cerr << endl;
		for(int i=0; i<7; ++i)
		{
			cerr << fa[i] << endl;
		}
		cerr << endl;
		
#endif
		
		// all coefficients have the same sign
		if (signsum+signzero == 7 || signsum+signzero == -7)
		{
#ifdef CHECKPLOT
			cerr << "C0" << endl;
#endif
			return 0;
		}
		
		if (fa[0] > fa[6] + fa[5] + fa[4] + fa[3] + fa[2] + fa[1])
		{
#ifdef CHECKPLOT
			cerr << "C1" << endl;
#endif
			return 0;
		}
		
		if (
			 (a[0]*a[1] < 0) && 
			 (fa[1] > 2*fa[2] + 3*fa[3] + 4*fa[4] + 5*fa[5] +6*fa[6]) && 
			 (a[0]*(a[0]+a[1]+a[2]+a[3])) > 0
		   )
		{
#ifdef CHECKPLOT
			cerr << "C2" << endl;
#endif
			return 0;
		}
		
		if (
			(a[0] > 0) && 
			(a[1] > 0) && 
			(a[0]+a[1] > fa[6]+fa[5]+fa[4]+fa[3]+fa[2])
		   )
		{
#ifdef CHECKPLOT
			cerr << "C3" << endl;
#endif
			return 0;
		}
		
		if (
			(a[0] > 0) && 
			(a[1] > 0) && 
			(a[2] > 0) && 
			(a[0]+a[1]+a[2] > fa[6]+fa[5]+fa[4]+fa[3])
			)
		{
#ifdef CHECKPLOT
			cerr << "C4" << endl;
#endif
		return 0;
		}

		if (
			(a[0] > 0) && 
			(a[1] > 0) && 
			(a[2] > 0) && 
			(a[3] > 0) && 
			(a[0]+a[1]+a[2]+a[3] > fa[6]+fa[5]+fa[4])
			)
		{
#ifdef CHECKPLOT
			cerr << "C5" << endl;
#endif
			return 0;
		}
		
		if (
			(a[0] > 0) && 
			(a[1] > 0) && 
			(a[2] > 0) && 
			(a[3] > 0) && 
			(a[4] > 0) && 
			(a[0]+a[1]+a[2]+a[3]+a[4] > fa[6]+fa[5])
			)
		{
#ifdef CHECKPLOT
			cerr << "C5" << endl;
#endif
			return 0;
		}
		
		if (
			(a[0] > 0) && 
			(a[1] > 0) && 
			(a[2] > 0) && 
			(a[3] > 0) && 
			(a[4] > 0) && 
			(a[5] > 0) && 
			(a[0]+a[1]+a[2]+a[3]+a[4]+a[5] > fa[6])
			)
		{
#ifdef CHECKPLOT
			cerr << "C6" << endl;
#endif
			return 0;
		}
		
		if (
			(a[0] < 0) && 
			(a[1] < 0) && 
			(a[0]+a[1] < -(fa[6]+fa[5]+fa[4]+fa[3]+fa[2]))
			)
		{
#ifdef CHECKPLOT
			cerr << "C7" << endl;
#endif
			return 0;
		}
		
		if (
			(a[0] < 0) && 
			(a[1] < 0) && 
			(a[2] < 0) && 
			(a[0]+a[1]+a[2] < -(fa[6]+fa[5]+fa[4]+fa[3]))
			)
		{
#ifdef CHECKPLOT
			cerr << "C8" << endl;
#endif
			return 0;
		}

		if (
			(a[0] < 0) && 
			(a[1] < 0) && 
			(a[2] < 0) && 
			(a[3] < 0) && 
			(a[0]+a[1]+a[2]+a[3] < -(fa[6]+fa[5]+fa[4]))
			)
		{
#ifdef CHECKPLOT
			cerr << "C9" << endl;
#endif
		return 0;
		}

		if (
			(a[0] < 0) && 
			(a[1] < 0) && 
			(a[2] < 0) && 
			(a[3] < 0) && 
			(a[4] < 0) && 
			(a[0]+a[1]+a[2]+a[3]+a[4] < -(fa[6]+fa[5]))
			)
		{
#ifdef CHECKPLOT
			cerr << "C10" << endl;
#endif
			return 0;
		}

		if (
			(a[0] < 0) && 
			(a[1] < 0) && 
			(a[2] < 0) && 
			(a[3] < 0) && 
			(a[4] < 0) && 
			(a[5] < 0) && 
			(a[0]+a[1]+a[2]+a[3]+a[4]+a[5] < -fa[6])
			)
		{
#ifdef CHECKPLOT
			cerr << "C11" << endl;
#endif
			return 0;
		}
		
	}
	
	
	
	//std::cout << coeffs[0] << " " << coeffs[1] << " " << coeffs[2] << " " << coeffs[3] << " "
	//		  << coeffs[4] << " " << coeffs[5] << " " << coeffs[6] << " " << leadCoeff << std::endl;
	
	// Find the roots
	//
	RootFinder rf;
	double rl[6];
	double im[6];
	int cnt = rf.rpoly(co, leadCoeff, rl, im);
	int count = 0;
	if (cnt != -1)
	{
		for (size_t i=0; i<cnt; ++i)
		{
			if (im[i] == 0.0)
				dt[count++] = rl[i];
		}
	}
	
	return count;
}

bool PrismVolume::EdgeEdgeProximity(Vector3d& x0, Vector3d& x1, Vector3d& x2, Vector3d& x3,
									  Vector3d& v0, Vector3d& v1, Vector3d& v2, Vector3d& v3,
									  double h, double& tRet)
{
        bool flag = false;
        tRet = 2;
	//double t[20];
	double t[7];
	int count = areEEProximate(x0, x1, x2, x3, v0, v1, v2, v3, h, t);

        assert(count <= 7);
        if(count > 7)
            std::cout<<"e-e num of roots: "<<count<<std::endl;
	
	for (int k=0; k<count; ++k)
	{
                if(t[k] == 0)
                    std::cout<<"e-e root is zero???\n";

		if (t[k] > 0.0 && t[k] <= 1.0)
		{
                        Vector3d p0 = x0 + t[k] * v0;
			Vector3d p1 = x1 + t[k] * v1;
			Vector3d p2 = x2 + t[k] * v2;
			Vector3d p3 = x3 + t[k] * v3;
			
			double r, s;
			double distance = edgeEdgeDist(p0, p1, p2, p3, r, s);
			
			Vector3d relVel = (v0 * (1.0 - r) + v1 * r) - (v2 * (1.0 - s) + v3 * s);
			Vector3d e1 = p1 - p0;
			Vector3d e2 = p3 - p2;
			Vector3d n  = e1.cross(e2);

			double Vn = relVel.dot(n);
			//double thiseps = Vn * Vn * 1e-8;
			//double thiseps = h * 1e-8;
			double thiseps = 1e-9;

			if (std::fabs(distance - h*h) < thiseps)
			{
				//tRet = max(t[k]-0.01, 0.0);
				//tRet = max(t[k]-0.02, 0.0);
                                if(t[k] < tRet)
				  tRet = t[k];
                                flag = true;
				//return true;
			}
		}
	}
	
        return flag;
	//return false;
	
}

int PrismVolume::areEEProximate(Eigen::Vector3d& x0, Eigen::Vector3d& x1, Eigen::Vector3d& x2, Eigen::Vector3d& x3,
							Eigen::Vector3d& v0, Eigen::Vector3d& v1, Eigen::Vector3d& v2, Eigen::Vector3d& v3, double h, double *dt)
{
  int count = areEdgeEdgeInt(x0,x1,x2,x3,v0,v1,v2,v3,h,dt);
  return count;
}

double PrismVolume::vertexFaceDist(Vector3d& p, Vector3d& a, Vector3d& b, Vector3d& c, double &t1, double &t2, double &t3)
{
	double ab[3], ac[3], ap[3], bp[3];
	
	ab[0] = b[0] - a[0];
	ab[1] = b[1] - a[1];
	ab[2] = b[2] - a[2];
	
	ac[0] = c[0] - a[0];
	ac[1] = c[1] - a[1];
	ac[2] = c[2] - a[2];
	
	ap[0] = p[0] - a[0];
	ap[1] = p[1] - a[1];
	ap[2] = p[2] - a[2];
	
	double d1 = ab[0]*ap[0] + ab[1]*ap[1] + ab[2]*ap[2];
	double d2 = ac[0]*ap[0] + ac[1]*ap[1] + ac[2]*ap[2];
	
	if ((d1 <= 0.0f) && (d2 <= 0.0f))
	{
		t1 = 1.0f;
		t2 = 0.0f;
		t3 = 0.0f;
		
		return ((p[0]-a[0])*(p[0]-a[0]) + (p[1]-a[1])*(p[1]-a[1]) + (p[2]-a[2])*(p[2]-a[2]));
	}
	
	bp[0] = p[0] - b[0];
	bp[1] = p[1] - b[1];
	bp[2] = p[2] - b[2];
	
	double d3 = ab[0]*bp[0] + ab[1]*bp[1] + ab[2]*bp[2];
	double d4 = ac[0]*bp[0] + ac[1]*bp[1] + ac[2]*bp[2];
	
	if ((d3 >= 0.0f) && (d4 <= d3))
	{
		t1 = 0.0f;
		t2 = 1.0f;
		t3 = 0.0f;
		
		return ((p[0]-b[0])*(p[0]-b[0]) + (p[1]-b[1])*(p[1]-b[1]) + (p[2]-b[2])*(p[2]-b[2]));
	}
	
	double vc = d1*d4 - d3*d2;
	
	if ((vc <= 0.0f) && (d1 >= 0.0f) && (d3 <= 0.0f))
	{
		double v = d1 / (d1 - d3);
		
		t1 = 1-v;
		t2 = v;
		t3 = 0;
		
		double vec[3];
		vec[0] = p[0] - (a[0]+v*ab[0]);
		vec[1] = p[1] - (a[1]+v*ab[1]);
		vec[2] = p[2] - (a[2]+v*ab[2]);
		
		return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	}
	
	double cp[3];
	cp[0] = p[0] - c[0];
	cp[1] = p[1] - c[1];
	cp[2] = p[2] - c[2];
	
	double d5 = ab[0]*cp[0] + ab[1]*cp[1] + ab[2]*cp[2];
	double d6 = ac[0]*cp[0] + ac[1]*cp[1] + ac[2]*cp[2];
	
	if ((d6 >= 0.0f) && (d5 <= d6))
	{
		t1 = 0;
		t2 = 0;
		t3 = 1;
		
		return ((p[0]-c[0])*(p[0]-c[0]) + (p[1]-c[1])*(p[1]-c[1]) + (p[2]-c[2])*(p[2]-c[2]));
	}
	
	double vb = d5*d2 - d1*d6;
	
	if ((vb <= 0.0f) && (d2 >= 0.0f) && (d6 <= 0.0f))
	{
		double w = d2 / (d2 - d6);
		
		t1 = 1-w;
		t2 = 0;
		t3 = w;
		
		double vec[3];
		vec[0] = p[0] - (a[0]+w*ac[0]);
		vec[1] = p[1] - (a[1]+w*ac[1]);
		vec[2] = p[2] - (a[2]+w*ac[2]);
		
		return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	}
	
	double va = d3*d6 - d5*d4;
	
	if ((va <= 0.0f) && ((d4-d3) >= 0.0f) && ((d5-d6) >= 0.0f))
	{
		double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		
		t1 = 0;
		t2 = 1-w;
		t3 = w;
		
		double vec[3];
		vec[0] = p[0] - (b[0]+w*(c[0]-b[0]));
		vec[1] = p[1] - (b[1]+w*(c[1]-b[1]));
		vec[2] = p[2] - (b[2]+w*(c[2]-b[2]));
		
		return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	}
	
	double denom = 1.0f / (va + vb + vc);
	double v = vb * denom;
	double w = vc * denom;
	double u = 1.0 - v - w;
	
	t1 = u;
	t2 = v;
	t3 = w;
	
	double vec[3];
	vec[0] = p[0] - (u*a[0] + v*b[0] + w*c[0]);
	vec[1] = p[1] - (u*a[1] + v*b[1] + w*c[1]);
	vec[2] = p[2] - (u*a[2] + v*b[2] + w*c[2]);
	
	return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

double PrismVolume::edgeEdgeDist(Vector3d &p1, Vector3d &q1,
								 Vector3d &p2, Vector3d &q2, double &s, double &t)
{
	double d1[3], d2[3], r[3], a, e, f;
	double c1[3], c2[3];

	d1[0] = q1[0] - p1[0];
	d1[1] = q1[1] - p1[1];
	d1[2] = q1[2] - p1[2];

	d2[0] = q2[0] - p2[0];
	d2[1] = q2[1] - p2[1];
	d2[2] = q2[2] - p2[2];

	r[0] = p1[0] - p2[0];
	r[1] = p1[1] - p2[1];
	r[2] = p1[2] - p2[2];

	a = d1[0]*d1[0] + d1[1]*d1[1] + d1[2]*d1[2];
	e = d2[0]*d2[0] + d2[1]*d2[1] + d2[2]*d2[2];
	f = d2[0]*r[0] + d2[1]*r[1] + d2[2]*r[2];

	// check if either or both segments degenerate into points
	//
	if ((a <= 1e-8) && (e <= 1e-8))
	{
		s = t = 0.0f;
		c1[0] = p1[0]; c1[1] = p1[1]; c1[2] = p1[2];
		c2[0] = p2[0]; c2[1] = p2[1]; c2[2] = p2[2];

		return ((c1[0]-c2[0])*(c1[0]-c2[0]) + (c1[1]-c2[1])*(c1[1]-c2[1]) + (c1[2]-c2[2])*(c1[2]-c2[2]));
	}

	if (a <= 1e-8)
	{
		// first segment degenerates into a point
		//
		s = 0.0f;
		t = f / e;
		if (t<0.0f) t = 0.0f;
		if (t>1.0f) t = 1.0f;
	}
	else
	{
		double c = d1[0]*r[0] + d1[1]*r[1] + d1[2]*r[2];

		if (e <= 1e-8)
		{
			// second segment degenerates into a point
			//
			t = 0.0f;
			s = -c / a;
			if (s<0.0f) s = 0.0f;
			if (s>1.0f) s = 1.0f;
		}
		else
		{
			// nondegenerate case
			//
			double b = d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2];
			double denom = a*e - b*b;

			if (denom != 0.0f)
			{
				s = (b*f - c*e) / denom;
				if (s<0.0f) s = 0.0f;
				if (s>1.0f) s = 1.0f;
			}
			else
				s = 0.0f;

			double tnom = b*s + f;
			if (tnom < 0.0f)
			{
				t = 0.0f;
				s = -c / a;
				if (s<0.0f) s = 0.0f;
				if (s>1.0f) s = 1.0f;
			}
			else if (tnom > e)
			{
				t = 1.0f;
				s = (b - c) / a;
				if (s<0.0f) s = 0.0f;
				if (s>1.0f) s = 1.0f;
			}
			else
				t = tnom / e;
		}
	}

        //cout<<"edge-edge dist, s: "<<s<<" --- t: "<<t<<endl;

	c1[0] = p1[0] + d1[0] * s;
	c1[1] = p1[1] + d1[1] * s;
	c1[2] = p1[2] + d1[2] * s;

	c2[0] = p2[0] + d2[0] * t;
	c2[1] = p2[1] + d2[1] * t;
	c2[2] = p2[2] + d2[2] * t;

	return ((c1[0]-c2[0])*(c1[0]-c2[0]) + (c1[1]-c2[1])*(c1[1]-c2[1]) + (c1[2]-c2[2])*(c1[2]-c2[2]));
}

int PrismVolume::areLineLineInt(Vector3d& x0, Vector3d& x1, Vector3d& x2, Vector3d& x3,
    Vector3d& v0, Vector3d& v1, Vector3d& v2, Vector3d& v3, double h, double *dt)
{
  
  Vector3d x10 = x1 - x0;
  Vector3d x20 = x3 - x2;
  Vector3d x30 = x3 - x0;
  Vector3d v10 = v1 - v0;
  Vector3d v20 = v3 - v2;
  Vector3d v30 = v3 - v0;

  Vector3d p = x10.cross(x20);
  Vector3d q = v10.cross(x20) + x10.cross(v20);
  Vector3d r = v10.cross(v20);

  double a = x30.dot(p);
  double b = x30.dot(q) + v30.dot(p);
  double c = x30.dot(r) + v30.dot(q);
  double d = v30.dot(r);
    
  // Coefficients are in order of decreasing powers (coeffs[0] is t^6 term)
  //
  double coeffs[7];
  coeffs[6] = ((a * a - p.dot(p) * h * h)); //
  coeffs[5] = ((2.0 * a * b - (2.0 * p.dot(q)) * h * h)); //
  coeffs[4] = ((b * b + 2.0 * a * c - (q.dot(q) + 2.0 * p.dot(r)) * h * h)); //
  coeffs[3] = ((2.0 * a * d + 2.0 * b * c - (2.0 * q.dot(r)) * h * h)); //
  coeffs[2] = ((c * c + 2.0 * b * d - r.dot(r) * h * h)); //
  coeffs[1] = ((2.0 * c * d)); //
  coeffs[0] = (d * d); //
  
  // check leading coeff
  //
  int leadCoeff = 6;
  double *co = &coeffs[0];
  while (coeffs[6-leadCoeff] == 0.0)
  {
    leadCoeff--;
    co++;
  }
  if(leadCoeff == -1)
    return 0;
  
  // TODO: fast polynomial root checks

  // Find the roots
  //
  RootFinder rf;
  double rl[6];
  double im[6];
  int cnt = rf.rpoly(co, leadCoeff, rl, im);
  int count = 0;
  if (cnt != -1)
  {
    for (size_t i=0; i<cnt; ++i)
    {
      if (im[i] == 0.0)
        dt[count++] = rl[i];
    }
  }
  
  return count;
}

int PrismVolume::arePointPointInt(Vector3d& x0, Vector3d& x1,
    Vector3d& v0, Vector3d& v1, double h, double *dt)
{
  Vector3d x10 = x1 - x0;
  Vector3d v10 = v1 - v0;
  
  double coeffs[3];
  coeffs[0] = v10.dot(v10);
  coeffs[1] = 2*v10.dot(x10);
  coeffs[2] = x10.dot(x10) - h*h;
  
  // check leading coeff
  //
  int leadCoeff = 2;
  double *co = &coeffs[0];
  while (coeffs[2-leadCoeff] == 0.0)
  {
    leadCoeff--;
    co++;
  }
  if(leadCoeff == -1)
    return 0;
  
  // TODO: fast polynomial root checks
  
  // Find the roots
  //
  RootFinder rf;
  double rl[2];
  double im[2];
  int cnt = rf.rpoly(co, leadCoeff, rl, im);
  int count = 0;
  if (cnt != -1)
  {
    for (size_t i=0; i<cnt; ++i)
    {
      if (im[i] == 0.0)
        dt[count++] = rl[i];
    }
  }
  
  return count;
}

int PrismVolume::arePointLineInt(Vector3d& x0, Vector3d& x1, Vector3d& x2,
    Vector3d& v0, Vector3d& v1, Vector3d& v2, double h, double *dt)
{
  Vector3d x10 = x1 - x0;
  Vector3d x20 = x2 - x0;
  Vector3d v10 = v1 - v0;
  Vector3d v20 = v2 - v0;
  
  Vector3d p = v20.cross(v10);
  Vector3d q = v20.cross(x10) + x20.cross(v10);
  Vector3d r = x20.cross(x10);

  double coeffs[5];
  coeffs[0] = p.dot(p);
  coeffs[1] = 2*p.dot(q);
  coeffs[2] = 2*p.dot(r) + q.dot(q) - v10.dot(v10)*h*h; 
  coeffs[3] = 2*q.dot(r) - 2*v10.dot(x10)*h*h; 
  coeffs[4] = r.dot(r) - x10.dot(x10)*h*h;
  
  // check leading coeff
  //
  int leadCoeff = 4;
  double *co = &coeffs[0];
  while (coeffs[4-leadCoeff] == 0.0)
  {
    leadCoeff--;
    co++;
  }
  if(leadCoeff == -1)
    return 0;
  
  // TODO: fast polynomial root checks
  
  // Find the roots
  //
  RootFinder rf;
  double rl[4];
  double im[4];
  int cnt = rf.rpoly(co, leadCoeff, rl, im);
  int count = 0;
  if (cnt != -1)
  {
    for (size_t i=0; i<cnt; ++i)
    {
      if (im[i] == 0.0)
        dt[count++] = rl[i];
    }
  }
  
  return count;
}

int PrismVolume::areEdgeEdgeInt(Vector3d& x0, Vector3d& x1, Vector3d& x2, Vector3d& x3,
    Vector3d& v0, Vector3d& v1, Vector3d& v2, Vector3d& v3, double h, double *dt)
{
  //cout<<"begin EE Prox"<<endl; 
  int count = 0;

  int count_tmp = 0;
  double dt_tmp[7];

  
  /*
  count_tmp = arePointPointInt(x0,x2,v0,v2,h,dt_tmp);
  for (size_t i=0; i<count_tmp; ++i)
  {
    if(dt_tmp[i]>=0 && dt_tmp[i]<=1)
    {
      dt[count] = dt_tmp[i];
      //std::cout<<"count_tmp: "<<count_tmp<<"---"<<dt_tmp[i]<<std::endl;
      count++;
    }
  }
  count_tmp = arePointPointInt(x0,x3,v0,v3,h,dt_tmp);
  for (size_t i=0; i<count_tmp; ++i)
  {
    if(dt_tmp[i]>=0 && dt_tmp[i]<=1)
    {
      dt[count] = dt_tmp[i];
      //std::cout<<"count_tmp: "<<count_tmp<<"---"<<dt_tmp[i]<<std::endl;
      count++;
    }
  }
  count_tmp = arePointPointInt(x1,x2,v1,v2,h,dt_tmp);
  for (size_t i=0; i<count_tmp; ++i)
  {
    if(dt_tmp[i]>=0 && dt_tmp[i]<=1)
    {
      dt[count] = dt_tmp[i];
      //std::cout<<"count_tmp: "<<count_tmp<<"---"<<dt_tmp[i]<<std::endl;
      count++;
    }
  }
  count_tmp = arePointPointInt(x1,x3,v1,v3,h,dt_tmp);
  for (size_t i=0; i<count_tmp; ++i)
  {
    if(dt_tmp[i]>=0 && dt_tmp[i]<=1)
    {
      dt[count] = dt_tmp[i];
      //std::cout<<"count_tmp: "<<count_tmp<<"---"<<dt_tmp[i]<<std::endl;
      count++;
    }
  }
  
  count_tmp = arePointLineInt(x0,x1,x2,v0,v1,v2,h,dt_tmp);
  for (size_t i=0; i<count_tmp; ++i)
  {
    if(dt_tmp[i]>=0 && dt_tmp[i]<=1)
    {
      dt[count] = dt_tmp[i];
      //std::cout<<"count_tmp: "<<count_tmp<<"---"<<dt_tmp[i]<<std::endl;
      count++; 
    }
  }
  count_tmp = arePointLineInt(x0,x1,x3,v0,v1,v3,h,dt_tmp);
  for (size_t i=0; i<count_tmp; ++i)
  {
    if(dt_tmp[i]>=0 && dt_tmp[i]<=1)
    {
      dt[count] = dt_tmp[i];
      //std::cout<<"count_tmp: "<<count_tmp<<"---"<<dt_tmp[i]<<std::endl;
      count++; 
    }
  }
  count_tmp = arePointLineInt(x2,x3,x0,v2,v3,v0,h,dt_tmp);
  for (size_t i=0; i<count_tmp; ++i)
  {
    if(dt_tmp[i]>=0 && dt_tmp[i]<=1)
    {
      dt[count] = dt_tmp[i];
      //std::cout<<"count_tmp: "<<count_tmp<<"---"<<dt_tmp[i]<<std::endl;
      count++; 
    }
  }
  count_tmp = arePointLineInt(x2,x3,x1,v2,v3,v1,h,dt_tmp);
  for (size_t i=0; i<count_tmp; ++i)
  {
    if(dt_tmp[i]>=0 && dt_tmp[i]<=1)
    {
      dt[count] = dt_tmp[i];
      //std::cout<<"count_tmp: "<<count_tmp<<"---"<<dt_tmp[i]<<std::endl;
      count++; 
    }
  }
  */

  count_tmp = areLineLineInt(x0,x1,x2,x3,v0,v1,v2,v3,h,dt_tmp);
  for (size_t i=0; i<count_tmp; ++i)
  {
    if(dt_tmp[i]>=0 && dt_tmp[i]<=1)
    {
      dt[count] = dt_tmp[i];
      //std::cout<<"count_tmp: "<<count_tmp<<"---"<<dt_tmp[i]<<std::endl;
      count++; 
    }
  }

  //cout<<"end EE Prox: "<<count<<endl; 
  //for (size_t i=0; i<count; ++i)
    //cout<<"t["<<i<<"]: "<<dt[i]<<endl;
  return count;
}

void PrismVolume::SegmentAABB(Vector3d &v1s, Vector3d &v1e, Vector3d &v2s, Vector3d &v2e,
                 Vector3d &emin, Vector3d &emax, double h)
{
    Point3d p[4];
    p[0] = v1s;
    p[1] = v1e;
    p[2] = v2s;
    p[3] = v2e;
		
    for(unsigned j = 0; j < 3; ++j) 
    {
        emin[j] = p[0][j];
        emax[j] = p[0][j];
    }
		
    for(unsigned j = 1; j < 4; ++j)
    {
        for(unsigned k = 0; k < 3; ++k) 
        {
            if (p[j][k] < emin[k])
                emin[k] = p[j][k];
            if (p[j][k] > emax[k])
                emax[k] = p[j][k];
        }
    }
		
    for (unsigned j=0; j<3; ++j)
    {
        emin[j] -= h/2; emin[j] -= 1e-10; 
        emax[j] += h/2; emax[j] += 1e-10;
    }
}

bool PrismVolume::TestSegmentSegment(Vector3d &e1v1s, Vector3d &e1v1e, Vector3d &e1v2s, Vector3d &e1v2e,
                        Vector3d &e2v1s, Vector3d &e2v1e, Vector3d &e2v2s, Vector3d &e2v2e, double h)
{
    Vector3d e1min, e1max, e2min, e2max;
    SegmentAABB(e1v1s, e1v1e, e1v2s, e1v2e, e1min, e1max, h);
    SegmentAABB(e2v1s, e2v1e, e2v2s, e2v2e, e2min, e2max, h);
    
    for (size_t i=0; i<3; ++i)
    {
        if (e1min[i] > e2max[i] || e1max[i] < e2min[i])
            return false;
    }

    return true;
}

bool PrismVolume::TestSegmentSubdivideTime(Vector3d &e1v1s, Vector3d &e1v1e, Vector3d &e1v2s, Vector3d &e1v2e,
                              Vector3d &e2v1s, Vector3d &e2v1e, Vector3d &e2v2s, Vector3d &e2v2e, double h)
{
    // time mid point m
    Vector3d e1v1m = (e1v1s + e1v1e)/2; Vector3d e1v2m = (e1v2s + e1v2e)/2;
    Vector3d e2v1m = (e2v1s + e2v1e)/2; Vector3d e2v2m = (e2v2s + e2v2e)/2;

    return( TestSegmentSegment(e1v1s, e1v1m, e1v2s, e1v2m, e2v1s, e2v1m, e2v2s, e2v2m, h) || 
            TestSegmentSegment(e1v1m, e1v1e, e1v2m, e1v2e, e2v1m, e2v1e, e2v2m, e2v2e, h) );
}

bool PrismVolume::TestSegmentSubdivideSpace(Vector3d &e1v1s, Vector3d &e1v1e, Vector3d &e1v2s, Vector3d &e1v2e,
                               Vector3d &e2v1s, Vector3d &e2v1e, Vector3d &e2v2s, Vector3d &e2v2e, double h)
{
    // space mid point v3
    Vector3d e1v3s = (e1v1s + e1v2s)/2; Vector3d e1v3e = (e1v1e + e1v2e)/2;
    Vector3d e2v3s = (e2v1s + e2v2s)/2; Vector3d e2v3e = (e2v1e + e2v2e)/2;

    return(TestSegmentSubdivideTime(e1v1s, e1v1e, e1v3s, e1v3e, e2v1s, e2v1e, e2v3s, e2v3e, h) ||
           TestSegmentSubdivideTime(e1v1s, e1v1e, e1v3s, e1v3e, e2v3s, e2v3e, e2v2s, e2v2e, h) ||
           TestSegmentSubdivideTime(e1v3s, e1v3e, e1v2s, e1v2e, e2v1s, e2v1e, e2v3s, e2v3e, h) ||
           TestSegmentSubdivideTime(e1v3s, e1v3e, e1v2s, e1v2e, e2v3s, e2v3e, e2v2s, e2v2e, h) );
}

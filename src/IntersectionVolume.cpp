// IntersectionVolume.cpp
//

#include "IntersectionVolume.h"
#include "PrismVolume.h"
#include "ColSurface.h"

#include <boost/foreach.hpp>
#include <Eigen/Geometry>
#include <unsupported/Eigen/SparseExtra>

using std::set;
using std::multimap;
using std::vector;
using std::make_pair;

using namespace Eigen;



IntersectionVolume::IntersectionVolume(Polyhedron &mesh, double h, Collisions &collisions, double periodicLength, int nv)
 : m_thickness(h), m_periodicLength(periodicLength), m_nv(nv)
{
        m_bdv = 0;
	m_x_cm = Vector3d::Zero();

	CollisionMap sortedCollisions;
	findFirstCollisions(mesh, collisions, sortedCollisions);

	computeVolumeNew(mesh, sortedCollisions);
}

void IntersectionVolume::findFirstCollisionsApprox(Polyhedron &mesh, Collisions &collisions, CollisionMap &sorted)
{
    CollisionMapIterator cmItr;
    for (CollisionsIterator ci=collisions.begin(); ci!=collisions.end(); ++ci)
    {
        size_t idx[4];
        if(ci->isVF())
        {
            idx[0] = ci->getTriangleVertex1();
            idx[1] = ci->getTriangleVertex2();
            idx[2] = ci->getTriangleVertex3();
            idx[3] = ci->getVertex();
        }
        else
        {
            idx[0] = ci->getFirstEdgeVertex1();
            idx[1] = ci->getFirstEdgeVertex2(); 
            idx[2] = ci->getSecondEdgeVertex1();
            idx[3] = ci->getSecondEdgeVertex2();
        }
        Vector3d r[4];
        for (size_t i=0; i<4; ++i)
            r[i] = mesh.getVertex(idx[i])->Pori() + ci->getTime() * (mesh.getVertex(idx[i])->P() - mesh.getVertex(idx[i])->Pori());
        double weights[4];
        PrismVolume::getBarycentricCoordinates(r[0], r[1], r[2], r[3], ci->isVF(), weights);
        for (size_t i=0; i<4; ++i)
        {
            cmItr = sorted.find(idx[i]);
            if (cmItr == sorted.end()) // No collision for this vertex exists
            {
                if (weights[i] > 0)
                    sorted.insert(make_pair(idx[i], ci));
            }
            else if (ci->getTime() < cmItr->second->getTime()) // One exists, but this one is sooner
            {
                if (weights[i] > 0)
                {
                    sorted.erase(cmItr);
                    sorted.insert(make_pair(idx[i], ci));
                }
            }
        }
    }
}

void IntersectionVolume::findFirstCollisions(Polyhedron &mesh, Collisions &collisions, CollisionMap &sorted)
{
  /*
	for (CollisionsIterator ci=collisions.begin(); ci!=collisions.end(); ++ci)
	{
		//cout << ci->getVertex() << " (" << ci->getTriangleVertex1() << " " << ci->getTriangleVertex2() << " " << ci->getTriangleVertex3() << ")" << "; " << ci->getTime() << endl;
		// Check vertex
		//
		CollisionMapIterator cmItr;
		cmItr = sorted.find(ci->getVertex());
		//if (cmItr == sorted.end() || (ci->getTime() < cmItr->second->getTime()))
			sorted.insert(make_pair(ci->getVertex(), ci));

		// Triangle vertices
		//
		cmItr = sorted.find(ci->getTriangleVertex1());
		//if (cmItr == sorted.end() || (ci->getTime() < cmItr->second->getTime()))
			sorted.insert(make_pair(ci->getTriangleVertex1(), ci));
		
		cmItr = sorted.find(ci->getTriangleVertex2());
		//if (cmItr == sorted.end() || (ci->getTime() < cmItr->second->getTime()))
			sorted.insert(make_pair(ci->getTriangleVertex2(), ci));
		
		cmItr = sorted.find(ci->getTriangleVertex3());
		//if (cmItr == sorted.end() || (ci->getTime() < cmItr->second->getTime()))
			sorted.insert(make_pair(ci->getTriangleVertex3(), ci));
	}
  */
	
	CollisionMapIterator cmItr;
	for (CollisionsIterator ci=collisions.begin(); ci!=collisions.end(); ++ci)
	{
          if(ci->isVF())
          {
		size_t idx[4] = { ci->getTriangleVertex1(), ci->getTriangleVertex2(), ci->getTriangleVertex3(), ci->getVertex() };
		Vector3d r[4];
		for (size_t i=0; i<4; ++i)
			r[i] = mesh.getVertex(idx[i])->Pori() + ci->getTime() * (mesh.getVertex(idx[i])->P() - mesh.getVertex(idx[i])->Pori());

		double weights[4];
		PrismVolume::getBarycentricCoordinates(r[0], r[1], r[2], r[3], ci->isVF(), weights);
		
		for (size_t i=0; i<4; ++i)
		{
			cmItr = sorted.find(idx[i]);
			if (cmItr == sorted.end()) // No collision for this vertex exists
			{
				bool good = true;
				for (size_t j=0; j<3; ++j)
				//for (size_t j=0; j<0; ++j)
				{
					if (i != j && weights[i] < weights[j])
					{
						good = false;
						break;
					}
				}
				//if (good)
				if (good || i==3)
					sorted.insert(make_pair(idx[i], ci));
			}
			else if (ci->getTime() < cmItr->second->getTime()) // One exists, but this one is sooner
			{
				bool good = true;
				for (size_t j=0; j<3; ++j)
				//for (size_t j=0; j<0; ++j)
				{
					if (i != j && weights[i] < weights[j])
					{
						good = false;
						break;
					}
				}
				//if (good)
				if (good || i==3)
				{
					sorted.erase(cmItr);
					sorted.insert(make_pair(idx[i], ci));
				}
			}
		}
          }
          else
          {
                size_t idx[4] = {ci->getFirstEdgeVertex1(), ci->getFirstEdgeVertex2(), ci->getSecondEdgeVertex1(), ci->getSecondEdgeVertex2()};
		Vector3d r[4];
		for (size_t i=0; i<4; ++i)
			r[i] = mesh.getVertex(idx[i])->Pori() + ci->getTime() * (mesh.getVertex(idx[i])->P() - mesh.getVertex(idx[i])->Pori());

		double weights[4];
		PrismVolume::getBarycentricCoordinates(r[0], r[1], r[2], r[3], ci->isVF(), weights);
                
                // first edge
                for (size_t i=0; i<2; ++i)
		{
			cmItr = sorted.find(idx[i]);
			if (cmItr == sorted.end()) // No collision for this vertex exists
			{
				bool good = true;
				for (size_t j=0; j<2; ++j)
				//for (size_t j=0; j<0; ++j)
				{
					if (i != j && weights[i] < weights[j])
					{
						good = false;
						break;
					}
				}
				if (good)
					sorted.insert(make_pair(idx[i], ci));
			}
			else if (ci->getTime() < cmItr->second->getTime()) // One exists, but this one is sooner
			{
				bool good = true;
				for (size_t j=0; j<2; ++j)
				//for (size_t j=0; j<0; ++j)
				{
					if (i != j && weights[i] < weights[j])
					{
						good = false;
						break;
					}
				}
				if (good)
				{
					sorted.erase(cmItr);
					sorted.insert(make_pair(idx[i], ci));
				}
			}
		}
                
                //second edge
                for (size_t i=2; i<4; ++i)
		{
			cmItr = sorted.find(idx[i]);
			if (cmItr == sorted.end()) // No collision for this vertex exists
			{
				bool good = true;
				for (size_t j=2; j<4; ++j)
				//for (size_t j=2; j<2; ++j)
				{
					if (i != j && weights[i] < weights[j])
					{
						good = false;
						break;
					}
				}
				if (good)
					sorted.insert(make_pair(idx[i], ci));
			}
			else if (ci->getTime() < cmItr->second->getTime()) // One exists, but this one is sooner
			{
				bool good = true;
				for (size_t j=2; j<4; ++j)
				//for (size_t j=2; j<2; ++j)
				{
					if (i != j && weights[i] < weights[j])
					{
						good = false;
						break;
					}
				}
				if (good)
				{
					sorted.erase(cmItr);
					sorted.insert(make_pair(idx[i], ci));
				}
			}
		}
          }
	}
}

double IntersectionVolume::getNormal(Polyhedron &mesh, Collision &c, unsigned int index[4], unsigned int vIndex, Vector3d &n)
{
	double l = 0.0;
	n = Vector3d::Zero();

	// 2 cases, Vertex-Triangle and Edge-Edge
	//
	if (c.isVF())
	{
		// 2 cases, Vertex-Triangle or Triangle-Vertex
		//
		if (vIndex == 3) // Vertex-Triangle, compute triangle normal
		{
			Vertex *vh0 = mesh.v_array[index[0]];
			Vertex *vh1 = mesh.v_array[index[1]];
			Vertex *vh2 = mesh.v_array[index[2]];
				
			Vector3d f0 = vh0->Pori() + c.getTime() * (vh0->P() - vh0->Pori());
			Vector3d f1 = vh1->Pori() + c.getTime() * (vh1->P() - vh1->Pori());
			Vector3d f2 = vh2->Pori() + c.getTime() * (vh2->P() - vh2->Pori());

			n = (f1 - f0).cross(f2 - f0);
		}
		else // Triangle-Vertex, fetch vertex normal
		{
			std::vector<std::pair<int, Halfedge *> > vp = mesh.get1RingH(index[vIndex]);
			for (std::vector<std::pair<int, Halfedge *> >::iterator p=vp.begin(); p!=vp.end(); ++p)
			{
				if (p->second->is_border())
					continue;

				Vertex *vh0 = mesh.v_array[p->second->facet()->vi[0]];
				Vertex *vh1 = mesh.v_array[p->second->facet()->vi[1]];
				Vertex *vh2 = mesh.v_array[p->second->facet()->vi[2]];

				Vector3d f0 = vh0->Pori() + c.getTime() * (vh0->P() - vh0->Pori());
				Vector3d f1 = vh1->Pori() + c.getTime() * (vh1->P() - vh1->Pori());
				Vector3d f2 = vh2->Pori() + c.getTime() * (vh2->P() - vh2->Pori());

				n += (f1 - f0).cross(f2 - f0);
			}
		}
	}
	else
	{
		unsigned int vtx1;
		unsigned int vtx2;
		if (vIndex == 0 || vIndex == 1)
		{
			vtx1 = 2;
			vtx2 = 3;
		}
		else
		{
			vtx1 = 0;
			vtx2 = 1;
		}

		Halfedge_Vertex_circulator he = mesh.getVertex(index[vtx1])->vertex_begin();
		do
		{
			if (he->next()->vertex()->getIndex() == index[vtx2])
				break;
			++he;
		} while (he != mesh.getVertex(index[vtx1])->vertex_begin());

		Vertex *vh0 = mesh.v_array[he->facet()->vi[0]];
		Vertex *vh1 = mesh.v_array[he->facet()->vi[1]];
		Vertex *vh2 = mesh.v_array[he->facet()->vi[2]];

		Vector3d f0 = vh0->Pori() + c.getTime() * (vh0->P() - vh0->Pori());
		Vector3d f1 = vh1->Pori() + c.getTime() * (vh1->P() - vh1->Pori());
		Vector3d f2 = vh2->Pori() + c.getTime() * (vh2->P() - vh2->Pori());

		n += (f1 - f0).cross(f2 - f0);
		
		if (!he->is_border())
		{
			vh0 = mesh.v_array[he->opposite()->facet()->vi[0]];
			vh1 = mesh.v_array[he->opposite()->facet()->vi[1]];
			vh2 = mesh.v_array[he->opposite()->facet()->vi[2]];

			f0 = vh0->Pori() + c.getTime() * (vh0->P() - vh0->Pori());
			f1 = vh1->Pori() + c.getTime() * (vh1->P() - vh1->Pori());
			f2 = vh2->Pori() + c.getTime() * (vh2->P() - vh2->Pori());

			n += (f1 - f0).cross(f2 - f0);
		}
	}
			
	l = n.norm();
	n.normalize();

	return l;
}

void IntersectionVolume::getNormalGradient(Vector3d &r0, Vector3d &r1, Vector3d &r2, bool flip, Matrix3d N[3])
{
	Vector3d e[3];
	e[0] = r2 - r1;
	e[1] = r0 - r2;
	e[2] = r1 - r0;

	for (size_t i=0; i<3; ++i)
	{
		N[i] = Matrix3d::Zero();
		N[i](0,1) = -e[i](2);
		N[i](0,2) =  e[i](1);
		N[i](1,0) =  e[i](2);
		N[i](1,2) = -e[i](0);
		N[i](2,0) = -e[i](1);
		N[i](2,1) =  e[i](0);
	}

//	Vector3d n = (r1 - r0).cross(r2 - r0);
//
//	double l = n.norm();
//	n /= l;
//
//	if (flip)
//		n = -n;
//
//	for (size_t i=0; i<3; ++i)
//		N[i] = ((e[i].cross(n)) / l) * n.transpose();
}

void IntersectionVolume::getNormalGradient(Polyhedron &mesh, Collision &c, unsigned int index[4], unsigned int vIndex, bool flip, Vector3d &n, double l, Vector3d &nGrad)
{
	nGrad = Vector3d::Zero();
	
	Matrix3d dN = (1.0 / l) * (Matrix3d::Identity() - n * n.transpose());

	// 2 cases, Vertex-Triangle and Edge-Edge
	//
	if (c.isVF())
	{
		// 2 cases, Vertex-Triangle or Triangle-Vertex
		//
		if (vIndex == 3) // Vertex-Triangle, compute triangle normal
		{
			Vertex *vh0 = mesh.v_array[index[0]];
			Vertex *vh1 = mesh.v_array[index[1]];
			Vertex *vh2 = mesh.v_array[index[2]];
			
			Vector3d f0 = vh0->Pori() + c.getTime() * (vh0->P() - vh0->Pori());
			Vector3d f1 = vh1->Pori() + c.getTime() * (vh1->P() - vh1->Pori());
			Vector3d f2 = vh2->Pori() + c.getTime() * (vh2->P() - vh2->Pori());
		
			Matrix3d N[3];
			getNormalGradient(f0, f1, f2, flip, N);

			for (size_t i=0; i<3; ++i)
				N[i] = dN * N[i];

			nGrad += N[0] * (vh0->P() - vh0->Pori());
			nGrad += N[1] * (vh1->P() - vh1->Pori());
			nGrad += N[2] * (vh2->P() - vh2->Pori());
		}
		else // Triangle-Vertex, fetch vertex normal
		{
			std::vector<std::pair<int, Halfedge *> > vp = mesh.get1RingH(index[vIndex]);
			for (std::vector<std::pair<int, Halfedge *> >::iterator p=vp.begin(); p!=vp.end(); ++p)
			{
				if (p->second->is_border())
					continue;

				Vertex *vh0 = mesh.v_array[p->second->facet()->vi[0]];
				Vertex *vh1 = mesh.v_array[p->second->facet()->vi[1]];
				Vertex *vh2 = mesh.v_array[p->second->facet()->vi[2]];
			
				Vector3d f0 = vh0->Pori() + c.getTime() * (vh0->P() - vh0->Pori());
				Vector3d f1 = vh1->Pori() + c.getTime() * (vh1->P() - vh1->Pori());
				Vector3d f2 = vh2->Pori() + c.getTime() * (vh2->P() - vh2->Pori());

				Matrix3d N[3];
				getNormalGradient(f0, f1, f2, flip, N);
			
				for (size_t i=0; i<3; ++i)
					N[i] = dN * N[i];
			
				nGrad += N[0] * (vh0->P() - vh0->Pori());
				nGrad += N[1] * (vh1->P() - vh1->Pori());
				nGrad += N[2] * (vh2->P() - vh2->Pori());
			}
		}
	}
	else
	{
		unsigned int vtx1;
		unsigned int vtx2;
		if (vIndex == 0 || vIndex == 1)
		{
			vtx1 = 2;
			vtx2 = 3;
		}
		else
		{
			vtx1 = 0;
			vtx2 = 1;
		}

		Halfedge_Vertex_circulator he = mesh.getVertex(index[vtx1])->vertex_begin();
		do
		{
			if (he->next()->vertex()->getIndex() == index[vtx2])
				break;
			++he;
		} while (he != mesh.getVertex(index[vtx1])->vertex_begin());

		Vertex *vh0 = mesh.v_array[he->facet()->vi[0]];
		Vertex *vh1 = mesh.v_array[he->facet()->vi[1]];
		Vertex *vh2 = mesh.v_array[he->facet()->vi[2]];

		Vector3d f0 = vh0->Pori() + c.getTime() * (vh0->P() - vh0->Pori());
		Vector3d f1 = vh1->Pori() + c.getTime() * (vh1->P() - vh1->Pori());
		Vector3d f2 = vh2->Pori() + c.getTime() * (vh2->P() - vh2->Pori());

		Matrix3d N[3];
		getNormalGradient(f0, f1, f2, flip, N);

		nGrad += N[0] * (vh0->P() - vh0->Pori());
		nGrad += N[1] * (vh1->P() - vh1->Pori());
		nGrad += N[2] * (vh2->P() - vh2->Pori());
		
		if (!he->is_border())
		{
			vh0 = mesh.v_array[he->opposite()->facet()->vi[0]];
			vh1 = mesh.v_array[he->opposite()->facet()->vi[1]];
			vh2 = mesh.v_array[he->opposite()->facet()->vi[2]];

			f0 = vh0->Pori() + c.getTime() * (vh0->P() - vh0->Pori());
			f1 = vh1->Pori() + c.getTime() * (vh1->P() - vh1->Pori());
			f2 = vh2->Pori() + c.getTime() * (vh2->P() - vh2->Pori());

			getNormalGradient(f0, f1, f2, flip, N);

			nGrad += N[0] * (vh0->P() - vh0->Pori());
			nGrad += N[1] * (vh1->P() - vh1->Pori());
			nGrad += N[2] * (vh2->P() - vh2->Pori());
		}
	}
}

/*
void IntersectionVolume::computeVolume(Polyhedron &mesh, CollisionMap &collisions)
{
	DynamicSparseMatrix<double> gradients(3 * mesh.size_of_vertices(), 1);

	m_volume = 0.0;
	for (CollisionMapIterator cmItr=collisions.begin(); cmItr!=collisions.end(); ++cmItr)
	{
		CollisionsIterator ci = cmItr->second;
		unsigned int index[4] = { ci->getTriangleVertex1(), ci->getTriangleVertex2(), ci->getTriangleVertex3(), ci->getVertex() };

		Vector3d x[4], v[4];
		for (size_t i=0; i<4; ++i)
		{
			x[i] = mesh.getVertex(index[i])->Pori();
			v[i] = mesh.getVertex(index[i])->P() - mesh.getVertex(index[i])->Pori();
		}

		unsigned int whichIndex = 3;
		for (size_t i=0; i<3; ++i)
		{
			if (cmItr->first == index[i])
			{
				whichIndex = i;
				break;
			}
		}

		Vector3d n;
		double l = getNormal(mesh, *ci, index, whichIndex, n);

		// Compute the positions at the moment of intersection
		//
		Vector3d r[4];
		for (size_t i=0; i<4; ++i)
			r[i] = x[i] + ci->getTime() * v[i];

		double weights[4];
		PrismVolume::getBarycentricCoordinates(r[0], r[1], r[2], r[3], ci->isVF(), weights);

		double A = getVoronoiArea(mesh, index[whichIndex]);
		double vn = weights[whichIndex] * v[whichIndex].dot(n);
		double V = (1.0 - ci->getTime()) * vn * A;

		double flip = 1.0;
		if (V > 0.0)
		{
			m_volume -= V;
			n = -n;
			vn = -vn;
			flip = -1.0;
		}
		else if (V < 0.0)
		{
			m_volume += V;
		}
		else // Volume is 0.0, insufficient information to determine flip
		{
			// Use relative velocity in normal direction to determine "correct" normal orientation,
			// since gradients may still be non-zero, and we need a consistent direction
			//
			if (ci->isVF())
			{
				double relVel = (v[3] - (weights[0] * v[0] + weights[1] * v[1] + weights[2] * v[2])).dot(n);
				if ((relVel > 0.0 && whichIndex == 3) || (relVel < 0.0 && whichIndex != 3))
				{
					n = -n;
					vn = -vn;
					flip = -1.0;
				}
			}
			else
			{
				double relVel = ((v[2] * weights[2] + v[3] * weights[3]) - (v[0] * weights[0] + v[1] * weights[1])).dot(n);
				if ((relVel > 0.0 && (whichIndex == 2 || whichIndex == 3)) || (relVel < 0.0 && (whichIndex == 0 || whichIndex == 1)))
				{
					n = -n;
					vn = -vn;
					flip = -1.0;
				}
			}
		}

		Vector3d tGrad[4];
		if (ci->isVF())
			PrismVolume::getTimeGradient(x[0], x[1], x[2], x[3], v[0], v[1], v[2], v[3], weights, ci->getTime(), m_thickness, tGrad);
		else
			PrismVolume::getEETimeGradient(x[0], x[1], x[2], x[3], v[0], v[1], v[2], v[3], weights, ci->getTime(), m_thickness, tGrad);

		Vector3d nGrad;
		getNormalGradient(mesh, *ci, index, whichIndex, (flip < 0.0), n, l, nGrad);

		// Use finite difference for the barycentric coordinate gradient
		//
		double wGrad = 0.0;
		{
			static const double eps = 1.0e-5;

			for (size_t j=0; j<4; ++j)
				r[j] = x[j] + (ci->getTime() + eps) * v[j];

			double wgts[4];
			PrismVolume::getBarycentricCoordinates(r[0], r[1], r[2], r[3], ci->isVF(), wgts);
			wGrad = (wgts[whichIndex] - weights[whichIndex]) / eps;
		}

		for (size_t i=0; i<4; ++i)
		{
			double g = ((1.0 - ci->getTime()) * (wGrad * v[whichIndex].transpose() * n + weights[whichIndex] * v[whichIndex].dot(nGrad)) - vn);

			for (size_t j=0; j<3; ++j)
				gradients.coeffRef(3*index[i]+j, 0) += A * (g * tGrad[i][j] + (i == whichIndex) * (1.0 - ci->getTime()) * weights[whichIndex] * n[j]);
		}
	}

	//m_volume *= 1.15;

	gradients.finalize();
	m_gradients = SparseVector<double>(gradients);
}
*/

/*
void IntersectionVolume::computeVolume(Polyhedron &mesh, CollisionMap &collisions)
{
	DynamicSparseMatrix<double> gradients(3 * mesh.size_of_vertices(), 1);

	m_volume = 0.0;
	for (CollisionMapIterator cmItr=collisions.begin(); cmItr!=collisions.end(); ++cmItr)
	{
		CollisionsIterator ci = cmItr->second;

		if (ci->isVF())
		{
			unsigned int index[4] = { ci->getTriangleVertex1(), ci->getTriangleVertex2(), ci->getTriangleVertex3(), ci->getVertex() };

			Vector3d x[4], v[4];
			for (size_t i=0; i<4; ++i)
			{
				x[i] = mesh.getVertex(index[i])->Pori();
				v[i] = mesh.getVertex(index[i])->P() - mesh.getVertex(index[i])->Pori();
			}

                        if(m_periodicLength > 0)
                        {
                            Vector3d transV(0,0,0);
                            PrismVolume::getTranslationVector(transV, 
                                    mesh.getVertex(index[3]), mesh.getVertex(index[0]), m_periodicLength);
                            x[3] = x[3] + transV;
                        }

			unsigned int whichIndex = 3;
			for (size_t i=0; i<3; ++i)
			{
				if (cmItr->first == index[i])
				{
					whichIndex = i;
					break;
				}
			}
                        
                        //double A = getVoronoiArea(mesh, index[whichIndex]);
			double A = 1.0;
			double time = 1.0;
			double V1 = getVolume(x[0], x[1], x[2], x[3], v[0], v[1], v[2], v[3], v[whichIndex], A, time);

			// Compute the positions at the moment of intersection
                        //
			Vector3d r[4];
			for (size_t i=0; i<4; ++i)
				r[i] = x[i] + time * v[i];

			double weights[4];
			PrismVolume::getBarycentricCoordinates(r[0], r[1], r[2], r[3], true, weights);
                        assert(weights[0]>=0 && weights[0] <= 1);
                        assert(weights[1]>=0 && weights[1] <= 1);
                        assert(weights[2]>=0 && weights[2] <= 1);
                        assert(weights[3]>=0 && weights[3] <= 1);

                        //std::cout<<"this is vf contact? "<<ci->isVF()<<". --- multipling weight: "
                        //<<weights[whichIndex]<<std::endl; 
			//V1 *= weights[whichIndex];
                        
                        //std::cout<<"this is vf contact? "<<ci->isVF()<<
                        //    "\n"<<whichIndex<<
                        //    "\n"<<weights[0]<<
                        //    "\n"<<weights[1]<<
                        //    "\n"<<weights[2]<<
                        //    "\n"<<weights[3]<<std::endl;

                        //std::cout<<"V1: "<<V1<<". time: "<<time<<std::endl;
                        //assert(V1 != 0.0);
			if (V1 > 0.0)
			{
				m_volume -= V1;
			}
			else
			{
				m_volume += V1;
			}

			for (size_t i=0; i<4; ++i)
				for (size_t j=0; j<3; ++j)
                                {
                                        if(i==whichIndex)
					  gradients.coeffRef(3 * index[i] + j, 0) += mesh.getVertex(index[i])->normal[j];//A*weights[whichIndex];
                                }
		}
		else
		{
			unsigned int index[4] = { ci->getFirstEdgeVertex1(), ci->getFirstEdgeVertex2(),
									  ci->getSecondEdgeVertex1(), ci->getSecondEdgeVertex2() };

			Vector3d x[4], v[4];
			for (size_t i=0; i<4; ++i)
			{
				x[i] = mesh.getVertex(index[i])->Pori();
				v[i] = mesh.getVertex(index[i])->P() - mesh.getVertex(index[i])->Pori();
			}

                        if(m_periodicLength > 0)
                        {
                            Vector3d transV(0,0,0);
                            PrismVolume::getTranslationVector(transV, 
                                    mesh.getVertex(index[0]), mesh.getVertex(index[2]), m_periodicLength);
                            x[0] = x[0] + transV;
                            x[1] = x[1] + transV;
                        }

			unsigned int whichIndex = 3;
			for (size_t i=0; i<3; ++i)
			{
				if (cmItr->first == index[i])
				{
					whichIndex = i;
					break;
				}
			}
                        
                        //double A = getVoronoiArea(mesh, index[whichIndex]);
			double A = 1.0;
			double time = 1.0;
			double V1 = getEEVolume(x[0], x[1], x[2], x[3], v[0], v[1], v[2], v[3], v[whichIndex], A, time);

			// Compute the positions at the moment of intersection
			//
			Vector3d r[4];
			for (size_t i=0; i<4; ++i)
				r[i] = x[i] + time * v[i];

			// TODO
			double weights[4];
			PrismVolume::getBarycentricCoordinates(r[0], r[1], r[2], r[3], false, weights);
                        assert(weights[0]>=0 && weights[0] <= 1);
                        assert(weights[1]>=0 && weights[1] <= 1);
                        assert(weights[2]>=0 && weights[2] <= 1);
                        assert(weights[3]>=0 && weights[3] <= 1);
			//weights[0] = (1.0 - weights[1]); weights[2] = (1.0 - weights[3]);

                        //std::cout<<"this is vf contact? "<<ci->isVF()<<". --- multipling weight: "
                        //<<weights[whichIndex]<<std::endl; 
			//V1 *= weights[whichIndex];
                        
                        //std::cout<<"this is vf contact? "<<ci->isVF()<<
                        //    "\n"<<whichIndex<<
                        //    "\n"<<weights[0]<<
                        //    "\n"<<weights[1]<<
                        //    "\n"<<weights[2]<<
                        //    "\n"<<weights[3]<<std::endl;

                        //assert(V1 != 0.0);
			if (V1 > 0.0)
			{
				m_volume -= V1;
			}
			else
			{
				m_volume += V1;
			}

			for (size_t i=0; i<4; ++i)
				for (size_t j=0; j<3; ++j)
                                {
                                        if(i==whichIndex)
					  gradients.coeffRef(3 * index[i] + j, 0) += mesh.getVertex(index[i])->normal[j];//A*weights[whichIndex];
                                }
		}
	}

	m_volume *= 1.00;
        m_volume = min(m_volume,-1e-10);

	gradients.finalize();
	m_gradients = SparseVector<double>(gradients);
}
*/

void IntersectionVolume::computeVolumeApprox(Polyhedron &mesh, CollisionMap &collisions)
{
	DynamicSparseMatrix<double> gradients(3 * mesh.size_of_vertices(), 1);
	m_volume = 0.0;
	for (CollisionMapIterator cmItr=collisions.begin(); cmItr!=collisions.end(); ++cmItr)
	{
		CollisionsIterator ci = cmItr->second;
		if (ci->isVF())
		{
			unsigned int index[4] = { ci->getTriangleVertex1(), ci->getTriangleVertex2(), ci->getTriangleVertex3(), ci->getVertex() };

			Vector3d x[4], v[4];
			for (size_t i=0; i<4; ++i)
			{
				x[i] = mesh.getVertex(index[i])->Pori();
				v[i] = mesh.getVertex(index[i])->P() - mesh.getVertex(index[i])->Pori();
			}

                        if(m_periodicLength > 0)
                        {
                            Vector3d transV(0,0,0);
                            PrismVolume::getTranslationVector(transV, 
                                    mesh.getVertex(index[3]), mesh.getVertex(index[0]), m_periodicLength);
                            x[3] = x[3] + transV;
                        }

			unsigned int whichIndex = 3;
			for (size_t i=0; i<3; ++i)
			{
				if (cmItr->first == index[i])
				{
					whichIndex = i;
					break;
				}
			}
                        
                        //double A = getVoronoiArea(mesh, index[whichIndex]);
			double A = 1.0;
			double time = 1.0;
			double V1 = getVolume(x[0], x[1], x[2], x[3], v[0], v[1], v[2], v[3], v[whichIndex], A, time);

			// Compute the positions at the moment of intersection
                        //
			Vector3d r[4];
			for (size_t i=0; i<4; ++i)
				r[i] = x[i] + time * v[i];

			double weights[4];
			PrismVolume::getBarycentricCoordinates(r[0], r[1], r[2], r[3], true, weights);
                        assert(weights[0]>=0 && weights[0] <= 1);
                        assert(weights[1]>=0 && weights[1] <= 1);
                        assert(weights[2]>=0 && weights[2] <= 1);
                        assert(weights[3]>=0 && weights[3] <= 1);

                        //std::cout<<"this is vf contact? "<<ci->isVF()<<". --- multipling weight: "
                        //<<weights[whichIndex]<<std::endl; 
			//V1 *= weights[whichIndex];
                        /*
                        std::cout<<"this is vf contact? "<<ci->isVF()<<
                            "\n"<<whichIndex<<
                            "\n"<<weights[0]<<
                            "\n"<<weights[1]<<
                            "\n"<<weights[2]<<
                            "\n"<<weights[3]<<std::endl;
                        */

                        //std::cout<<"V1: "<<V1<<". time: "<<time<<std::endl;
                        //assert(V1 != 0.0);
			if (V1 > 0.0)
			{
				m_volume -= V1;
			}
			else
			{
				m_volume += V1;
			}

			for (size_t i=0; i<4; ++i)
				for (size_t j=0; j<3; ++j)
                                {
                                        if(i==whichIndex)
					  gradients.coeffRef(3 * index[i] + j, 0) = mesh.getVertex(index[i])->normal[j]*A;//*weights[whichIndex];
                                }
		}
		else
		{
			unsigned int index[4] = { ci->getFirstEdgeVertex1(), ci->getFirstEdgeVertex2(),
									  ci->getSecondEdgeVertex1(), ci->getSecondEdgeVertex2() };

			Vector3d x[4], v[4];
			for (size_t i=0; i<4; ++i)
			{
				x[i] = mesh.getVertex(index[i])->Pori();
				v[i] = mesh.getVertex(index[i])->P() - mesh.getVertex(index[i])->Pori();
			}

                        if(m_periodicLength > 0)
                        {
                            Vector3d transV(0,0,0);
                            PrismVolume::getTranslationVector(transV, 
                                    mesh.getVertex(index[0]), mesh.getVertex(index[2]), m_periodicLength);
                            x[0] = x[0] + transV;
                            x[1] = x[1] + transV;
                        }

			unsigned int whichIndex = 3;
			for (size_t i=0; i<3; ++i)
			{
				if (cmItr->first == index[i])
				{
					whichIndex = i;
					break;
				}
			}
                        
                        //double A = getVoronoiArea(mesh, index[whichIndex]);
			double A = 1.0;
			double time = 1.0;
			double V1 = getEEVolume(x[0], x[1], x[2], x[3], v[0], v[1], v[2], v[3], v[whichIndex], A, time);

			// Compute the positions at the moment of intersection
			//
			Vector3d r[4];
			for (size_t i=0; i<4; ++i)
				r[i] = x[i] + time * v[i];

			// TODO
			double weights[4];
			PrismVolume::getBarycentricCoordinates(r[0], r[1], r[2], r[3], false, weights);
                        assert(weights[0]>=0 && weights[0] <= 1);
                        assert(weights[1]>=0 && weights[1] <= 1);
                        assert(weights[2]>=0 && weights[2] <= 1);
                        assert(weights[3]>=0 && weights[3] <= 1);
			//weights[0] = (1.0 - weights[1]); weights[2] = (1.0 - weights[3]);

                        //std::cout<<"this is vf contact? "<<ci->isVF()<<". --- multipling weight: "
                        //<<weights[whichIndex]<<std::endl; 
			//V1 *= weights[whichIndex];
                        /*
                        std::cout<<"this is vf contact? "<<ci->isVF()<<
                            "\n"<<whichIndex<<
                            "\n"<<weights[0]<<
                            "\n"<<weights[1]<<
                            "\n"<<weights[2]<<
                            "\n"<<weights[3]<<std::endl;
                        */

                        //assert(V1 != 0.0);
			if (V1 > 0.0)
			{
				m_volume -= V1;
			}
			else
			{
				m_volume += V1;
			}

			for (size_t i=0; i<4; ++i)
				for (size_t j=0; j<3; ++j)
                                {
                                        if(i==whichIndex)
					  gradients.coeffRef(3 * index[i] + j, 0) = mesh.getVertex(index[i])->normal[j]*A;//*weights[whichIndex];
                                }
		}
	}

	m_volume *= 1.01;
        m_volume = min(m_volume,-1e-10);

	gradients.finalize();
	m_gradients = SparseVector<double>(gradients);
}

void IntersectionVolume::computeVolume(Polyhedron &mesh, CollisionMap &collisions)
{
	DynamicSparseMatrix<double> gradients(3 * mesh.size_of_vertices(), 1);

	m_volume = 0.0;
	for (CollisionMapIterator cmItr=collisions.begin(); cmItr!=collisions.end(); ++cmItr)
	{
		CollisionsIterator ci = cmItr->second;

		if (ci->isVF())
		{
			unsigned int index[4] = { ci->getTriangleVertex1(), ci->getTriangleVertex2(), ci->getTriangleVertex3(), ci->getVertex() };

			Vector3d x[4], v[4];
			for (size_t i=0; i<4; ++i)
			{
				x[i] = mesh.getVertex(index[i])->Pori();
				v[i] = mesh.getVertex(index[i])->P() - mesh.getVertex(index[i])->Pori();
			}

                        if(m_periodicLength > 0)
                        {
                            Vector3d transV(0,0,0);
                            PrismVolume::getTranslationVector(transV, 
                                    mesh.getVertex(index[3]), mesh.getVertex(index[0]), m_periodicLength);
                            x[3] = x[3] + transV;
                        }

			unsigned int whichIndex = 3;
			for (size_t i=0; i<3; ++i)
			{
				if (cmItr->first == index[i])
				{
					whichIndex = i;
					break;
				}
			}

			//double A = getVoronoiArea(mesh, index[whichIndex]);
			double A = 1.0;
			double time = 1.0;
			double V1 = getVolume(x[0], x[1], x[2], x[3], v[0], v[1], v[2], v[3], v[whichIndex], A, time);

			// Compute the positions at the moment of intersection
                        //
			Vector3d r[4];
			for (size_t i=0; i<4; ++i)
				r[i] = x[i] + time * v[i];

			double weights[4];
			PrismVolume::getBarycentricCoordinates(r[0], r[1], r[2], r[3], true, weights);
                        assert(weights[0]>=0 && weights[0] <= 1);
                        assert(weights[1]>=0 && weights[1] <= 1);
                        assert(weights[2]>=0 && weights[2] <= 1);
                        assert(weights[3]>=0 && weights[3] <= 1);

                        //std::cout<<"this is vf contact? "<<ci->isVF()<<". --- multipling weight: "<<weights[whichIndex]<<std::endl; 
			V1 *= weights[whichIndex];

			Vector3d normal = (r[1] - r[0]).cross(r[2] - r[0]);
			Vector3d n(normal);
			double l = n.norm();
			n.normalize();


			// NOTE: flip should only be -1.0 if hitting the backside of a face, which should
			// only happen on a surface with boundary, maybe assert this?
			//
                        //std::cout<<"V1: "<<V1<<". time: "<<time<<std::endl;
                        //assert(V1 != 0.0);
			double flip = 1.0;
			if (V1 > 0.0)
			{
				m_volume -= V1;
                                flip = -1.0;
				//if (whichIndex == 3)
				//	flip = -1.0;
			}
			else if (V1 < 0.0)
			{
				m_volume += V1;
				//if (whichIndex != 3)
				//	flip = -1.0;
			}
			else // Volume is 0.0, insufficient information to determine flip
			{
				// Use relative velocity in normal direction to determine "correct" normal orientation,
				// since gradients may still be non-zero, and we need a consistent direction
				//
				double relVel = (v[3] - (weights[0] * v[0] + weights[1] * v[1] + weights[2] * v[2])).dot(n);
				if (relVel > 0.0)
					flip = -1.0;
			}
                        Vector3d rel_vel = v[3] - (weights[0] * v[0] + weights[1] * v[1] + weights[2] * v[2]);

			// The three edge vectors
			//
			Vector3d e[3];
			e[0] = r[2] - r[1];
			e[1] = r[0] - r[2];
			e[2] = r[1] - r[0];

			Vector3d t[4];
                        // TODO: fix analytic time gradient when vertex face contact is 
                        // degenerated to vertex-vertex or vertex-edge contact.
                        if(1)
                        //if(weights[0]<0 || weights[1]<0 || weights[2]<0 || weights[3]<0)
                        {
                          // finite difference for time gradient
                          for (size_t i=0; i<4; ++i)
                          {
                            for (size_t j=0; j<3; ++j)
                            { 
                              const double deps = 1e-8;
                              Vector3d veps[4];
                              for (size_t ii=0; ii<4; ++ii)
                              {
                                for (size_t jj=0; jj<3; ++jj)
                                {
                                  veps[ii][jj] = v[ii][jj];
                                }
                              }
                              veps[i][j] = veps[i][j] + deps;
                              double timeeps = 1.0;
                              double V1eps = getVolume(x[0], x[1], x[2], x[3], veps[0], veps[1], veps[2], veps[3], veps[whichIndex], A, timeeps);
                              t[i][j] = (timeeps - time)/deps;
                            }
                          }
                          // end of finite difference for time gradient
                        }
                        else
                        {
                          PrismVolume::getTimeGradient(x[0], x[1], x[2], x[3], v[0], v[1], v[2], v[3], weights, time, m_thickness, t);
                        }
#ifdef DEBUGVG
                        // finite difference for time gradient
                        Vector3d teps[4];
                        for (size_t i=0; i<4; ++i)
                        {
                          for (size_t j=0; j<3; ++j)
                          { 
                            const double deps = 1e-8;
                            Vector3d veps[4];
                            for (size_t ii=0; ii<4; ++ii)
                            {
                              for (size_t jj=0; jj<3; ++jj)
                              {
                                veps[ii][jj] = v[ii][jj];
                              }
                            }
                            veps[i][j] = veps[i][j] + deps;
                            double timeeps = 1.0;
                            double V1eps = getVolume(x[0], x[1], x[2], x[3], veps[0], veps[1], veps[2], veps[3], veps[whichIndex], A, timeeps);
                            teps[i][j] = (timeeps - time)/deps;
                            t[i][j] = (timeeps - time)/deps;
                          }
                          cout<<"FV teps["<<i<<"] norm is : "<<teps[i].norm()<<endl;
                          cout<<"FV t["<<i<<"] norm is : "<<t[i].norm()<<endl;
                          cout<<"FV t["<<i<<"] difference norm is: "<<(teps[i] - t[i]).norm()<<endl;
                          for (size_t ii=0; ii<4; ++ii)
                            cout<<"FV weight of vertex: "<<ii<<" is: "<<weights[ii]<<endl;
                        }
                        // end of finite difference for time gradient
#endif
			Vector3d w[4][4];
			PrismVolume::getBarycentricGradients(r[0], r[1], r[2], r[3], w);
#ifdef DEBUGVG
                        // finite difference for w gradient
                        Vector3d weps[4][4];
			for (size_t i=0; i<4; ++i)
                        {
                                for (size_t j=0; j<3; ++j)
                                {
                                  double weightseps[4];
                                  const double deps = 1e-8;
                                  Vector3d reps[4];
                                  for (size_t ii=0; ii<4; ++ii)
                                  {
                                    for (size_t jj=0; jj<3; ++jj)
                                    {
                                      reps[ii][jj] = r[ii][jj];
                                    }
                                  }
                                  reps[i][j] += deps;
                                  PrismVolume::getBarycentricCoordinates(reps[0], reps[1], reps[2], reps[3], true, weightseps);
                                  
                                  for (size_t ii=0; ii<4; ++ii)
                                    weps[ii][i][j] = (weightseps[ii] - weights[ii]) / deps;
                                }
                        }
                        for (size_t i=0; i<4; ++i)
                        {
                          for (size_t j=0; j<4; ++j)
                          {
                            cout<<"FV w component: "<<i<<", "<<j<<". Difference norm: "<<(weps[i][j]-w[i][j]).norm()<<endl;
                          }
                        }
                        // end of finite difference for w gradient
#endif
			double weight = weights[whichIndex];
			Vector3d vGradient[4];
			for (size_t i=0; i<4; ++i)
			{
				Matrix3d N = Matrix3d::Zero();
				for (size_t j=0; j<3; ++j)
				{
					if (i == j)
						N += ((e[j].cross(n) / l) * n.transpose()) * (time * Matrix3d::Identity() + t[i] * v[j].transpose());
					else
						N += ((e[j].cross(n) / l) * n.transpose()) * (t[i] * v[j].transpose());
				}
#ifdef DEBUGVG
                                // finite difference for dn/dx
				Matrix3d Neps = Matrix3d::Zero();
				for (size_t j=0; j<3; ++j)
                                {
                                  const double deps = 1e-8;
                                  Vector3d veps[4];
                                  for (size_t ii=0; ii<4; ++ii)
                                  {
                                    for (size_t jj=0; jj<3; ++jj)
                                    {
                                      veps[ii][jj] = v[ii][jj];
                                    }
                                  }
                                  veps[i][j] += deps;
                                  
                                  double timeeps = 1.0;
			          double V1eps = getVolume(x[0], x[1], x[2], x[3], veps[0], veps[1], veps[2], veps[3], veps[whichIndex], A, timeeps);
                                  Vector3d reps[4];
                                  for (size_t ii=0; ii<4; ++ii)
                                          reps[ii] = x[ii] + timeeps * veps[ii];
                                  Vector3d normaleps = (reps[1] - reps[0]).cross(reps[2] - reps[0]);
                                  Vector3d neps(normaleps);
                                  neps.normalize();
                                  Vector3d dndx = (neps - n)/deps;
                                  for (size_t ii=0; ii<3; ++ii)
                                    Neps(ii,j) = dndx(ii);
                                }
                                cout<<"FV N: \n"<<N<<endl;
                                cout<<"FV Neps: \n"<<Neps<<endl;
                                cout<<"FV norm dndx difference: "<<(Neps-N).norm()<<endl;
                                // end of finite difference for dn/dx
#endif
				Vector3d wGrad = Vector3d::Zero();
                                /*
				for (size_t j=0; j<4; ++j)
				{
					if (i == j)
						wGrad += ((t[i] * v[j].transpose()) + time * Matrix3d::Identity()).transpose() * w[whichIndex][j];
					else
						wGrad += (t[i] * v[j].transpose()).transpose() * w[whichIndex][j];
				}
                                */
                                for (size_t j=0; j<4; ++j)
				{
					if (i == j)
						wGrad += ((v[j] * t[i].transpose()) + time * Matrix3d::Identity()).transpose() * w[whichIndex][j];
					else
						wGrad += (v[j] * t[i].transpose()).transpose() * w[whichIndex][j];
				}
#ifdef DEBUGVG
                                // finite difference for weight gradient wGrad
				Vector3d wGradeps = Vector3d::Zero();
				for (size_t j=0; j<3; ++j)
                                {
                                  double weightseps[4];
                                  const double deps = 1e-8;
                                  Vector3d veps[4];
                                  for (size_t ii=0; ii<4; ++ii)
                                  {
                                    for (size_t jj=0; jj<3; ++jj)
                                    {
                                      veps[ii][jj] = v[ii][jj];
                                    }
                                  }
                                  veps[i][j] += deps;
			          double timeeps = 1.0;
			          double V1eps = getVolume(x[0], x[1], x[2], x[3], veps[0], veps[1], veps[2], veps[3], veps[whichIndex], A, timeeps);
                                  Vector3d reps[4];
                                  for (size_t ii=0; ii<4; ++ii)
                                          reps[ii] = x[ii] + timeeps * veps[ii];
                                  PrismVolume::getBarycentricCoordinates(reps[0], reps[1], reps[2], reps[3], true, weightseps);
                                  wGradeps[j] = (weightseps[whichIndex]-weights[whichIndex])/deps;
                                  cout<<"FV testing wGrad: "<<whichIndex<<" -- "<<(weightseps[whichIndex]-weights[whichIndex])<<endl;
                                }
                                cout<<"FV wGradeps and wGrad difference norm: "<<(wGradeps-wGrad).norm()<<endl;
                                cout<<"FV wGradeps \n"<<wGradeps<<endl;
                                cout<<"FV wGrad \n"<<wGrad<<endl;
                                // end of finite difference for weight gradient wGrad
#endif
				vGradient[i] = wGrad * (1.0 - time) * v[whichIndex].dot(n) -
					weight * t[i] * v[whichIndex].dot(n) +
					weight * (1.0 - time) * N.transpose() * v[whichIndex];

				if (whichIndex == i)
					vGradient[i] += weight * (1.0 - time) * n;

				// TODO: Make the weights negative for triangle vertices. Should be the same, right?
				//
				//if (whichIndex == 0 || whichIndex == 1 || whichIndex == 2)
				//	vGradient[i] = -1.0 * vGradient[i];
                                // gradient direction is opposite to relative speed
                                // TODO what if vGradient.dot(rel_vel) < eps
                                if(i == 3)
                                {
                                  if(vGradient[i].dot(rel_vel) > 0)
                                    vGradient[i] = -vGradient[i];
                                }
                                else
                                {
                                  if(vGradient[i].dot(rel_vel) < 0)
                                    vGradient[i] = -vGradient[i];
                                }
			}

			for (size_t i=0; i<4; ++i)
				for (size_t j=0; j<3; ++j)
                                {
					gradients.coeffRef(3 * index[i] + j, 0) += vGradient[i][j] * A;
					//gradients.coeffRef(3 * index[i] + j, 0) += flip * vGradient[i][j] * A;
                                        //if(i==whichIndex)
					  //gradients.coeffRef(3 * index[i] + j, 0) += mesh.getVertex(index[i])->normal[j]*A;
#ifdef DEBUGVG
                                        // finite difference for vGradient 
                                        double weightseps[4];
                                        const double deps = 1e-8;
			                Vector3d veps[4];
                                        for (size_t ii=0; ii<4; ++ii)
                                        {
                                          for (size_t jj=0; jj<3; ++jj)
                                          {
                                            veps[ii][jj] = v[ii][jj];
                                          }
                                        }
                                        veps[i][j] += deps;

                                        double timeeps = 1.0;
                                        Vector3d reps[4];
                                        double V1eps = getVolume(x[0], x[1], x[2], x[3], veps[0], veps[1], veps[2], veps[3], veps[whichIndex], A, timeeps);
                                        for (size_t ii=0; ii<4; ++ii)
                                          reps[ii] = x[ii] + timeeps * veps[ii];
                                        PrismVolume::getBarycentricCoordinates(reps[0], reps[1], reps[2], reps[3], true, weightseps);
                                        V1eps *= weightseps[whichIndex];
                                        V1eps = flip*V1eps;
                                        double V1tmp = flip*V1;
                                        double grad_diff = (V1eps - V1tmp)/deps;
                                        /*
                                        cout<<"FV V1eps --- "<<V1eps<<endl;
                                        cout<<"FV V1tmp --- "<<V1tmp<<endl;
                                        cout<<"FV V1 ------ "<<V1<<endl;
                                        cout<<"FV volume gradient respect to vertex: "<<i<<" ,component: "<<j<<" is: "<<grad_diff<<" --- finite diff"<<endl; 
                                        cout<<"FV volume gradient respect to vertex: "<<i<<" ,component: "<<j<<" is: "<<(flip * vGradient[i][j] * A)<<" --- analytic"<<endl; 
                                        cout<<"FV diff of volume gradient respect to vertex: "<<i<<" ,component: "<<j<<" is: "<<grad_diff - (flip * vGradient[i][j] * A)<<endl; 
                                        for (size_t ii=0; ii<4; ++ii)
                                          cout<<"FV weight of vertex: "<<ii<<" is: "<<weights[ii]<<endl;
                                        */
                                        //gradients.coeffRef(3 * index[i] + j, 0) += grad_diff;
                                        // end of finite difference for vGradient
#endif
                                }
		}
		else
		{
			unsigned int index[4] = { ci->getFirstEdgeVertex1(), ci->getFirstEdgeVertex2(),
									  ci->getSecondEdgeVertex1(), ci->getSecondEdgeVertex2() };

			Vector3d x[4], v[4];
			for (size_t i=0; i<4; ++i)
			{
				x[i] = mesh.getVertex(index[i])->Pori();
				v[i] = mesh.getVertex(index[i])->P() - mesh.getVertex(index[i])->Pori();
			}
                        
                        if(m_periodicLength > 0)
                        {
                            Vector3d transV(0,0,0);
                            PrismVolume::getTranslationVector(transV, 
                                    mesh.getVertex(index[0]), mesh.getVertex(index[2]), m_periodicLength);
                            x[0] = x[0] + transV;
                            x[1] = x[1] + transV;
                        }

			unsigned int whichIndex = 3;
			for (size_t i=0; i<3; ++i)
			{
				if (cmItr->first == index[i])
				{
					whichIndex = i;
					break;
				}
			}

			//double A = getVoronoiArea(mesh, index[whichIndex]);
			double A = 1.0;
			double time = 1.0;
			double V1 = getEEVolume(x[0], x[1], x[2], x[3], v[0], v[1], v[2], v[3], v[whichIndex], A, time);

			// Compute the positions at the moment of intersection
			//
			Vector3d r[4];
			for (size_t i=0; i<4; ++i)
				r[i] = x[i] + time * v[i];

			// TODO
			double weights[4];
			PrismVolume::getBarycentricCoordinates(r[0], r[1], r[2], r[3], false, weights);
                        assert(weights[0]>=0 && weights[0] <= 1);
                        assert(weights[1]>=0 && weights[1] <= 1);
                        assert(weights[2]>=0 && weights[2] <= 1);
                        assert(weights[3]>=0 && weights[3] <= 1);
			//weights[0] = (1.0 - weights[1]); weights[2] = (1.0 - weights[3]);

                        //std::cout<<"this is vf contact? "<<ci->isVF()<<". --- multipling weight: "<<weights[whichIndex]<<std::endl; 
			V1 *= weights[whichIndex];

			Vector3d normal = (r[1] - r[0]).cross(r[3] - r[2]);
			Vector3d n(normal);
			double l = n.norm();
			n.normalize();

                        //assert(V1 != 0.0);
			double flip = 1.0;
			if (V1 > 0.0)
			{
				m_volume -= V1;
                                flip = -1.0;
				//if (whichIndex == 2 || whichIndex == 3)
				//	flip = -1.0;
			}
			else if (V1 < 0.0)
			{
				m_volume += V1;
				//if (whichIndex == 0 || whichIndex == 1)
				//	flip = -1.0;
			}
			else // Volume is 0.0, insufficient information to determine flip
			{
				// Use relative velocity in normal direction to determine "correct" normal orientation
				//
				//double relVel = (((1.0 - weights[0]) * v[0] + weights[1] * v[1]) - (1.0 - weights[2]) * v[2] + weights[3] * v[3]).dot(n);
				double relVel = ( (weights[0] * v[0] + weights[1] * v[1]) - (weights[2] * v[2] + weights[3] * v[3]) ).dot(n);
				if (relVel > 0.0)
					flip = -1.0;
			}
                        Vector3d rel_vel = ( (weights[0]*v[0] + weights[1]*v[1]) - (weights[2]*v[2] + weights[3]*v[3]) );

			// The four edge vectors
			//
			Vector3d e[4];
                        e[0] = r[3] - r[2];
			e[1] = r[2] - r[3];
                        e[2] = r[0] - r[1];
			e[3] = r[1] - r[0];

			Vector3d t[4];
                        if(1)
                        //if(weights[0]<0 || weights[1]<0 || weights[2]<0 || weights[3]<0)
                        {
                          // finite difference for time gradient
                          for (size_t i=0; i<4; ++i)
                          {
                            for (size_t j=0; j<3; ++j)
                            { 
                              const double deps = 1e-8;
                              Vector3d veps[4];
                              for (size_t ii=0; ii<4; ++ii)
                              {
                                for (size_t jj=0; jj<3; ++jj)
                                {
                                  veps[ii][jj] = v[ii][jj];
                                }
                              }
                              veps[i][j] = veps[i][j] + deps;
                              double timeeps = 1.0;
                              double V1eps = getEEVolume(x[0], x[1], x[2], x[3], veps[0], veps[1], veps[2], veps[3], veps[whichIndex], A, timeeps);
                              t[i][j] = (timeeps - time)/deps;
                            }
                          }
                          // end of finite difference for time gradient
                          /*
                          cout<<"weight[0]: "<<weights[0]<<endl;
                          cout<<"r[0]: "<<r[0]<<endl;
                          cout<<"weight[1]: "<<weights[1]<<endl;
                          cout<<"r[1]: "<<r[1]<<endl;
                          cout<<"weight[2]: "<<weights[2]<<endl;
                          cout<<"r[2]: "<<r[2]<<endl;
                          cout<<"weight[3]: "<<weights[3]<<endl;
                          cout<<"r[3]: "<<r[3]<<endl;
                          */
                        }
                        else
                        {
                          PrismVolume::getEETimeGradient(x[0], x[1], x[2], x[3], v[0], v[1], v[2], v[3], weights, time, m_thickness, t);
                        }
#ifdef DEBUGVG
                        // finite difference for time gradient
                        Vector3d teps[4];
                        for (size_t i=0; i<4; ++i)
                        {
                          for (size_t j=0; j<3; ++j)
                          { 
                            const double deps = 1e-8;
                            Vector3d veps[4];
                            for (size_t ii=0; ii<4; ++ii)
                            {
                              for (size_t jj=0; jj<3; ++jj)
                              {
                                veps[ii][jj] = v[ii][jj];
                              }
                            }
                            veps[i][j] = veps[i][j] + deps;
                            double timeeps = 1.0;
                            double V1eps = getEEVolume(x[0], x[1], x[2], x[3], veps[0], veps[1], veps[2], veps[3], veps[whichIndex], A, timeeps);
                            teps[i][j] = (timeeps - time)/deps;
                          }
                          cout<<"EE teps["<<i<<"] norm is : "<<teps[i].norm()<<endl;
                          cout<<"EE t["<<i<<"] norm is : "<<t[i].norm()<<endl;
                          cout<<"EE t["<<i<<"] difference norm is: "<<(teps[i] - t[i]).norm()<<endl;
                          for (size_t ii=0; ii<4; ++ii)
                            cout<<"EE weight of vertex: "<<ii<<" is: "<<weights[ii]<<endl;
                        }
                        // end of finite difference for time gradient
#endif                        
			Vector3d w[4][4];
			PrismVolume::getEEBarycentricGradients(r[0], r[1], r[2], r[3], w);
#ifdef DEBUGVG
                        // finite difference for w gradient
                        Vector3d weps[4][4];
			for (size_t i=0; i<4; ++i)
                        {
                                for (size_t j=0; j<3; ++j)
                                {
                                  double weightseps[4];
                                  const double deps = 1e-8;
                                  Vector3d reps[4];
                                  for (size_t ii=0; ii<4; ++ii)
                                  {
                                    for (size_t jj=0; jj<3; ++jj)
                                    {
                                      reps[ii][jj] = r[ii][jj];
                                    }
                                  }
                                  reps[i][j] += deps;
                                  PrismVolume::getBarycentricCoordinates(reps[0], reps[1], reps[2], reps[3], false, weightseps);
                                  
                                  for (size_t ii=0; ii<4; ++ii)
                                    weps[ii][i][j] = (weightseps[ii] - weights[ii]) / deps;
                                }
                        }
                        for (size_t i=0; i<4; ++i)
                        {
                          for (size_t j=0; j<4; ++j)
                          {
                            cout<<"EE w component: "<<i<<", "<<j<<". Difference norm: "<<(weps[i][j]-w[i][j]).norm()<<endl;
                          }
                        }
                        // end of finite difference for w gradient
#endif
			double weight = weights[whichIndex];
			Vector3d vGradient[4];
			for (size_t i=0; i<4; ++i)
			{
				Matrix3d N = Matrix3d::Zero();
				for (size_t j=0; j<4; ++j)
				{
					if (i == j)
						N += ((e[j].cross(n) / l) * n.transpose()) * (time * Matrix3d::Identity() + t[i] * v[j].transpose());
					else
						N += ((e[j].cross(n) / l) * n.transpose()) * (t[i] * v[j].transpose());
				}
#ifdef DEBUGVG
                                // finite difference for dn/dx
				Matrix3d Neps = Matrix3d::Zero();
				for (size_t j=0; j<3; ++j)
                                {
                                  const double deps = 1e-8;
                                  Vector3d veps[4];
                                  for (size_t ii=0; ii<4; ++ii)
                                  {
                                    for (size_t jj=0; jj<3; ++jj)
                                    {
                                      veps[ii][jj] = v[ii][jj];
                                    }
                                  }
                                  veps[i][j] += deps;
                                  
                                  double timeeps = 1.0;
			          double V1eps = getEEVolume(x[0], x[1], x[2], x[3], veps[0], veps[1], veps[2], veps[3], veps[whichIndex], A, timeeps);
                                  Vector3d reps[4];
                                  for (size_t ii=0; ii<4; ++ii)
                                          reps[ii] = x[ii] + timeeps * veps[ii];
                                  Vector3d normaleps = (reps[1] - reps[0]).cross(reps[3] - reps[2]);
                                  Vector3d neps(normaleps);
                                  neps.normalize();
                                  Vector3d dndx = (neps - n)/deps;
                                  for (size_t ii=0; ii<3; ++ii)
                                    Neps(ii,j) = dndx(ii);
                                }
                                cout<<"EE N: \n"<<N<<endl;
                                cout<<"EE Neps: \n"<<Neps<<endl;
                                cout<<"EE norm dndx difference: "<<(Neps-N).norm()<<endl;
                                // end of finite difference for dn/dx
#endif
				Vector3d wGrad = Vector3d::Zero();
                                /*
				for (size_t j=0; j<4; ++j)
				{
					if (i == j)
						wGrad += ((t[i] * v[j].transpose()) + time * Matrix3d::Identity()).transpose() * w[whichIndex][j];
					else
						wGrad += (t[i] * v[j].transpose()).transpose() * w[whichIndex][j];
				}
                                */
                                for (size_t j=0; j<4; ++j)
				{
					if (i == j)
						wGrad += ((v[j] * t[i].transpose()) + time * Matrix3d::Identity()).transpose() * w[whichIndex][j];
					else
						wGrad += (v[j] * t[i].transpose()).transpose() * w[whichIndex][j];
				}
#ifdef DEBUGVG
                                // finite difference for weight gradient wGrad
				Vector3d wGradeps = Vector3d::Zero();
				for (size_t j=0; j<3; ++j)
                                {
                                  double weightseps[4];
                                  const double deps = 1e-8;
                                  Vector3d veps[4];
                                  for (size_t ii=0; ii<4; ++ii)
                                  {
                                    for (size_t jj=0; jj<3; ++jj)
                                    {
                                      veps[ii][jj] = v[ii][jj];
                                    }
                                  }
                                  veps[i][j] += deps;
			          double timeeps = 1.0;
			          double V1eps = getEEVolume(x[0], x[1], x[2], x[3], veps[0], veps[1], veps[2], veps[3], veps[whichIndex], A, timeeps);
                                  Vector3d reps[4];
                                  for (size_t ii=0; ii<4; ++ii)
                                          reps[ii] = x[ii] + timeeps * veps[ii];
                                  PrismVolume::getBarycentricCoordinates(reps[0], reps[1], reps[2], reps[3], false, weightseps);
                                  wGradeps[j] = (weightseps[whichIndex]-weights[whichIndex])/deps;
                                  cout<<"EE testing wGrad: "<<whichIndex<<" -- "<<(weightseps[whichIndex]-weights[whichIndex])<<endl;
                                }
                                cout<<"EE wGradeps and wGrad difference norm: "<<(wGradeps-wGrad).norm()<<endl;
                                cout<<"EE wGradeps \n"<<wGradeps<<endl;
                                cout<<"EE wGrad \n"<<wGrad<<endl;
                                // end of finite difference for weight gradient wGrad
#endif
				vGradient[i] = wGrad * (1.0 - time) * v[whichIndex].dot(n) -
							   weight * t[i] * v[whichIndex].dot(n) +
							   weight * (1.0 - time) * N.transpose() * v[whichIndex];
				//vGradient[i] = -weight * t[i] * v[whichIndex].dot(n) +
			        //			  weight * (1.0 - time) * N.transpose() * v[whichIndex];

				if (whichIndex == i)
					vGradient[i] += weight * (1.0 - time) * n;

				//if (whichIndex == 2 || whichIndex == 3)
				//if (whichIndex == 0 || whichIndex == 1)
				//	vGradient[i] = -1.0 * vGradient[i];
                                // gradient direction is opposite to relative speed
                                // TODO what if vGradient.dot(rel_vel) < eps
                                if(i == 0 || i == 1)
                                {
                                  if(vGradient[i].dot(rel_vel) > 0)
                                    vGradient[i] = -vGradient[i];
                                }
                                else
                                {
                                  if(vGradient[i].dot(rel_vel) < 0)
                                    vGradient[i] = -vGradient[i];
                                }
			}

			for (size_t i=0; i<4; ++i)
				for (size_t j=0; j<3; ++j)
                                {
				        gradients.coeffRef(3 * index[i] + j, 0) += vGradient[i][j] * A;
				        //gradients.coeffRef(3 * index[i] + j, 0) += flip * vGradient[i][j] * A;
                                        //if(i==whichIndex)
					  //gradients.coeffRef(3 * index[i] + j, 0) += mesh.getVertex(index[i])->normal[j]*A;
#ifdef DEBUGVG
                                        // finite difference for vGradient 
                                        double weightseps[4];
                                        const double deps = 1e-8;
			                Vector3d veps[4];
                                        for (size_t ii=0; ii<4; ++ii)
                                        {
                                          for (size_t jj=0; jj<3; ++jj)
                                          {
                                            veps[ii][jj] = v[ii][jj];
                                          }
                                        }
                                        veps[i][j] += deps;

                                        double timeeps = 1.0;
                                        Vector3d reps[4];
                                        double V1eps = getEEVolume(x[0], x[1], x[2], x[3], veps[0], veps[1], veps[2], veps[3], veps[whichIndex], A, timeeps);
                                        for (size_t ii=0; ii<4; ++ii)
                                          reps[ii] = x[ii] + timeeps * veps[ii];
                                        PrismVolume::getBarycentricCoordinates(reps[0], reps[1], reps[2], reps[3], false, weightseps);
                                        V1eps *= weightseps[whichIndex];
                                        V1eps = flip*V1eps;
                                        double V1tmp = flip*V1;
                                        double grad_diff = (V1eps - V1tmp)/deps;
                                        /*
                                        cout<<"EE V1eps --- "<<V1eps<<endl;
                                        cout<<"EE V1tmp --- "<<V1tmp<<endl;
                                        cout<<"EE V1 ------ "<<V1<<endl;
                                        cout<<"EE weighteps ------ "<<weightseps[whichIndex]<<endl;
                                        cout<<"EE weight ------ "<<weights[whichIndex]<<endl;
                                        cout<<"EE volume gradient respect to vertex: "<<i<<" ,component: "<<j<<" is: "<<grad_diff<<" --- finite diff"<<endl; 
                                        cout<<"EE volume gradient respect to vertex: "<<i<<" ,component: "<<j<<" is: "<<(flip * vGradient[i][j] * A)<<" --- analytic"<<endl; 
                                        cout<<"EE collision time respect to vertex: "<<i<<" ,component: "<<j<<" is: "<<timeeps<<" --- finite diff"<<endl; 
                                        cout<<"EE collision time respect to vertex: "<<i<<" ,component: "<<j<<" is: "<<time<<" --- analytic"<<endl; 
                                        cout<<"EE diff of volume gradient respect to vertex: "<<i<<" ,component: "<<j<<" is: "<<grad_diff - (flip * vGradient[i][j] * A)<<endl; 
                                        for (size_t ii=0; ii<4; ++ii)
                                          cout<<"EE weight of vertex: "<<ii<<" is: "<<weights[ii]<<endl;
                                        */
                                        //gradients.coeffRef(3 * index[i] + j, 0) += grad_diff;
                                        // end of finite difference for vGradient
#endif
                                }
		}
	}

	m_volume *= 1.00;
        m_volume = min(m_volume,-1e-10);

	gradients.finalize();
	m_gradients = SparseVector<double>(gradients);
}

void IntersectionVolume::computeVolumeNew(Polyhedron &mesh, CollisionMap &collisions)
{
	DynamicSparseMatrix<double> gradients(3 * mesh.size_of_vertices(), 1);

	m_volume = 0.0;
        double time_constant = 0.0002;
	for (CollisionMapIterator cmItr=collisions.begin(); cmItr!=collisions.end(); ++cmItr)
	{
		CollisionsIterator ci = cmItr->second;

		if (ci->isVF())
		{
			unsigned int index[4] = { ci->getTriangleVertex1(), ci->getTriangleVertex2(), ci->getTriangleVertex3(), ci->getVertex() };

			Vector3d x[4], v[4];
			for (size_t i=0; i<4; ++i)
			{
				x[i] = mesh.getVertex(index[i])->Pori();
				v[i] = mesh.getVertex(index[i])->P() - mesh.getVertex(index[i])->Pori();
				if(mesh.getVertex(index[i])->meshTag>=m_nv)
                                    m_bdv = 1;
			}

                        if(m_periodicLength > 0)
                        {
                            Vector3d transV(0,0,0);
                            PrismVolume::getTranslationVector(transV, 
                                    mesh.getVertex(index[3]), mesh.getVertex(index[0]), m_periodicLength);
                            x[3] = x[3] + transV;
                        }

			unsigned int whichIndex = 3;
			for (size_t i=0; i<3; ++i)
			{
				if (cmItr->first == index[i])
				{
					whichIndex = i;
					break;
				}
			}

			//double A = getVoronoiArea(mesh, index[whichIndex]);
			double A = 1.0;
			double time = 1.0;
                        double v_partial = 0.0;
			double V1 = getVolume(x[0], x[1], x[2], x[3], v[0], v[1], v[2], v[3], v[whichIndex], A, time, time_constant, v_partial);
                        double v_partial_time = sqrt(time_constant*time_constant + v_partial*v_partial);

			// Compute the positions at the moment of intersection
			Vector3d r[4];
			for (size_t i=0; i<4; ++i)
				r[i] = x[i] + time * v[i];

			double weights[4];
			PrismVolume::getBarycentricCoordinates(r[0], r[1], r[2], r[3], true, weights);
                        assert(weights[0]>=0 && weights[0] <= 1);
                        assert(weights[1]>=0 && weights[1] <= 1);
                        assert(weights[2]>=0 && weights[2] <= 1);
                        assert(weights[3]>=0 && weights[3] <= 1);

                        //std::cout<<"this is vf contact? "<<ci->isVF()<<". --- multipling weight: "<<weights[whichIndex]<<std::endl; 
			//V1 *= weights[whichIndex];

			Vector3d normal = (r[1] - r[0]).cross(r[2] - r[0]);
			Vector3d n(normal);
			double l = n.norm();
			n.normalize();


			// NOTE: flip should only be -1.0 if hitting the backside of a face, which should
			// only happen on a surface with boundary, maybe assert this?
			//
                        //std::cout<<"V1: "<<V1<<". time: "<<time<<std::endl;
                        //assert(V1 != 0.0);
			double flip = 1.0;
			if (V1 > 0.0)
			{
				m_volume -= V1;
                                flip = -1.0;
				//if (whichIndex == 3)
				//	flip = -1.0;
			}
			else if (V1 < 0.0)
			{
				m_volume += V1;
				//if (whichIndex != 3)
				//	flip = -1.0;
			}
			else // Volume is 0.0, insufficient information to determine flip
			{
				// Use relative velocity in normal direction to determine "correct" normal orientation,
				// since gradients may still be non-zero, and we need a consistent direction
				//
				double relVel = (v[3] - (weights[0] * v[0] + weights[1] * v[1] + weights[2] * v[2])).dot(n);
				if (relVel > 0.0)
					flip = -1.0;
			}
                        Vector3d rel_vel = v[3] - (weights[0] * v[0] + weights[1] * v[1] + weights[2] * v[2]);

			// The three edge vectors
			//
			Vector3d e[3];
			e[0] = r[2] - r[1];
			e[1] = r[0] - r[2];
			e[2] = r[1] - r[0];

			Vector3d t[4];
                        // TODO: fix analytic time gradient when vertex face contact is 
                        // degenerated to vertex-vertex or vertex-edge contact.
                        if(1)
                        //if(weights[0]<0 || weights[1]<0 || weights[2]<0 || weights[3]<0)
                        //if(weights[0]<=0 || weights[1]<=0 || weights[2]<=0)
                        {
                          // finite difference for time gradient
                          for (size_t i=0; i<4; ++i)
                          {
                            for (size_t j=0; j<3; ++j)
                            { 
                              //const double deps = 1e-8;
                              const double deps = 1e-12 + fabs(v[i][j]*1e-8);
                              Vector3d veps[4];
                              for (size_t ii=0; ii<4; ++ii)
                              {
                                for (size_t jj=0; jj<3; ++jj)
                                {
                                  veps[ii][jj] = v[ii][jj];
                                }
                              }
                              veps[i][j] = veps[i][j] + deps;
                              double timeeps = 1.0;
                              double v_partialeps = 0.0;
                              double V1eps = getVolume(x[0], x[1], x[2], x[3], veps[0], veps[1], veps[2], veps[3], veps[whichIndex], A, timeeps, time_constant, v_partialeps);
                              t[i][j] = (timeeps - time)/deps;
                            }
                          }
                          // end of finite difference for time gradient
                        }
                        else
                        {
                          PrismVolume::getTimeGradient(x[0], x[1], x[2], x[3], v[0], v[1], v[2], v[3], weights, time, m_thickness, t);
                        }
#ifdef DEBUGVG
                        // finite difference for time gradient
                        Vector3d teps[4];
                        for (size_t i=0; i<4; ++i)
                        {
                          for (size_t j=0; j<3; ++j)
                          { 
                            const double deps = 1e-8;
                            Vector3d veps[4];
                            for (size_t ii=0; ii<4; ++ii)
                            {
                              for (size_t jj=0; jj<3; ++jj)
                              {
                                veps[ii][jj] = v[ii][jj];
                              }
                            }
                            veps[i][j] = veps[i][j] + deps;
                            double timeeps = 1.0;
                            double v_partialeps = 0.0;
                            double V1eps = getVolume(x[0], x[1], x[2], x[3], veps[0], veps[1], veps[2], veps[3], veps[whichIndex], A, timeeps, time_constant, v_partialeps);
                            teps[i][j] = (timeeps - time)/deps;
                          }
                          cout<<"FV teps["<<i<<"] norm is : "<<teps[i].norm()<<endl;
                          cout<<"FV t["<<i<<"] norm is : "<<t[i].norm()<<endl;
                          cout<<"FV t["<<i<<"] difference norm is: "<<(teps[i] - t[i]).norm()<<endl;
                          for (size_t ii=0; ii<4; ++ii)
                            cout<<"FV weight of vertex: "<<ii<<" is: "<<weights[ii]<<endl;
                        }
                        // end of finite difference for time gradient
#endif
			Vector3d w[4][4];
			PrismVolume::getBarycentricGradients(r[0], r[1], r[2], r[3], w);
#ifdef DEBUGVG
                        // finite difference for w gradient
                        Vector3d weps[4][4];
			for (size_t i=0; i<4; ++i)
                        {
                                for (size_t j=0; j<3; ++j)
                                {
                                  double weightseps[4];
                                  const double deps = 1e-8;
                                  Vector3d reps[4];
                                  for (size_t ii=0; ii<4; ++ii)
                                  {
                                    for (size_t jj=0; jj<3; ++jj)
                                    {
                                      reps[ii][jj] = r[ii][jj];
                                    }
                                  }
                                  reps[i][j] += deps;
                                  PrismVolume::getBarycentricCoordinates(reps[0], reps[1], reps[2], reps[3], true, weightseps);
                                  
                                  for (size_t ii=0; ii<4; ++ii)
                                    weps[ii][i][j] = (weightseps[ii] - weights[ii]) / deps;
                                }
                        }
                        for (size_t i=0; i<4; ++i)
                        {
                          for (size_t j=0; j<4; ++j)
                          {
                            cout<<"FV w component: "<<i<<", "<<j<<". Difference norm: "<<(weps[i][j]-w[i][j]).norm()<<endl;
                          }
                        }
                        // end of finite difference for w gradient
#endif
			double weight = weights[whichIndex];
			Vector3d vGradient[4];
			for (size_t i=0; i<4; ++i)
			{
				Matrix3d N = Matrix3d::Zero();
				for (size_t j=0; j<3; ++j)
				{
					if (i == j)
						N += ((e[j].cross(n) / l) * n.transpose()) * (time * Matrix3d::Identity() + t[i] * v[j].transpose());
					else
						N += ((e[j].cross(n) / l) * n.transpose()) * (t[i] * v[j].transpose());
				}
#ifdef DEBUGVG
                                // finite difference for dn/dx
				Matrix3d Neps = Matrix3d::Zero();
				for (size_t j=0; j<3; ++j)
                                {
                                  const double deps = 1e-8;
                                  Vector3d veps[4];
                                  for (size_t ii=0; ii<4; ++ii)
                                  {
                                    for (size_t jj=0; jj<3; ++jj)
                                    {
                                      veps[ii][jj] = v[ii][jj];
                                    }
                                  }
                                  veps[i][j] += deps;
                                  
                                  double timeeps = 1.0;
                                  double v_partialeps = 0.0;
			          double V1eps = getVolume(x[0], x[1], x[2], x[3], veps[0], veps[1], veps[2], veps[3], veps[whichIndex], A, timeeps, time_constant, v_partialeps);
                                  Vector3d reps[4];
                                  for (size_t ii=0; ii<4; ++ii)
                                          reps[ii] = x[ii] + timeeps * veps[ii];
                                  Vector3d normaleps = (reps[1] - reps[0]).cross(reps[2] - reps[0]);
                                  Vector3d neps(normaleps);
                                  neps.normalize();
                                  Vector3d dndx = (neps - n)/deps;
                                  for (size_t ii=0; ii<3; ++ii)
                                    Neps(ii,j) = dndx(ii);
                                }
                                cout<<"FV N: \n"<<N<<endl;
                                cout<<"FV Neps: \n"<<Neps<<endl;
                                cout<<"FV norm dndx difference: "<<(Neps-N).norm()<<endl;
                                // end of finite difference for dn/dx
#endif
				Vector3d wGrad = Vector3d::Zero();
                                /*
				for (size_t j=0; j<4; ++j)
				{
					if (i == j)
						wGrad += ((t[i] * v[j].transpose()) + time * Matrix3d::Identity()).transpose() * w[whichIndex][j];
					else
						wGrad += (t[i] * v[j].transpose()).transpose() * w[whichIndex][j];
				}
                                */
                                for (size_t j=0; j<4; ++j)
				{
					if (i == j)
						wGrad += ((v[j] * t[i].transpose()) + time * Matrix3d::Identity()).transpose() * w[whichIndex][j];
					else
						wGrad += (v[j] * t[i].transpose()).transpose() * w[whichIndex][j];
				}
#ifdef DEBUGVG
                                // finite difference for weight gradient wGrad
				Vector3d wGradeps = Vector3d::Zero();
				for (size_t j=0; j<3; ++j)
                                {
                                  double weightseps[4];
                                  const double deps = 1e-8;
                                  Vector3d veps[4];
                                  for (size_t ii=0; ii<4; ++ii)
                                  {
                                    for (size_t jj=0; jj<3; ++jj)
                                    {
                                      veps[ii][jj] = v[ii][jj];
                                    }
                                  }
                                  veps[i][j] += deps;
			          double timeeps = 1.0;
                                  double v_partialeps;
			          double V1eps = getVolume(x[0], x[1], x[2], x[3], veps[0], veps[1], veps[2], veps[3], veps[whichIndex], A, timeeps, time_constant, v_partialeps);
                                  Vector3d reps[4];
                                  for (size_t ii=0; ii<4; ++ii)
                                          reps[ii] = x[ii] + timeeps * veps[ii];
                                  PrismVolume::getBarycentricCoordinates(reps[0], reps[1], reps[2], reps[3], true, weightseps);
                                  wGradeps[j] = (weightseps[whichIndex]-weights[whichIndex])/deps;
                                  cout<<"FV testing wGrad: "<<whichIndex<<" -- "<<(weightseps[whichIndex]-weights[whichIndex])<<endl;
                                }
                                cout<<"FV wGradeps and wGrad difference norm: "<<(wGradeps-wGrad).norm()<<endl;
                                cout<<"FV wGradeps \n"<<wGradeps<<endl;
                                cout<<"FV wGrad \n"<<wGrad<<endl;
                                // end of finite difference for weight gradient wGrad
#endif
                                
                                vGradient[i] = t[i]*v_partial_time + v_partial/v_partial_time*(time-1.0)*N.transpose()*v[whichIndex];

				if (whichIndex == i)
                                {
					//vGradient[i] += v_partial/v_partial_time*(time-1.0)*n;
                                        Vector3d normal_tmp;
                                        normal_tmp[0] = mesh.getVertex(index[i])->normal[0];
                                        normal_tmp[1] = mesh.getVertex(index[i])->normal[1];
                                        normal_tmp[2] = mesh.getVertex(index[i])->normal[2];
                                        if(normal_tmp.dot(n)<=0)
                                            vGradient[i] -= fabs(v_partial/v_partial_time*(time-1.0))*n;
                                        else
                                            vGradient[i] += fabs(v_partial/v_partial_time*(time-1.0))*n;
                                }
                                
                                /*
                                if( (x[0][0]>-4.93) && (x[0][0]<-3.93) && (x[0][1]>0.9) && (x[0][1]<1.9) && (x[0][2]>-2.33) && (x[0][2]<-1.33) )
                                {
                                    cout<<"index: "<<i<<endl;
                                    cout<<"global index: "<<index[i]<<endl;
                                    cout<<"which index: "<<whichIndex<<endl;
                                    cout<<"collision time: "<<time<<endl;
                                    cout<<"weights: "<<weights[0]<<", "<<weights[1]<<", "<<weights[2]<<", "<<weights[3]<<endl;
                                    cout<<"ti: "<<t[i][0]<<", "<<t[i][1]<<", "<<t[i][2]<<endl;
                                    cout<<"vi: "<<v[i][0]<<", "<<v[i][1]<<", "<<v[i][2]<<endl;
                                    cout<<"v_partial: "<<v_partial<<endl;
                                    cout<<"v_partial_time: "<<v_partial_time<<endl;
                                }
                                */

                                /*
                                Vector3d normal_tmp;
                                normal_tmp[0] = mesh.getVertex(index[i])->normal[0];
                                normal_tmp[1] = mesh.getVertex(index[i])->normal[1];
                                normal_tmp[2] = mesh.getVertex(index[i])->normal[2];
                                //if(whichIndex == i)
                                    //std::cout<<"normal dot grad: "<<normal_tmp.dot(vGradient[i])<<std::endl;
                                if( (normal_tmp.dot(vGradient[i])<=0) && (whichIndex==i) )
                                    //vGradient[i] = -vGradient[i];
                                    vGradient[i] = normal_tmp;
                                */
                                // gradient direction is opposite to relative speed
                                //if(i == 3)
                                //{
                                //  if(vGradient[i].dot(rel_vel) > 0)
                                //    vGradient[i] = -vGradient[i];
                                //}
                                //else
                                //{
                                //  if(vGradient[i].dot(rel_vel) < 0)
                                //    vGradient[i] = -vGradient[i];
                                //}
			}

			for (size_t i=0; i<4; ++i)
				for (size_t j=0; j<3; ++j)
                                {
					gradients.coeffRef(3 * index[i] + j, 0) += vGradient[i][j];
#ifdef DEBUGVG
                                        // finite difference for vGradient 
                                        double weightseps[4];
                                        const double deps = 1e-8;
			                Vector3d veps[4];
                                        for (size_t ii=0; ii<4; ++ii)
                                        {
                                          for (size_t jj=0; jj<3; ++jj)
                                          {
                                            veps[ii][jj] = v[ii][jj];
                                          }
                                        }
                                        veps[i][j] += deps;

                                        double timeeps = 1.0;
                                        double v_partialeps = 0.0;
                                        Vector3d reps[4];
                                        double V1eps = getVolume(x[0], x[1], x[2], x[3], veps[0], veps[1], veps[2], veps[3], veps[whichIndex], A, timeeps, time_constant, v_partialeps);
                                        for (size_t ii=0; ii<4; ++ii)
                                          reps[ii] = x[ii] + timeeps * veps[ii];
                                        PrismVolume::getBarycentricCoordinates(reps[0], reps[1], reps[2], reps[3], true, weightseps);
                                        V1eps *= weightseps[whichIndex];
                                        V1eps = flip*V1eps;
                                        double V1tmp = flip*V1;
                                        double grad_diff = (V1eps - V1tmp)/deps;
                                        /*
                                        cout<<"FV V1eps --- "<<V1eps<<endl;
                                        cout<<"FV V1tmp --- "<<V1tmp<<endl;
                                        cout<<"FV V1 ------ "<<V1<<endl;
                                        cout<<"FV volume gradient respect to vertex: "<<i<<" ,component: "<<j<<" is: "<<grad_diff<<" --- finite diff"<<endl; 
                                        cout<<"FV volume gradient respect to vertex: "<<i<<" ,component: "<<j<<" is: "<<(flip * vGradient[i][j] * A)<<" --- analytic"<<endl; 
                                        cout<<"FV diff of volume gradient respect to vertex: "<<i<<" ,component: "<<j<<" is: "<<grad_diff - (flip * vGradient[i][j] * A)<<endl; 
                                        for (size_t ii=0; ii<4; ++ii)
                                          cout<<"FV weight of vertex: "<<ii<<" is: "<<weights[ii]<<endl;
                                        */
                                        //gradients.coeffRef(3 * index[i] + j, 0) += grad_diff;
                                        // end of finite difference for vGradient
#endif
                                }
		}
		else
		{
			unsigned int index[4] = { ci->getFirstEdgeVertex1(), ci->getFirstEdgeVertex2(),
									  ci->getSecondEdgeVertex1(), ci->getSecondEdgeVertex2() };

			Vector3d x[4], v[4];
			for (size_t i=0; i<4; ++i)
			{
				x[i] = mesh.getVertex(index[i])->Pori();
				v[i] = mesh.getVertex(index[i])->P() - mesh.getVertex(index[i])->Pori();
                                if(mesh.getVertex(index[i])->meshTag>=m_nv)
                                    m_bdv = 1;
			}
                        
                        if(m_periodicLength > 0)
                        {
                            Vector3d transV(0,0,0);
                            PrismVolume::getTranslationVector(transV, 
                                    mesh.getVertex(index[0]), mesh.getVertex(index[2]), m_periodicLength);
                            x[0] = x[0] + transV;
                            x[1] = x[1] + transV;
                        }

			unsigned int whichIndex = 3;
			for (size_t i=0; i<3; ++i)
			{
				if (cmItr->first == index[i])
				{
					whichIndex = i;
					break;
				}
			}

			//double A = getVoronoiArea(mesh, index[whichIndex]);
			double A = 1.0;
			double time = 1.0;
                        double v_partial = 0.0;
			double V1 = getEEVolume(x[0], x[1], x[2], x[3], v[0], v[1], v[2], v[3], v[whichIndex], A, time, time_constant, v_partial);
                        double v_partial_time = sqrt(time_constant*time_constant + v_partial*v_partial);

			// Compute the positions at the moment of intersection
			//
			Vector3d r[4];
			for (size_t i=0; i<4; ++i)
				r[i] = x[i] + time * v[i];

			// TODO
			double weights[4];
			PrismVolume::getBarycentricCoordinates(r[0], r[1], r[2], r[3], false, weights);
                        assert(weights[0]>=0 && weights[0] <= 1);
                        assert(weights[1]>=0 && weights[1] <= 1);
                        assert(weights[2]>=0 && weights[2] <= 1);
                        assert(weights[3]>=0 && weights[3] <= 1);

                        //std::cout<<"this is vf contact? "<<ci->isVF()<<". --- multipling weight: "<<weights[whichIndex]<<std::endl; 
			//V1 *= weights[whichIndex];

			Vector3d normal = (r[1] - r[0]).cross(r[3] - r[2]);
			Vector3d n(normal);
			double l = n.norm();
			n.normalize();

                        //assert(V1 != 0.0);
			double flip = 1.0;
			if (V1 > 0.0)
			{
				m_volume -= V1;
                                flip = -1.0;
				//if (whichIndex == 2 || whichIndex == 3)
				//	flip = -1.0;
			}
			else if (V1 < 0.0)
			{
				m_volume += V1;
				//if (whichIndex == 0 || whichIndex == 1)
				//	flip = -1.0;
			}
			else // Volume is 0.0, insufficient information to determine flip
			{
				// Use relative velocity in normal direction to determine "correct" normal orientation
				//
				//double relVel = (((1.0 - weights[0]) * v[0] + weights[1] * v[1]) - (1.0 - weights[2]) * v[2] + weights[3] * v[3]).dot(n);
				double relVel = ( (weights[0] * v[0] + weights[1] * v[1]) - (weights[2] * v[2] + weights[3] * v[3]) ).dot(n);
				if (relVel > 0.0)
					flip = -1.0;
			}
                        Vector3d rel_vel = ( (weights[0]*v[0] + weights[1]*v[1]) - (weights[2]*v[2] + weights[3]*v[3]) );

			// The four edge vectors
			//
			Vector3d e[4];
                        e[0] = r[3] - r[2];
			e[1] = r[2] - r[3];
                        e[2] = r[0] - r[1];
			e[3] = r[1] - r[0];

			Vector3d t[4];
                        if(1)
                        //if(weights[0]<0 || weights[1]<0 || weights[2]<0 || weights[3]<0)
                        {
                          // finite difference for time gradient
                          for (size_t i=0; i<4; ++i)
                          {
                            for (size_t j=0; j<3; ++j)
                            { 
                              //const double deps = 1e-8;
                              const double deps = 1e-12 + fabs(v[i][j]*1e-8);
                              Vector3d veps[4];
                              for (size_t ii=0; ii<4; ++ii)
                              {
                                for (size_t jj=0; jj<3; ++jj)
                                {
                                  veps[ii][jj] = v[ii][jj];
                                }
                              }
                              veps[i][j] = veps[i][j] + deps;
                              double timeeps = 1.0;
                              double v_partialeps = 0.0;
                              double V1eps = getEEVolume(x[0], x[1], x[2], x[3], veps[0], veps[1], veps[2], veps[3], veps[whichIndex], A, timeeps, time_constant, v_partialeps);
                              t[i][j] = (timeeps - time)/deps;
                            }
                          }
                          // end of finite difference for time gradient
                          /*
                          cout<<"weight[0]: "<<weights[0]<<endl;
                          cout<<"r[0]: "<<r[0]<<endl;
                          cout<<"weight[1]: "<<weights[1]<<endl;
                          cout<<"r[1]: "<<r[1]<<endl;
                          cout<<"weight[2]: "<<weights[2]<<endl;
                          cout<<"r[2]: "<<r[2]<<endl;
                          cout<<"weight[3]: "<<weights[3]<<endl;
                          cout<<"r[3]: "<<r[3]<<endl;
                          */
                        }
                        else
                        {
                          PrismVolume::getEETimeGradient(x[0], x[1], x[2], x[3], v[0], v[1], v[2], v[3], weights, time, m_thickness, t);
                        }
#ifdef DEBUGVG
                        // finite difference for time gradient
                        Vector3d teps[4];
                        for (size_t i=0; i<4; ++i)
                        {
                          for (size_t j=0; j<3; ++j)
                          { 
                            const double deps = 1e-8;
                            Vector3d veps[4];
                            for (size_t ii=0; ii<4; ++ii)
                            {
                              for (size_t jj=0; jj<3; ++jj)
                              {
                                veps[ii][jj] = v[ii][jj];
                              }
                            }
                            veps[i][j] = veps[i][j] + deps;
                            double timeeps = 1.0;
                            double v_partialeps = 0.0;
                            double V1eps = getEEVolume(x[0], x[1], x[2], x[3], veps[0], veps[1], veps[2], veps[3], veps[whichIndex], A, timeeps, time_constant, v_partialeps);
                            teps[i][j] = (timeeps - time)/deps;
                          }
                          cout<<"EE teps["<<i<<"] norm is : "<<teps[i].norm()<<endl;
                          cout<<"EE t["<<i<<"] norm is : "<<t[i].norm()<<endl;
                          cout<<"EE t["<<i<<"] difference norm is: "<<(teps[i] - t[i]).norm()<<endl;
                          for (size_t ii=0; ii<4; ++ii)
                            cout<<"EE weight of vertex: "<<ii<<" is: "<<weights[ii]<<endl;
                        }
                        // end of finite difference for time gradient
#endif                        
			Vector3d w[4][4];
			PrismVolume::getEEBarycentricGradients(r[0], r[1], r[2], r[3], w);
#ifdef DEBUGVG
                        // finite difference for w gradient
                        Vector3d weps[4][4];
			for (size_t i=0; i<4; ++i)
                        {
                                for (size_t j=0; j<3; ++j)
                                {
                                  double weightseps[4];
                                  const double deps = 1e-8;
                                  Vector3d reps[4];
                                  for (size_t ii=0; ii<4; ++ii)
                                  {
                                    for (size_t jj=0; jj<3; ++jj)
                                    {
                                      reps[ii][jj] = r[ii][jj];
                                    }
                                  }
                                  reps[i][j] += deps;
                                  PrismVolume::getBarycentricCoordinates(reps[0], reps[1], reps[2], reps[3], false, weightseps);
                                  
                                  for (size_t ii=0; ii<4; ++ii)
                                    weps[ii][i][j] = (weightseps[ii] - weights[ii]) / deps;
                                }
                        }
                        for (size_t i=0; i<4; ++i)
                        {
                          for (size_t j=0; j<4; ++j)
                          {
                            cout<<"EE w component: "<<i<<", "<<j<<". Difference norm: "<<(weps[i][j]-w[i][j]).norm()<<endl;
                          }
                        }
                        // end of finite difference for w gradient
#endif
			double weight = weights[whichIndex];
			Vector3d vGradient[4];
			for (size_t i=0; i<4; ++i)
			{
				Matrix3d N = Matrix3d::Zero();
				for (size_t j=0; j<4; ++j)
				{
					if (i == j)
						N += ((e[j].cross(n) / l) * n.transpose()) * (time * Matrix3d::Identity() + t[i] * v[j].transpose());
					else
						N += ((e[j].cross(n) / l) * n.transpose()) * (t[i] * v[j].transpose());
				}
#ifdef DEBUGVG
                                // finite difference for dn/dx
				Matrix3d Neps = Matrix3d::Zero();
				for (size_t j=0; j<3; ++j)
                                {
                                  const double deps = 1e-8;
                                  Vector3d veps[4];
                                  for (size_t ii=0; ii<4; ++ii)
                                  {
                                    for (size_t jj=0; jj<3; ++jj)
                                    {
                                      veps[ii][jj] = v[ii][jj];
                                    }
                                  }
                                  veps[i][j] += deps;
                                  
                                  double timeeps = 1.0;
                                  double v_partialeps = 0.0;
			          double V1eps = getEEVolume(x[0], x[1], x[2], x[3], veps[0], veps[1], veps[2], veps[3], veps[whichIndex], A, timeeps, time_constant, v_partialeps);
                                  Vector3d reps[4];
                                  for (size_t ii=0; ii<4; ++ii)
                                          reps[ii] = x[ii] + timeeps * veps[ii];
                                  Vector3d normaleps = (reps[1] - reps[0]).cross(reps[3] - reps[2]);
                                  Vector3d neps(normaleps);
                                  neps.normalize();
                                  Vector3d dndx = (neps - n)/deps;
                                  for (size_t ii=0; ii<3; ++ii)
                                    Neps(ii,j) = dndx(ii);
                                }
                                cout<<"EE N: \n"<<N<<endl;
                                cout<<"EE Neps: \n"<<Neps<<endl;
                                cout<<"EE norm dndx difference: "<<(Neps-N).norm()<<endl;
                                // end of finite difference for dn/dx
#endif
				Vector3d wGrad = Vector3d::Zero();
                                /*
				for (size_t j=0; j<4; ++j)
				{
					if (i == j)
						wGrad += ((t[i] * v[j].transpose()) + time * Matrix3d::Identity()).transpose() * w[whichIndex][j];
					else
						wGrad += (t[i] * v[j].transpose()).transpose() * w[whichIndex][j];
				}
                                */
                                for (size_t j=0; j<4; ++j)
				{
					if (i == j)
						wGrad += ((v[j] * t[i].transpose()) + time * Matrix3d::Identity()).transpose() * w[whichIndex][j];
					else
						wGrad += (v[j] * t[i].transpose()).transpose() * w[whichIndex][j];
				}
#ifdef DEBUGVG
                                // finite difference for weight gradient wGrad
				Vector3d wGradeps = Vector3d::Zero();
				for (size_t j=0; j<3; ++j)
                                {
                                  double weightseps[4];
                                  const double deps = 1e-8;
                                  Vector3d veps[4];
                                  for (size_t ii=0; ii<4; ++ii)
                                  {
                                    for (size_t jj=0; jj<3; ++jj)
                                    {
                                      veps[ii][jj] = v[ii][jj];
                                    }
                                  }
                                  veps[i][j] += deps;
			          double timeeps = 1.0;
                                  double v_partialeps = 0.0;
			          double V1eps = getEEVolume(x[0], x[1], x[2], x[3], veps[0], veps[1], veps[2], veps[3], veps[whichIndex], A, timeeps, time_constant, v_partialeps);
                                  Vector3d reps[4];
                                  for (size_t ii=0; ii<4; ++ii)
                                          reps[ii] = x[ii] + timeeps * veps[ii];
                                  PrismVolume::getBarycentricCoordinates(reps[0], reps[1], reps[2], reps[3], false, weightseps);
                                  wGradeps[j] = (weightseps[whichIndex]-weights[whichIndex])/deps;
                                  cout<<"EE testing wGrad: "<<whichIndex<<" -- "<<(weightseps[whichIndex]-weights[whichIndex])<<endl;
                                }
                                cout<<"EE wGradeps and wGrad difference norm: "<<(wGradeps-wGrad).norm()<<endl;
                                cout<<"EE wGradeps \n"<<wGradeps<<endl;
                                cout<<"EE wGrad \n"<<wGrad<<endl;
                                // end of finite difference for weight gradient wGrad
#endif

                                vGradient[i] = t[i]*v_partial_time + v_partial/v_partial_time*(time-1.0)*N.transpose()*v[whichIndex];

				if (whichIndex == i)
                                {
					//vGradient[i] += v_partial/v_partial_time*(time-1.0)*n;
                                        Vector3d normal_tmp;
                                        normal_tmp[0] = mesh.getVertex(index[i])->normal[0];
                                        normal_tmp[1] = mesh.getVertex(index[i])->normal[1];
                                        normal_tmp[2] = mesh.getVertex(index[i])->normal[2];
                                        if(normal_tmp.dot(n)<=0)
                                            vGradient[i] -= fabs(v_partial/v_partial_time*(time-1.0))*n;
                                        else
                                            vGradient[i] += fabs(v_partial/v_partial_time*(time-1.0))*n;
                                }

                                /*
                                Vector3d normal_tmp;
                                normal_tmp[0] = mesh.getVertex(index[i])->normal[0];
                                normal_tmp[1] = mesh.getVertex(index[i])->normal[1];
                                normal_tmp[2] = mesh.getVertex(index[i])->normal[2];
                                //if(whichIndex == i)
                                    //std::cout<<"normal dot grad: "<<normal_tmp.dot(vGradient[i])<<std::endl;
                                if( (normal_tmp.dot(vGradient[i])<=0) && (whichIndex==i) )
                                    //vGradient[i] = -vGradient[i];
                                    vGradient[i] = normal_tmp;
                                */
                                // gradient direction is opposite to relative speed
                                //if(i == 0 || i == 1)
                                //{
                                //  if(vGradient[i].dot(rel_vel) > 0)
                                //    vGradient[i] = -vGradient[i];
                                //}
                                //else
                                //{
                                //  if(vGradient[i].dot(rel_vel) < 0)
                                //    vGradient[i] = -vGradient[i];
                                //}
			}

			for (size_t i=0; i<4; ++i)
				for (size_t j=0; j<3; ++j)
                                {
				        gradients.coeffRef(3 * index[i] + j, 0) += vGradient[i][j];
				        //gradients.coeffRef(3 * index[i] + j, 0) += flip * vGradient[i][j] * A;
                                        //if(i==whichIndex)
					  //gradients.coeffRef(3 * index[i] + j, 0) += mesh.getVertex(index[i])->normal[j]*A;
#ifdef DEBUGVG
                                        // finite difference for vGradient 
                                        double weightseps[4];
                                        const double deps = 1e-8;
			                Vector3d veps[4];
                                        for (size_t ii=0; ii<4; ++ii)
                                        {
                                          for (size_t jj=0; jj<3; ++jj)
                                          {
                                            veps[ii][jj] = v[ii][jj];
                                          }
                                        }
                                        veps[i][j] += deps;

                                        double timeeps = 1.0;
                                        double v_partialeps;
                                        Vector3d reps[4];
                                        double V1eps = getEEVolume(x[0], x[1], x[2], x[3], veps[0], veps[1], veps[2], veps[3], veps[whichIndex], A, timeeps, time_constant, v_partialeps);
                                        for (size_t ii=0; ii<4; ++ii)
                                          reps[ii] = x[ii] + timeeps * veps[ii];
                                        PrismVolume::getBarycentricCoordinates(reps[0], reps[1], reps[2], reps[3], false, weightseps);
                                        V1eps *= weightseps[whichIndex];
                                        V1eps = flip*V1eps;
                                        double V1tmp = flip*V1;
                                        double grad_diff = (V1eps - V1tmp)/deps;
                                        /*
                                        cout<<"EE V1eps --- "<<V1eps<<endl;
                                        cout<<"EE V1tmp --- "<<V1tmp<<endl;
                                        cout<<"EE V1 ------ "<<V1<<endl;
                                        cout<<"EE weighteps ------ "<<weightseps[whichIndex]<<endl;
                                        cout<<"EE weight ------ "<<weights[whichIndex]<<endl;
                                        cout<<"EE volume gradient respect to vertex: "<<i<<" ,component: "<<j<<" is: "<<grad_diff<<" --- finite diff"<<endl; 
                                        cout<<"EE volume gradient respect to vertex: "<<i<<" ,component: "<<j<<" is: "<<(flip * vGradient[i][j] * A)<<" --- analytic"<<endl; 
                                        cout<<"EE collision time respect to vertex: "<<i<<" ,component: "<<j<<" is: "<<timeeps<<" --- finite diff"<<endl; 
                                        cout<<"EE collision time respect to vertex: "<<i<<" ,component: "<<j<<" is: "<<time<<" --- analytic"<<endl; 
                                        cout<<"EE diff of volume gradient respect to vertex: "<<i<<" ,component: "<<j<<" is: "<<grad_diff - (flip * vGradient[i][j] * A)<<endl; 
                                        for (size_t ii=0; ii<4; ++ii)
                                          cout<<"EE weight of vertex: "<<ii<<" is: "<<weights[ii]<<endl;
                                        */
                                        //gradients.coeffRef(3 * index[i] + j, 0) += grad_diff;
                                        // end of finite difference for vGradient
#endif
                                }
		}
	}

	m_volume *= 1.01;
        m_volume = min(m_volume,-1e-10);

	gradients.finalize();
	m_gradients = SparseVector<double>(gradients);
}

double IntersectionVolume::getVoronoiArea(Polyhedron &mesh, unsigned int idx)
{
	double A = 0.0;
	Halfedge_Vertex_circulator hvcir1 = mesh.getVertex(idx)->vertex_begin();
	Halfedge_Vertex_circulator hvcir2 = mesh.getVertex(idx)->vertex_begin(); ++hvcir2;
	do
	{
		Vector3d v1 = hvcir1->opposite()->vertex()->Pori() - mesh.getVertex(idx)->Pori();
		A += v1.cross((hvcir2->opposite()->vertex()->Pori() - mesh.getVertex(idx)->Pori())).norm();
		++hvcir1;
		++hvcir2;
	} while (hvcir1 != mesh.getVertex(idx)->vertex_begin());
	A = A / 6.0;

	//return 1.0;
	return A;
}

double IntersectionVolume::getVolume(Vector3d &x0,Vector3d &x1,Vector3d &x2,Vector3d &x3,
								     Vector3d &v0,Vector3d &v1,Vector3d &v2,Vector3d &v3, Vector3d &v, double A, double &t)
{
	double V = 0.0;

	if (m_thickness)
	{
		if (!PrismVolume::VertexFaceProximity(x0,x1,x2,x3,v0,v1,v2,v3,m_thickness,t))
                {
			t = 1.0;
                        cout<<"finite diff no time root for v-f with thickness!!!!"<<endl;
                }
	}
	else
	{
		if (!PrismVolume::VertexFaceIntersection(x0,x1,x2,x3,v0,v1,v2,v3,t))
                {
			t = 1.0;
                        cout<<"finite diff no time root for v-f!!!!"<<endl;
                }
	}
	
	Vector3d x0t = x0 + v0 * t;
	Vector3d x1t = x1 + v1 * t;
	Vector3d x2t = x2 + v2 * t;

	Vector3d n = (x1t - x0t).cross(x2t - x0t);
	n.normalize();

	double d = v.dot(n);
        //std::cout<<"d: "<<d<<". vnorm: "<<v.norm()<<std::endl;

        //double t_tmp = max(t - 0.01, 0.0);
        double t_tmp = max(t - 0.01, 0.0);
        V = (1.0 - t_tmp) * d * A;

	//V = (1.0 - t) * d * A;
	//V = d * A;
	
	return V;
}

double IntersectionVolume::getVolume(Vector3d &x0,Vector3d &x1,Vector3d &x2,Vector3d &x3,
								     Vector3d &v0,Vector3d &v1,Vector3d &v2,Vector3d &v3, Vector3d &v, double A, double &t, double time_constant, double &v_partial)
{
	double V = 0.0;

	if (m_thickness)
	{
		if (!PrismVolume::VertexFaceProximity(x0,x1,x2,x3,v0,v1,v2,v3,m_thickness,t))
                {
			t = 1.0;
                        cout<<"finite diff no time root for v-f with thickness!!!!"<<endl;
                }
	}
	else
	{
		if (!PrismVolume::VertexFaceIntersection(x0,x1,x2,x3,v0,v1,v2,v3,t))
                {
			t = 1.0;
                        cout<<"finite diff no time root for v-f!!!!"<<endl;
                }
	}
	
	Vector3d x0t = x0 + v0 * t;
	Vector3d x1t = x1 + v1 * t;
	Vector3d x2t = x2 + v2 * t;

	Vector3d n = (x1t - x0t).cross(x2t - x0t);
	n.normalize();

	double d = v.dot(n);

        v_partial = d*A;
        double t_tmp = max(t - 0.01, 0.0);
        V = (1.0 - t_tmp) * sqrt(time_constant*time_constant + v_partial*v_partial);
	
	return V;
}

double IntersectionVolume::getEEVolume(Vector3d &x0,Vector3d &x1,Vector3d &x2,Vector3d &x3,
								       Vector3d &v0,Vector3d &v1,Vector3d &v2,Vector3d &v3, Vector3d &v, double A, double &t)
{
	double V = 0.0;

	if (m_thickness)
	{
		//assert(0);
		if (!PrismVolume::EdgeEdgeProximity(x0,x1,x2,x3,v0,v1,v2,v3,m_thickness,t))
                {
			t = 1.0;
                        cout<<"finite diff no time root for e-e with thickness!!!!"<<endl;
                }
	}
	else
	{
		if (!PrismVolume::EdgeEdgeIntersection(x0,x1,x2,x3,v0,v1,v2,v3,t))
                {
			t = 1.0;
                        cout<<"finite diff no time root for e-e!!!!"<<endl;
                }
	}
	
	Vector3d x0t = x0 + v0 * t;
	Vector3d x1t = x1 + v1 * t;
	Vector3d x2t = x2 + v2 * t;
	Vector3d x3t = x3 + v3 * t;

	Vector3d n = (x1t - x0t).cross(x3t - x2t);
	n.normalize();

	double d = v.dot(n);
        /*
        cout<<"v is ---------\n"<<v<<endl;
        cout<<"v norm ---------\n"<<v.norm()<<endl;
        cout<<"n is ---------\n"<<n<<endl;
        cout<<"n norm ---------\n"<<n.norm()<<endl;
        cout<<"d is ---------"<<d<<endl;
        */

        //double t_tmp = max(t - 0.01, 0.0);
        double t_tmp = max(t - 0.01, 0.0);
        V = (1.0 - t_tmp) * d * A;
        
        //cout<<"t_tmp is ---------"<<t_tmp<<endl;
        //cout<<"V is ---------"<<V<<endl;

	//V = (1.0 - t) * d * A;
	//V = d * A;
	
	return V;
}

double IntersectionVolume::getEEVolume(Vector3d &x0,Vector3d &x1,Vector3d &x2,Vector3d &x3,
								       Vector3d &v0,Vector3d &v1,Vector3d &v2,Vector3d &v3, Vector3d &v, double A, double &t, double time_constant, double &v_partial)
{
	double V = 0.0;

	if (m_thickness)
	{
		//assert(0);
		if (!PrismVolume::EdgeEdgeProximity(x0,x1,x2,x3,v0,v1,v2,v3,m_thickness,t))
                {
			t = 1.0;
                        cout<<"finite diff no time root for e-e with thickness!!!!"<<endl;
                }
	}
	else
	{
		if (!PrismVolume::EdgeEdgeIntersection(x0,x1,x2,x3,v0,v1,v2,v3,t))
                {
			t = 1.0;
                        cout<<"finite diff no time root for e-e!!!!"<<endl;
                }
	}
	
	Vector3d x0t = x0 + v0 * t;
	Vector3d x1t = x1 + v1 * t;
	Vector3d x2t = x2 + v2 * t;
	Vector3d x3t = x3 + v3 * t;

	Vector3d n = (x1t - x0t).cross(x3t - x2t);
	n.normalize();

	double d = v.dot(n);

        v_partial = d*A;
        double t_tmp = max(t - 0.01, 0.0);
        V = (1.0 - t_tmp) * sqrt(time_constant*time_constant + v_partial*v_partial);
	
	return V;
}

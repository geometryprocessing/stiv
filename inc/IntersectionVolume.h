// IntersectionVolume.h
//

#ifndef _INTERSECTIONVOLUME_H_
#define _INTERSECTIONVOLUME_H_

#include <vector>
#include <set>
#include <map>

#include "Collision.h"
#include "Common.h"

#include <CGAL/Point_3.h>

#include <Eigen/Core>
#include <Eigen/Sparse>

typedef std::multimap<unsigned int, CollisionsIterator> CollisionMap;
typedef std::multimap<unsigned int, CollisionsIterator>::iterator CollisionMapIterator;

class ColSurface;



class IntersectionVolume
{
public:
    IntersectionVolume(Polyhedron &mesh, double h, Collisions &collisions, double periodicLength, int nv);

    double getVolume()
    { return m_volume; }
    
    Eigen::SparseVector<double>& getGradients()
    { return m_gradients; }

    int m_bdv;
protected:

	double getNormal(Polyhedron &mesh, Collision &c, unsigned int index[4], unsigned int vIndex, Eigen::Vector3d &n);
	void getNormalGradient(Polyhedron &mesh, Collision &c, unsigned int index[4], unsigned int vIndex, bool flip, Eigen::Vector3d &n, double l, Eigen::Vector3d &nGrad);
	void getNormalGradient(Eigen::Vector3d &r0, Eigen::Vector3d &r1, Eigen::Vector3d &r2, bool flip, Eigen::Matrix3d N[3]);

	void findFirstCollisions(Polyhedron &mesh, Collisions &collisions, CollisionMap &sorted);
	void findFirstCollisionsApprox(Polyhedron &mesh, Collisions &collisions, CollisionMap &sorted);

	void computeVolume(Polyhedron &mesh, CollisionMap &collisions);
	void computeVolumeApprox(Polyhedron &mesh, CollisionMap &collisions);
	void computeVolumeNew(Polyhedron &mesh, CollisionMap &collisions);

	double getVoronoiArea(Polyhedron &mesh, unsigned int idx);

	void getVolumeGradient(Polyhedron &mesh, Eigen::Vector3d &v, double oldVolume, double A, Collision &c, unsigned int idx, bool perturb, Eigen::Vector3d &gradient);

	double getVolume(Eigen::Vector3d &x0,Eigen::Vector3d &x1,Eigen::Vector3d &x2,Eigen::Vector3d &x3,
				     Eigen::Vector3d &v0,Eigen::Vector3d &v1,Eigen::Vector3d &v2,Eigen::Vector3d &v3,
					 Eigen::Vector3d &v, double A, double &t);
        
        double getVolume(Eigen::Vector3d &x0,Eigen::Vector3d &x1,Eigen::Vector3d &x2,Eigen::Vector3d &x3,
				     Eigen::Vector3d &v0,Eigen::Vector3d &v1,Eigen::Vector3d &v2,Eigen::Vector3d &v3,
					 Eigen::Vector3d &v, double A, double &t, double time_constant, double &v_partial);

	double getEEVolume(Eigen::Vector3d &x0,Eigen::Vector3d &x1,Eigen::Vector3d &x2,Eigen::Vector3d &x3,
					   Eigen::Vector3d &v0,Eigen::Vector3d &v1,Eigen::Vector3d &v2,Eigen::Vector3d &v3, Eigen::Vector3d &v, double A, double &t);
        
        double getEEVolume(Eigen::Vector3d &x0,Eigen::Vector3d &x1,Eigen::Vector3d &x2,Eigen::Vector3d &x3,
					   Eigen::Vector3d &v0,Eigen::Vector3d &v1,Eigen::Vector3d &v2,Eigen::Vector3d &v3, Eigen::Vector3d &v, double A, double &t, double time_constant, double &v_partial);

        double m_volume;

	double m_thickness;
        double m_periodicLength;

	Eigen::SparseVector<double> m_gradients;

        Eigen::VectorXd m_handleGradients;
	
	Eigen::Vector3d m_x_cm;
        int m_nv;

private:

};

typedef std::vector<IntersectionVolume> IntersectionVolumes;
typedef std::vector<IntersectionVolume>::iterator IntersectionVolumesIterator;

#endif


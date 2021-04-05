#pragma once
#include <vector>
#include "Common.h"
#include "ColSurface.h"
#include <CGAL/bounding_box.h>

typedef enum
{
	//SURF_SUBDIVISION,
	//SURF_LAPLACIAN,
	//SURF_ARAP,
	SURF_SIMPLE,
	//SURF_FASTLSM
} SURF_TYPE;

/**
 * Mesh manager manages all meshes created, which one is active.
 * There is no reason to have more than one manager, so I am using singleton
 * pattern here.
 */
class /*STRUCTURES_API*/ MeshManager
{
public:
        MeshManager() 
	{
            s = 0;
	}
        ~MeshManager() 
        {
            if(s)
                delete s;
        }
public:
	ColSurface& setMesh(Polyhedron *p, SURF_TYPE st, int nv, int order=0, int np=0, int pp=0);
        ColSurface* getMesh();
private:
	ColSurface* s;
};

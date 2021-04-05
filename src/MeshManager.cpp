#include <boost/foreach.hpp>
#include "MeshManager.h"

#include <CGAL/Subdivision_method_3.h>

//#include "LoopSubdivisionSurface.h"
//#include "LaplacianEditing.h"
//#include "ArapEditing.h"
#include "SimpleSurface.h"
//#include "FFDEditing.h"

//MeshManager* MeshManager::instance = NULL;

ColSurface& MeshManager::setMesh(Polyhedron *p, SURF_TYPE st, int nv, int order, int np, int pp) 
{
	
	if (s)
		delete s;
		
	switch (st) 
	{
		//case SURF_SUBDIVISION:
		//	s = new LoopSubdivisionSurface(p);
		//	break;
                /*
		case SURF_LAPLACIAN:
			s = new LaplacianEditing(p);
			break;
		case SURF_ARAP:
			s = new ArapEditing(p);
			break;
                */
		case SURF_SIMPLE:
			s = new SimpleSurface(p, nv, order, np, pp);
			break;
                /*
		case SURF_FASTLSM:
			s = new FFDEditing(p);
			break;
                */ 
		default:
			assert(0);
			break;
	}
	
    return *s;
}

ColSurface* MeshManager::getMesh() 
{
	return s;
}

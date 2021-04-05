// SimpleSurface.h
//

#ifndef SIMPLESURFACE_H_
#define SIMPLESURFACE_H_

#include "Common.h"
#include "ColSurface.h"

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

class SimpleSurface : public ColSurface
{
public:
	
	SimpleSurface(Polyhedron *p, int nv, int order=0, int np=0, int pp=0);
	~SimpleSurface() { }
protected:
	void initialize();
};

#endif


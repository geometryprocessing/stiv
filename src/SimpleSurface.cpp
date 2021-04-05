// SimpleSurface.h
//

#include "SimpleSurface.h"



SimpleSurface::SimpleSurface(Polyhedron *p, int nv, int order, int np, int pp)
 : ColSurface(p)
{
        m_nv = nv;
        m_p = order;
        m_np = np;
        m_pp = pp;
	initialize();
}

void SimpleSurface::initialize()
{
	m_slaveMesh = m_controlMesh;
	ColSurface::initialize();
	computeNormals(m_slaveMesh);
}

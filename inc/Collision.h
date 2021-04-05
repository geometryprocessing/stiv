// Collision.h
//

#ifndef _COLLISION_H_
#define _COLLISION_H_

#include <vector>

#include "Common.h"

class Collision
{
public:
    Collision(unsigned int v0,
              unsigned int v1, unsigned int v2, unsigned int v3,
              double t, bool vf)
	 : m_v0(v0), m_v1(v1), m_v2(v2), m_v3(v3),
       m_time(t), m_vf(vf)
    { }

	Collision()
	 : m_vf(true)
    { }

	bool isVF() const
	{ return m_vf; }

	bool isEE() const
	{ return !m_vf; }
	
	double getTime() const
	{ return m_time; }

	void set(unsigned int v0, unsigned int v1, unsigned int v2, unsigned int v3)
	{ m_v0 = v0; m_v1 = v1; m_v2 = v2; m_v3 = v3; }

	unsigned int getFirstEdgeVertex1() const
	{ return m_v0; }
	
	unsigned int getFirstEdgeVertex2() const
	{ return m_v1; }
	
	unsigned int getSecondEdgeVertex1() const
	{ return m_v2; }
	
	unsigned int getSecondEdgeVertex2() const
	{ return m_v3; }

    unsigned int getVertex() const
    { return m_v0; }
    
    unsigned int getTriangleVertex1() const
    { return m_v1; }
    
    unsigned int getTriangleVertex2() const
    { return m_v2; }
    
    unsigned int getTriangleVertex3() const
    { return m_v3; }

	void setFace(unsigned int face)
	{ m_face = face; }

	unsigned int getFace() const
	{ return m_face; }

protected:
    unsigned int m_v0;
    unsigned int m_v1;
    unsigned int m_v2;
    unsigned int m_v3;

	unsigned int m_face;

	bool m_vf;

    double m_time;

private:

};

typedef std::vector<Collision> Collisions;
typedef std::vector<Collision>::iterator CollisionsIterator;

class CollisionEdge
{
public:
    CollisionEdge(unsigned int e1, unsigned int e2, // edge indices
              double t)
	: m_e1(e1), m_e2(e2), m_time(t)
    { }
	
	CollisionEdge()
    { }
	
	double getTime() const
	{ return m_time; }
	
    unsigned int getEdge1() const
    { return m_e1; }
    
    unsigned int getEdge2() const
    { return m_e2; }
	
protected:
    unsigned int m_e1;
	unsigned int m_e2;
    double m_time;
	
private:
	
};

typedef std::vector<CollisionEdge> CollisionEdges;
typedef std::vector<CollisionEdge>::iterator CollisionEdgesIterator;

#endif


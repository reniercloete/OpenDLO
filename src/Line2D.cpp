// Line2D.cpp
// Copyright (c) 2022, Renier Cloete
// This program is released under the BSD license. See the file LICENSE.txt for details.

#include "Line2D.h"

#include <vector>
#include <algorithm>
#include "Vector2D.h"
#include "Constants.h"

//required input point must be colinear with the line
bool on_segment( const CVector2D& p, const CLine2D& l )
{
	//if a point is on the line, the sum of the vectors formed by the point to the line endpoints must be equal
	CVector2D va = p - l.P1();
	CVector2D vb = p - l.P2();
	double ma = va.Length();
	double mb = vb.Length();
	double ml = (l.P2() - l.P1()).Length();
	double s = ma + mb;
	bool r = s <= ml + EPSILON;
	return r;
}

void MakeUnique( std::vector<CVector2D>& List )
{
	size_t i = 0;
	while ( i < List.size() )
	{
		size_t j = i + 1;
		while ( j < List.size() )
		{
			if ( (List[i] - List[j]).LengthSquared() < EPSILON )
				List.erase( List.begin() + j );
			else
				++j;
		}
		++i;
	}
}

bool Colinear( const CLine2D& la,
			   const CLine2D& lb )
{
	/*CVector2D oa, ob, da, db; //origin and direction vectors
	double sa, sb; //scalar values
	oa = la.P1();
	da = la.P2() - la.P1();
	ob = lb.P1();
	db = lb.P2() - lb.P1();

	return (abs( da*db ) < EPSILON && abs( (ob - oa)*da ) < EPSILON);*/

	const CVector2D& va = la.Vector();
	const CVector2D& vb = lb.Vector();

	return
		( abs( va.x - vb.x ) < EPSILON &&
		  abs( va.y - vb.y ) < EPSILON) ||
		( abs( va.x + vb.x ) < EPSILON &&
		  abs( va.y + vb.y ) < EPSILON);
}

//compute using vector math
//returns 0 points if the lines do not intersect or overlap
//returns 1 point if the lines intersect
//returns 2 points if the lines overlap, contain the points where overlapping start starts and stop
void Intersect( const CLine2D& la, 
				const CLine2D& lb,
				std::vector<CVector2D>& Points )
{
	Points.clear();

#define LEPSILON 11e-12

	CVector2D oa, ob, da, db; //origin and direction vectors
	double sa, sb; //scalar values
	oa = la.P1();
	da = la.P2() - la.P1();
	ob = lb.P1();
	db = lb.P2() - lb.P1();

	if ( abs( da*db ) < LEPSILON && abs( (ob - oa)*da ) < LEPSILON ) //if colinear
	{
		if ( on_segment( lb.P1(), la ) && on_segment( lb.P2(), la ) )
		{
			Points.push_back( lb.P1() );
			Points.push_back( lb.P2() );
			MakeUnique( Points );
			return;
		}

		if ( on_segment( la.P1(), lb ) && on_segment( la.P2(), lb ) )
		{
			Points.push_back( la.P1() );
			Points.push_back( la.P2() );
			MakeUnique( Points );
			return;
		}

		if ( on_segment( la.P1(), lb ) )
			Points.push_back( la.P1() );

		if ( on_segment( la.P2(), lb ) )
			Points.push_back( la.P2() );

		if ( on_segment( lb.P1(), la ) )
			Points.push_back( lb.P1() );

		if ( on_segment( lb.P2(), la ) )
			Points.push_back( lb.P2() );

		MakeUnique( Points );
		return;
	}

	if ( abs( da*db ) < LEPSILON && abs( (ob - oa)*da ) > LEPSILON )
	{
		MakeUnique( Points );
		return;
	}

	//math trick db cross db == 0, which is a single scalar in 2D
	//crossing both sides with vector db gives
	sa = ((ob - oa)*db) / (da*db);

	//crossing both sides with vector da gives
	sb = ((oa - ob)*da) / (db*da);

	if ( (0 < sa || abs( sa ) < LEPSILON) && (sa < 1 || abs( sa - 1 ) < LEPSILON) &&
		 (0 < sb || abs( sb ) < LEPSILON) && (sb < 1 || abs( sb - 1 ) < LEPSILON) )
	{
		Points.push_back( oa + da * sa );
		MakeUnique( Points );
		return;
	}

	MakeUnique( Points );
	return;
}

CLine2D::CLine2D() :
	mp1( 0, 0 ), 
	mp2( 0, 0 )
{
	Update();
}

CLine2D::CLine2D( const CPoint2D& p1, 
				  const CPoint2D& p2 ) :
	mp1( p1 ), 
	mp2( p2 ),
	mSlope( 0 )
{
	Update();
};

void 
CLine2D::Set( const CPoint2D& p1,
			  const CPoint2D&  p2,
			  bool DoUpdate )
{
	mp1 = p1;
	mp2 = p2;
	mSlope = 0;

	if ( DoUpdate )
		Update();
}

void 
CLine2D::Update()
{
	if ( P1().x < P2().x )
		mMin.x = P1().x;
	else
		mMin.x = P2().x;
	if ( P1().y < P2().y )
		mMin.y = P1().y;
	else
		mMin.y = P2().y;

	if ( P1().x > P2().x )
		mMax.x = P1().x;
	else
		mMax.x = P2().x;
	if ( P1().y > P2().y )
		mMax.y = P1().y;
	else
		mMax.y = P2().y;

	mVector = P2() - P1();
	mVector.Normalize();

	double dx = mp2.x - mp1.x;
	if ( abs( dx ) > LEPSILON )
		mSlope = (mp2.y - mp1.y) / dx;
	else
		mSlope = DBL_MAX;
}

double 
CLine2D::DistanceTo( const CPoint2D& Point ) const
{
	double Result = 0;

	CVector2D v1 = P2() - P1();
	v1.Normalize();

	CVector2D v2 = Point - P1();
	double l = v2.Length();
	v2.Normalize();

	double SinTheta = abs( v2*v1 );
	Result = l*SinTheta;

	if ( abs(Result) < EPSILON )
	{
		double l = (P1() - P2()).Length();
		double l1 = (P1() - Point).Length();
		double l2 = (P2() - Point).Length();
		if ( l1 > l || l2 > l )
			Result = (std::min)( l1, l2 );
	}

	return Result;
}

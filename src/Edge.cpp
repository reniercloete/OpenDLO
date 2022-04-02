// Edge.cpp
// Copyright (c) 2022, Renier Cloete
// This program is released under the BSD license. See the file LICENSE.txt for details.

#include "Edge.h"

#include "Constants.h"

#include <assert.h>

std::vector<CNode>* CEdge::mNodes = nullptr;

CEdge::CEdge( size_t aN1,
			  size_t aN2,
			  eEdgeType aType,
			  double aLength,
			  double aMpPosx,
			  double aMpNegx,
			  double aMpPosy,
			  double aMpNegy ) :
	N1( aN1 ),
	N2( aN2 ),
	Type( aType ),
	Length( aLength ),
	MpPosx( aMpPosx ),
	MpNegx( aMpNegx ),
	MpPosy( aMpPosy ),
	MpNegy( aMpNegy )
{
	assert( mNodes );

	Line = CLine2D( (*mNodes)[N1 - 1].Point, (*mNodes)[N2 - 1].Point );
	Added = false;
	Removeable = false;
	UDLVectorCalculated = false;
	Delete = false;

	++Counter;
	mID = Counter;
	YieldRatio = 0;

	CVector2D v = ((*mNodes)[N2 - 1].Point - (*mNodes)[N1 - 1].Point) / Length;

	double c = v.Dot( kXAxis );
	double s = v.Dot( kYAxis );

	MpPos = MpPosx*c*c + MpPosy*s*s;
	MpNeg = MpNegx*c*c + MpNegy*s*s;
}

const CPoint2D& 
CEdge::GetPoint( size_t Index )
{
	if (Index == 0)
		return mNodes->at( N1 - 1 ).Point;

	return mNodes->at( N2 - 1 ).Point;
}

void 
CEdge::SetNodes( std::vector<CNode>* Nodes )
{
	mNodes = Nodes;
}

void
CEdge::GetUDLLoadVector( std::array<double, 3>& UDLVector,
						 const CPoly2D& Outline )
{
	if ( !UDLVectorCalculated )
	{
		mUDLVector[0] = 0;
		mUDLVector[1] = 0;
		mUDLVector[2] = 0;

		CVector2D Ray( 0, 1e6 );
		std::vector<CPoint2D> Intersections;
		CPoint2D p, c;

		CPoint2D p1 = (*mNodes)[N1 - 1].Point;
		CPoint2D p2 = (*mNodes)[N2 - 1].Point;

		if ( abs( p1.x - p2.x ) > EPSILON )
		{
			bool left = false;
			if ( p1.x > p2.x )
			{
				left = true;
				std::swap( p1, p2 );
			}

			CLine2D Line1( p1 - Ray, p1 + Ray );
			CLine2D Line2( p2 - Ray, p2 + Ray );
			CLine2D Line3( p1, p2 );

			std::vector<CPoly2D> Polies = Outline.ClipRight( Line2 );
			Polies = CPoly2D::GetByPoint( (p1 + p2) / 2, Polies ).ClipLeft( Line1 );
			Polies = CPoly2D::GetByPoint( (p1 + p2) / 2, Polies ).ClipRight( Line3 );

			CPoly2D Poly = CPoly2D::GetByPoint( (p1 + p2) / 2, Polies );

			std::vector<CPoint2D> Original = Poly.GetPoints();

			for ( size_t j = 0; j < Original.size(); ++j )
			{
				Line3.Set( Original[j] - Ray, Original[j] + Ray, false );

				Intersections = Poly.IntersectWith( Line3 );
			}

			Original = Poly.GetPoints();
			size_t j = 0;
			while ( j < Original.size() )
			{
				Line3.Set( Original[j] - Ray, Original[j] + Ray, false );

				Intersections = Poly.IntersectWith( Line3 );

				bool Remove = false;
				for ( size_t k = 0; k < Intersections.size() - 1; ++k )
				{
					if ( Intersections[k].y - Original[j].y < 0 )
					{
						p = (Intersections[k] + Intersections[k + 1]) / 2;
						if ( !(Poly.PointOnPoly( p ) > -1) && !Poly.PointInPoly( p ) )
						{
							Remove = true;
							break;
						}
					}
				}
				if ( Remove )
					Original.erase( Original.begin() + j );
				else
					++j;
			}

			if ( Original.size() )
			{
				Poly.SetPoints( Original );
				CPoly2D TestPoly;

				double A = Poly.Area();
				c = Poly.Centroid();

				p1 = (*mNodes)[N1 - 1].Point;
				p2 = (*mNodes)[N2 - 1].Point;

				Line3.Set( p1, p2, false );

				p = (p1 + p2) / 2;

				double dn = Line3.DistanceTo( c );
				double dt = Line3.Vector().Dot( c - p );

				mUDLVector[0] = A*dn;
				mUDLVector[1] = A*dt;
				mUDLVector[2] = A;
			}
		}

		if ( Type != eEdgeType::FREE &&
			 Type != eEdgeType::SYMMETRY )
		{
			mUDLVector[1] = 0;
			mUDLVector[2] = 0;
		}

		UDLVectorCalculated = true;
	}

	UDLVector = mUDLVector;
}

int CEdge::Counter = 0;

int
CEdge::DOF()
{
	if ( Type == eEdgeType::FREE ||
		 Type == eEdgeType::SYMMETRY )
	{
		return 3;
	}

	return 1;
}

void
CEdge::GetCompatibilityMatrix( std::vector<std::vector<double>>& Matrix,
							   bool ApplBoundaryConditions )
{
	CVector2D v = ((*mNodes)[N2 - 1].Point - (*mNodes)[N1 - 1].Point) / Length;

	double c = v.Dot( kXAxis );
	double s = v.Dot( kYAxis );

	//Node1
	Matrix[0][0] = c;
	Matrix[0][1] = -s;
	Matrix[0][2] = 0;

	Matrix[1][0] = s;
	Matrix[1][1] = c;
	Matrix[1][2] = 0;

	Matrix[2][0] = 0;
	Matrix[2][1] = Length / 2;
	Matrix[2][2] = 1;

	//Node2
	Matrix[3][0] = -c;
	Matrix[3][1] = s;
	Matrix[3][2] = 0;

	Matrix[4][0] = -s;
	Matrix[4][1] = -c;
	Matrix[4][2] = 0;

	Matrix[5][0] = 0;
	Matrix[5][1] = Length / 2;
	Matrix[5][2] = -1;

	if ( ApplBoundaryConditions &&
		 Type != eEdgeType::FREE &&
		 Type != eEdgeType::SYMMETRY )
	{
		//Node1
		Matrix[0][1] = 0;
		Matrix[0][2] = 0;

		Matrix[1][1] = 0;
		Matrix[1][2] = 0;

		Matrix[2][1] = 0;
		Matrix[2][2] = 0;

		//Node2
		Matrix[3][1] = 0;
		Matrix[3][2] = 0;

		Matrix[4][1] = 0;
		Matrix[4][2] = 0;

		Matrix[5][1] = 0;
		Matrix[5][2] = 0;
	}
}
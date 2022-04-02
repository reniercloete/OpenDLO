// Poly2D.h
// Copyright (c) 2022, Renier Cloete
// This program is released under the BSD license. See the file LICENSE.txt for details.

#pragma once

#include "Point2D.h"
#include "Line2D.h"
#include "Enums.h"

class CPoly2D
{
public:
	CPoly2D();
	virtual ~CPoly2D() {};

	size_t AddPoint( const CPoint2D& p,
					 bool Check = false );

	CPoint2D& GetPoint( size_t index );
	eEdgeType GetEdgeType( size_t Index );
	void SetEdgeType( size_t Index,
					  eEdgeType Type );

	const CPoint2D& GetPoint( size_t index ) const;

	const std::vector<CPoint2D>& GetPoints() { return mPoints; }
	void SetPoints( const std::vector<CPoint2D>& Points );

	void InsertPoint( int index, const CPoint2D& p );
	size_t GetNumPoints() const;

	std::vector<CPoint2D> IntersectWith( const CLine2D& line );
	void GetOrderedIntersections( const CLine2D& line,
								  std::vector<CPoint2D>& Intersections ) const;

	CPoint2D mMin;
	CPoint2D mMax;

	CPoint2D Centroid() const;
	double Area() const;
	void Reverse();
	void MakeAntiClockWise();

	std::vector<CPoly2D> ClipLeft( const CLine2D& Line ) const;
	std::vector<CPoly2D> ClipRight( const CLine2D& Line ) const;

	bool PointInPoly( const CPoint2D& p ) const;
	int PointOnPoly( const CPoint2D& p ) const;

	static CPoly2D GetByPoint( const CPoint2D& Point,
							   const std::vector<CPoly2D>& Polies );

protected:
	std::vector<CPoint2D> mPoints;
	std::vector<eEdgeType> mEdgeTypes;
};


// Line2D.h
// Copyright (c) 2022, Renier Cloete
// This program is released under the BSD license. See the file LICENSE.txt for details.

#pragma once

#include "Point2D.h"
#include "Vector2D.h"

#include <vector>

class CLine2D
{
public:
	CLine2D();
	CLine2D( const CPoint2D& _p1, const CPoint2D& _p2 );
	virtual ~CLine2D() {};

	double DistanceTo( const CPoint2D& Point ) const;
	const CVector2D& Vector() const
	{
		return mVector;
	};

	const CPoint2D& Min() const
	{
		return mMin;
	};

	const CPoint2D& Max() const
	{
		return mMax;
	};

	void Update();
	void Set( const CPoint2D& p1,
			  const CPoint2D&  p2,
			  bool DoUpate = true );

	const CPoint2D& P1() const
	{
		return mp1;
	};
	const CPoint2D& P2() const
	{
		return mp2;
	};

	double Slope() const { return mSlope; };

private:
	CPoint2D mp1;
	CPoint2D mp2;
	CPoint2D mMax;
	CPoint2D mMin;

	CVector2D mVector;
	double mSlope;
};

bool Colinear( const CLine2D& la,
			   const CLine2D& lb );
void Intersect( const CLine2D& la, 
				const CLine2D& lb,
				std::vector<CVector2D>& Points );

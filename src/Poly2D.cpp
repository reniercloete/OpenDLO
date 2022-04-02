// Poly2D.cpp
// Copyright (c) 2022, Renier Cloete
// This program is released under the BSD license. See the file LICENSE.txt for details.

#include "Poly2D.h"

#include <assert.h>
#include <algorithm>
#include "Constants.h"

/*#include <gl/gl.h>
#include <gl/glu.h>*/

CPoly2D::CPoly2D()
{
	mMin = CPoint2D( DBL_MAX, DBL_MAX );
	mMax = CPoint2D( -DBL_MAX, -DBL_MAX );
}

size_t
CPoly2D::AddPoint( const CPoint2D& p,
				   bool Check )
{
	bool Add = true;
	size_t Result = -1;
	if ( Check )
	{
		double d;
		for ( size_t i = 0; i < mPoints.size(); ++i )
		{
			d = (mPoints[i] - p).LengthSquared();
			if ( d < 1e-20 )
			{
				Add = false;
				Result = i;
				break;
			}
		}
	}

	if ( Add )
	{
		if ( p.x < mMin.x )
			mMin.x = p.x;
		if ( p.x > mMax.x )
			mMax.x = p.x;

		if ( p.y < mMin.y )
			mMin.y = p.y;
		if ( p.y > mMax.y )
			mMax.y = p.y;

		Result = mPoints.size();
		mPoints.push_back( p );
		mEdgeTypes.push_back( eEdgeType::FREE );
	}

	return Result;
}

eEdgeType
CPoly2D::GetEdgeType( size_t Index )
{
	return mEdgeTypes[Index];
}

void
CPoly2D::SetEdgeType( size_t Index,
					  eEdgeType Type )
{
	mEdgeTypes[Index] = Type;
}

CPoint2D&
CPoly2D::GetPoint( size_t index )
{
	return mPoints[index];
}

const CPoint2D&
CPoly2D::GetPoint( size_t index ) const
{
	return mPoints[index];
}

void
CPoly2D::SetPoints( const std::vector<CPoint2D>& Points )
{
	mPoints = Points;
}

void
CPoly2D::InsertPoint( int index, const CPoint2D& p )
{
	mPoints.insert( mPoints.begin() + index, p );
}

void AddUnique( const CPoint2D& p,
				std::vector<CPoint2D>& Vector )
{
	auto iter =std::find_if( Vector.begin(), Vector.end(),
					  [p]( const CPoint2D& l )
	{
		return (l - p).Length() < EPSILON;
	}
	);

	if ( iter == Vector.end() )
		Vector.push_back( p );
}

void
CPoly2D::GetOrderedIntersections( const CLine2D& line,
								  std::vector<CPoint2D>& Intersections ) const
{
	bool Result = true;
	std::vector<CPoint2D> Points;

	CLine2D pline;
	for ( int i = 0; i < (int)mPoints.size() - 1; ++i )
	{
		pline.Set( mPoints[i], mPoints[i + 1], false );
		Intersect( pline, line, Points );

		for ( int j = 0; j<(int)Points.size(); ++j )
			AddUnique( Points[j], Intersections );
	}

	pline.Set( mPoints[mPoints.size() - 1], mPoints[0], false );
	Intersect( pline, line, Points );

	for ( int j = 0; j<(int)Points.size(); ++j )
		AddUnique( Points[j], Intersections );

	AddUnique( line.P1(), Intersections );
	AddUnique( line.P2(), Intersections );

	std::sort( Intersections.begin(), Intersections.end(),
			   [line]( const CPoint2D& l, const CPoint2D& r )
	{
		return (line.P1() - l).Length() < (line.P1() - r).Length();
	}
	);
}

std::vector<CPoint2D>
CPoly2D::IntersectWith( const CLine2D& line )
{
	std::vector<CPoint2D> result;
	std::vector<CPoint2D> Points;

	CLine2D pline;
	for ( int i = 0; i<(int)mPoints.size() - 1; ++i )
	{
		pline.Set( mPoints[i], mPoints[i + 1], false );
		Intersect( pline, line, Points );

		for ( int j = 0; j<(int)Points.size(); ++j )
		{
			CPoint2D& p = Points[j];
			if ( p.DistanceTo( mPoints[i] ) < EPSILON )
				mPoints[i].mark = 1;
			else if ( p.DistanceTo( mPoints[i + 1] ) < EPSILON )
				mPoints[i + 1].mark = 1;
			else
			{
				InsertPoint( i + 1, p );
				mPoints[i + 1].mark = 1;
			}

			std::vector<CPoint2D>::iterator iter =
				std::find_if( result.begin(), result.end(),
							  [p]( const CPoint2D& l )
			{
				return (l - p).Length() < EPSILON;
			}
			);
			if ( iter == result.end() )
				result.push_back( p );
		}
	}

	if ( mPoints.size() )
	{
		pline.Set( mPoints[mPoints.size() - 1], mPoints[0], false );
		Intersect( pline, line, Points );

		for ( int j = 0; j < (int)Points.size(); ++j )
		{
			CPoint2D& p = Points[j];
			if ( p.DistanceTo( mPoints[mPoints.size() - 1] ) < EPSILON )
				mPoints[mPoints.size() - 1].mark = 1;
			else if ( p.DistanceTo( mPoints[0] ) < EPSILON )
				mPoints[0].mark = 1;
			else
			{
				mPoints.push_back( p );
				mPoints.back().mark = 1;
			}

			std::vector<CPoint2D>::iterator iter =
				std::find_if( result.begin(), result.end(),
							  [p]( const CPoint2D& l )
			{
				return (l - p).Length() < EPSILON;
			}
			);
			if ( iter == result.end() )
				result.push_back( p );
		}
	}

	std::sort( result.begin(), result.end(),
			   [line]( const CPoint2D& l, const CPoint2D& r )
	{
		return (line.P1() - l).Length() < (line.P1() - r).Length();
	}
	);

	return result;
}

size_t
CPoly2D::GetNumPoints() const
{
	return mPoints.size();
}

void 
CPoly2D::Reverse()
{
	std::reverse( mPoints.begin(), mPoints.end() );
}

double
CPoly2D::Area() const
{
	CPoint2D p1, p2;

	double a = 0;

	for ( int i = 0; i<mPoints.size() - 1; i++ )
	{
		p1 = mPoints[i];
		p2 = mPoints[i + 1];

		a = a + (p1[0] * p2[1] - p2[0] * p1[1]);
	}

	p1 = mPoints[mPoints.size() - 1];
	p2 = mPoints[0];

	a = a + (p1[0] * p2[1] - p2[0] * p1[1]);

	a *= 0.5;

	return a;
}

CPoint2D
CPoly2D::Centroid() const
{
	CPoint2D result( 0, 0 );
	CPoint2D p1, p2;

	double a = 0;

	for ( int i = 0; i<mPoints.size() - 1; i++ )
	{
		p1 = mPoints[i];
		p2 = mPoints[i + 1];

		a = a + (p1[0] * p2[1] - p2[0] * p1[1]);
	}

	p1 = mPoints[mPoints.size() - 1];
	p2 = mPoints[0];

	a = a + (p1[0] * p2[1] - p2[0] * p1[1]);

	a *= 0.5;

	for ( int i = 0; i<mPoints.size() - 1; i++ )
	{
		p1 = mPoints[i];
		p2 = mPoints[i + 1];

		result[0] = result[0] + (p1[0] + p2[0])*(p1[0] * p2[1] - p2[0] * p1[1]);
		result[1] = result[1] + (p1[1] + p2[1])*(p1[0] * p2[1] - p2[0] * p1[1]);
	}

	p1 = mPoints[mPoints.size() - 1];
	p2 = mPoints[0];

	result[0] = result[0] + (p1[0] + p2[0])*(p1[0] * p2[1] - p2[0] * p1[1]);
	result[1] = result[1] + (p1[1] + p2[1])*(p1[0] * p2[1] - p2[0] * p1[1]);

	result[0] = result[0] / (6 * a);
	result[1] = result[1] / (6 * a);

	return result;
}

std::vector<CPoly2D>
CPoly2D::ClipLeft( const CLine2D& Line ) const
{
	CPoly2D poly = *this;
	std::vector<CPoint2D> points = poly.IntersectWith( Line );

	std::vector<CPoly2D> polies;
	CPoint2D LastPoint;

	if ( points.size() > 1 )
	{
		do
		{
			CPoint2D p_begin = points[0];
			points.erase( points.begin() );

			int sz = static_cast<int>(poly.GetNumPoints());
			int poly_index = -1;
			int poly_start = -1;

			for ( int i = 0; i < sz; ++i )
			{
				if ( p_begin.DistanceTo( poly.GetPoint( i ) ) < EPSILON )
				{
					poly_index = i;
					break;
				}
			}

			assert( poly_index != -1 );

			poly_start = poly_index;
			poly.GetPoint( poly_index ).mark = 0;

			CPoly2D new_poly;
			new_poly.AddPoint( p_begin );

			++poly_index;
			if ( poly_index >= sz )
				poly_index = 0;
			do
			{
				CPoint2D& p = poly.GetPoint( poly_index );

				if ( p.mark == 1 )
				{
					int point_index = -1;
					for ( int i = 0; i < points.size(); ++i )
					{
						if ( points[i].DistanceTo( p ) < EPSILON )
						{
							point_index = i;
							break;
						}
					}

					assert( point_index != -1 );
					p.mark = 0;

					new_poly.AddPoint( p );

					if ( point_index == 0 )
					{
						LastPoint = points.front();
						points.erase( points.begin() );
						break;
					}

					--point_index;

					poly_index = -1;

					for ( int i = 0; i < sz; ++i )
					{
						if ( points[point_index].DistanceTo( poly.GetPoint( i ) ) < EPSILON )
						{
							poly_index = i;
							break;
						}
					}

					assert( poly_index != -1 );

					poly.GetPoint( poly_index ).mark = 0;
					new_poly.AddPoint( points[point_index] );
					points.erase( points.begin() + point_index );
					points.erase( points.begin() + point_index );
				}
				else
					new_poly.AddPoint( p );

				++poly_index;

				if ( poly_index >= sz )
					poly_index = 0;
			} while ( poly_index != poly_start );

			if ( new_poly.GetNumPoints() > 2 )
				polies.push_back( new_poly );
			else if ( points.size() )
				points.insert( points.begin(), LastPoint );

		} while ( points.size() > 0 );
	}
	else
	{
		int sz = static_cast<int>(poly.GetNumPoints());
		CVector2D vl = Line.P2() - Line.P1();
		vl.Normalize();

		bool zero = true;

		for ( int i = 0; i < sz; ++i )
		{
			CVector2D v = poly.GetPoint( i ) - Line.P1();
			double d = vl*v;
			if ( d<0 && abs( d )>EPSILON )
			{
				zero = false;
				break;
			}
		}

		if ( !zero )
		{
			CPoly2D new_poly;
			for ( int i = 0; i < sz; ++i )
				new_poly.AddPoint( poly.GetPoint( i ) );
			polies.push_back( new_poly );
		}
	}

	for ( size_t i = 0; i < polies.size(); ++i )
		polies[i].MakeAntiClockWise();

	return polies;
}

void
CPoly2D::MakeAntiClockWise()
{
	double A = Area();

	if ( A < 0 )
		Reverse();
}

std::vector<CPoly2D>
CPoly2D::ClipRight( const CLine2D& Line ) const
{
	CPoly2D poly = *this;
	poly.Reverse();
	std::vector<CPoint2D> points = poly.IntersectWith( Line );

	std::vector<CPoly2D> polies;
	CPoint2D LastPoint;

	if ( points.size() > 1 )
	{
		do
		{
			CPoint2D p_begin = points[0];
			points.erase( points.begin() );

			/*bool skip = false;
			if ( points.size() )
			{
				CPoint2D mid = (p_begin + points[0]) / 2;
				if
			}*/

			int sz = static_cast<int>(poly.GetNumPoints());
			int poly_index = -1;
			int poly_start = -1;

			for ( int i = 0; i < sz; ++i )
			{
				if ( p_begin.DistanceTo( poly.GetPoint( i ) ) < EPSILON )
				{
					poly_index = i;
					break;
				}
			}

			assert( poly_index != -1 );

			poly_start = poly_index;
			poly.GetPoint( poly_index ).mark = 0;

			CPoly2D new_poly;
			new_poly.AddPoint( p_begin );

			++poly_index;
			if ( poly_index >= sz )
				poly_index = 0;
			do
			{
				CPoint2D& p = poly.GetPoint( poly_index );

				if ( p.mark == 1 )
				{
					int point_index = -1;
					for ( int i = 0; i < points.size(); ++i )
					{
						if ( points[i].DistanceTo( p ) < EPSILON )
						{
							point_index = i;
							break;
						}
					}

					assert( point_index != -1 );
					p.mark = 0;

					new_poly.AddPoint( p );

					if ( point_index == 0 )
					{
						LastPoint = points.front();
						points.erase( points.begin() );
						break;
					}

					--point_index;

					poly_index = -1;

					for ( int i = 0; i < sz; ++i )
					{
						if ( points[point_index].DistanceTo( poly.GetPoint( i ) ) < EPSILON )
						{
							poly_index = i;
							break;
						}
					}

					assert( poly_index != -1 );

					poly.GetPoint( poly_index ).mark = 0;
					new_poly.AddPoint( points[point_index] );
					points.erase( points.begin() + point_index );
					points.erase( points.begin() + point_index );
				}
				else
					new_poly.AddPoint( p );

				++poly_index;

				if ( poly_index >= sz )
					poly_index = 0;
			} while ( poly_index != poly_start );

			if ( new_poly.GetNumPoints() > 2 )
				polies.push_back( new_poly );
			else if( points.size() )
				points.insert( points.begin(), LastPoint );

		} while ( points.size() > 0 );
	}
	else
	{
		int sz = static_cast<int>(poly.GetNumPoints());
		CVector2D vl = Line.P2() - Line.P1();
		vl.Normalize();

		bool zero = true;

		for ( int i = 0; i < sz; ++i )
		{
			CVector2D v = poly.GetPoint( i ) - Line.P1();
			double d = vl*v;
			if ( d>0 && abs( d )>EPSILON )
			{
				zero = false;
				break;
			}
		}

		if ( !zero )
		{
			CPoly2D new_poly;
			for ( int i = 0; i < sz; ++i )
				new_poly.AddPoint( poly.GetPoint( i ) );
			polies.push_back( new_poly );
		}
	}

	for ( size_t i = 0; i < polies.size(); ++i )
		polies[i].MakeAntiClockWise();
	return polies;
}

bool 
CPoly2D::PointInPoly( const CPoint2D& p ) const
{
	int polyCorners = static_cast<int>(mPoints.size());

	int   i, j = polyCorners - 1;
	bool  oddNodes = false;

	for ( i = 0; i<polyCorners; i++ )
	{
		if ( (mPoints[i].y< p.y && mPoints[j].y >= p.y
			   || mPoints[j].y< p.y && mPoints[i].y >= p.y)
			 && (mPoints[i].x <= p.x || mPoints[j].x <= p.x) )
		{
			oddNodes ^= (mPoints[i].x + (p.y - mPoints[i].y) / (mPoints[j].y - mPoints[i].y)*(mPoints[j].x - mPoints[i].x)<p.x);
		}
		j = i;
	}

	return oddNodes;
}

int
CPoly2D::PointOnPoly( const CPoint2D& p ) const
{
	CLine2D Line;
	for ( size_t i = 1; i < mPoints.size(); ++i )
	{

		Line.Set( mPoints[i - 1], mPoints[i] );

		if ( Line.DistanceTo( p ) < EPSILON )
			return static_cast<int>(i)-1;
	}

	if ( mPoints.size() > 1 )
	{
		Line.Set( mPoints[mPoints.size()-1], mPoints[0] );

		if ( Line.DistanceTo( p ) < EPSILON )
			return static_cast<int>(mPoints.size())-1;
	}

	return -1;
}

CPoly2D 
CPoly2D::GetByPoint( const CPoint2D& Point,
					 const std::vector<CPoly2D>& Polies )
{
	CPoly2D Result;
	for ( const auto& Poly : Polies )
	{
		if ( Poly.PointOnPoly( Point )>-1 || Poly.PointInPoly( Point ) )
		{
			Result = Poly;
			break;
		}
	}

	return Result;
}
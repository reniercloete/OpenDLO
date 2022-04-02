// Point2D.cpp
// Copyright (c) 2022, Renier Cloete
// This program is released under the BSD license. See the file LICENSE.txt for details.

#include "Point2D.h"
#include "Constants.h"

#include <assert.h>
#include <math.h>
#include <cmath>

CPoint2D
CPoint2D::operator+( const CPoint2D& other ) const
{
	CPoint2D result = *this;
	result.x += other.x;
	result.y += other.y;

	return result;
}

CPoint2D
CPoint2D::operator-() const
{
	CPoint2D result = *this;
	result.x = -x;
	result.y = -y;

	return result;
}

CPoint2D
CPoint2D::operator-( const CPoint2D& other ) const
{
	CPoint2D result = *this;
	result.x -= other.x;
	result.y -= other.y;

	return result;
}

double
CPoint2D::operator*( const CPoint2D& other ) const
{
	double u1 = x; double u2 = y;
	double v1 = other.x; double v2 = other.y;

	return  u1*v2 - u2*v1;
}

CPoint2D
CPoint2D::operator*( const double& other ) const
{
	CPoint2D result = *this;

	result.x *= other;
	result.y *= other;

	return result;
}

bool
CPoint2D::operator==( const CPoint2D& other ) const
{
	return 
		abs( x - other.x ) < EPSILON &&
		abs( y - other.y ) < EPSILON;
}

CPoint2D 
CPoint2D::operator/( const double& other ) const
{
	CPoint2D result = *this;

	result.x /= other;
	result.y /= other;

	return result;
}

double
CPoint2D::Dot( const CPoint2D& other ) const
{
	double u1 = x; double u2 = y;
	double v1 = other.x; double v2 = other.y;

	return u1*v1 + u2*v2;
}

double
CPoint2D::Length()
{
	return sqrt( x*x + y*y );
}

double
CPoint2D::LengthSquared()
{
	return x*x + y*y;
}

double
CPoint2D::DistanceTo( const CPoint2D& other )
{
	return operator-( other ).Length();
}

void
CPoint2D::Normalize()
{
	double l = Length();
	if ( l>EPSILON )
	{
		x = x / l;
		y = y / l;
	}
}

double&
CPoint2D::operator[]( int index )
{
	assert( index >= 0 && index <2 );
	if ( index == 0 )
		return x;

	return y;
}

const double&
CPoint2D::operator[]( int index ) const
{
	assert( index >= 0 && index <2 );
	if ( index == 0 )
		return x;

	return y;
}

void
CPoint2D::RotateBy( double angle_deg )
{
	double len = Length();

	double theta = 0;
	if ( abs( x>EPSILON ) )
		theta = atan( y / x );
	else
		theta = PI / 2;

	x = len*cos( angle_deg / 180 * PI + theta );
	y = len*sin( angle_deg / 180 * PI + theta );
}
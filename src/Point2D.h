// Point2D.h
// Copyright (c) 2022, Renier Cloete
// This program is released under the BSD license. See the file LICENSE.txt for details.

#pragma once

class CPoint2D
{
public:
	CPoint2D( double _x, double _y ) :x( _x ), y( _y ), mark( 0 ) {};
	CPoint2D() :x( 0 ), y( 0 ), mark( 0 ) {};
	virtual ~CPoint2D() {};

	double x;
	double y;
	int mark;

	CPoint2D operator+( const CPoint2D& other ) const;
	CPoint2D operator-( const CPoint2D& other ) const;
	CPoint2D operator-() const;
	double operator*( const CPoint2D& other ) const;
	CPoint2D operator*( const double& other ) const;
	CPoint2D operator/( const double& other ) const;
	bool operator==( const CPoint2D& other ) const;
	double& operator[]( int index );
	const double& operator[]( int index ) const;

	double Dot( const CPoint2D& other ) const;
	double DistanceTo( const CPoint2D& other );

	void RotateBy( double angle_deg );

	double Length();
	double LengthSquared();
	void Normalize();
};


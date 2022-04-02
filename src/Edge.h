// Edge.h
// Copyright (c) 2022, Renier Cloete
// This program is released under the BSD license. See the file LICENSE.txt for details.

#pragma once

#include "Line2D.h"
#include "Poly2D.h"
#include "Node.h"
#include "Enums.h"

#include <vector>
#include <array>

class CEdge
{
public:

	CEdge( size_t aN1,
		   size_t aN2,
		   eEdgeType aType,
		   double aLength,
		   double aMpPosx,
		   double aMpNegx,
		   double aMpPosy,
		   double aMpNegy );

	void GetCompatibilityMatrix( std::vector<std::vector<double>>& Matrix,
								 bool ApplBoundaryConditions );
	void GetUDLLoadVector( std::array<double, 3>& UDLVector,
						   const CPoly2D& Outline );

	size_t N1, N2;
	eEdgeType Type;
	double Length;
	double MpPosx, MpNegx, MpPosy, MpNegy;
	double MpPos, MpNeg;
	CLine2D Line;
	double YieldRatio;
	bool Added;
	bool Removeable;
	bool Delete;

	int DOF();
	size_t ID() { return mID; };

	static void SetNodes( std::vector<CNode>* Nodes );
	const CNode& GetNode( size_t Index ) { return mNodes->at( Index ); }
	const CPoint2D& GetPoint( size_t Index );

private:
	std::array<double, 3> mUDLVector;
	bool UDLVectorCalculated;
	size_t mID;
	static int Counter;

	static std::vector<CNode>* mNodes;
};
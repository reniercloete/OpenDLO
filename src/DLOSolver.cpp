// DLOSolver.cpp
// Copyright (c) 2022, Renier Cloete
// This program is released under the BSD license. See the file LICENSE.txt for details.

#include "DLOSolver.h"

#include "Domain.h"

#include <algorithm>

CDLOSolver::~CDLOSolver()
{
}

void 
CDLOSolver::fCalculateNodalForces( std::map<size_t, std::array<double, 3>>& Forces,
								   double* rowDual )
{
	std::vector<CNode>& Nodes = mDomain->mNodes;

	std::array<double, 3> Values = { 0,0,0 };
	for ( size_t i = 0; i < Nodes.size(); ++i )
		Forces[i + 1] = Values;

	for ( size_t i = 0; i < Nodes.size(); ++i )
	{
		Forces[i + 1][0] = rowDual[3 * (i + 1) - 3];
		Forces[i + 1][1] = rowDual[3 * (i + 1) - 2];
		Forces[i + 1][2] = rowDual[3 * (i + 1) - 1];
	}
}

bool 
CDLOSolver::fNewViolatedEdges( double Lambda,
							   double* rowDual )
{
	const double kYieldZero = 1e-6;
	std::vector<CEdge*>& Edges = mDomain->mEdges;

	double LiveLoad = mDomain->mLiveLoad;
	double DeadLoad = mDomain->mDeadLoad;

	std::map<size_t, std::array<double, 3>> Forces;
	fCalculateNodalForces( Forces, rowDual );

	std::vector<double> Row = { 0,0,0 };
	std::vector<std::vector<double>> B;
	B.resize( 6, Row );
	std::array<double, 3> EdgefL, EdgefD;

	bool Result = false;

	std::vector<CEdge*> NewEdges, OldEdges;

	for ( auto Edge : Edges )
	{
		if ( Edge->Removeable )
		{
			Edge->GetCompatibilityMatrix( B, false );
			Edge->GetUDLLoadVector( EdgefL, mDomain->mPoly );
			EdgefD = EdgefL;

			double Mn = 0;

			for ( size_t Row = 0; Row < 3; ++Row )
				Mn += B[Row][0] * Forces[Edge->N1][Row];

			for ( size_t Row = 3; Row < 6; ++Row )
				Mn += B[Row][0] * Forces[Edge->N2][Row - 3];

			Mn += Lambda*EdgefL[0] * LiveLoad;
			Mn += EdgefD[0] * DeadLoad;

			if ( Mn < 0 )
				Edge->YieldRatio = abs( Mn / (Edge->MpNeg*Edge->Length) );
			else
				Edge->YieldRatio = abs( Mn / (Edge->MpPos*Edge->Length) );

			if ( Edge->YieldRatio - 1.0 > kYieldZero )
			{
				if ( !Edge->Added )
					NewEdges.push_back( Edge );
			}
		}
	}

	Result = NewEdges.size();

	if ( Result )
	{
		std::sort( NewEdges.begin(), NewEdges.end(),
				   []( const CEdge* l, const CEdge* r ) -> bool
		{
			return l->YieldRatio > r->YieldRatio;
		} );

		size_t Fraction = 5;

		size_t NumberToUse = fGetEdgeCount()*Fraction/100;
		if ( NumberToUse == 0 )
			NumberToUse = NewEdges.size();
		if ( NumberToUse > NewEdges.size() )
			NumberToUse = NewEdges.size();

		int Count = 0;
		for ( auto Edge : NewEdges )
		{
			if ( Count < NumberToUse )
			{
				Edge->Added = true;
				++Count;
			}
			else
				break;
		}
	}

	return Result;
}

double 
CDLOSolver::Solve( CDomain* Domain )
{
	double Result = 0;
	double Old;
	int IterationCount = 1;
	int SameCount = 0;

	mDomain = Domain;

	double *dualRow = fSolve( Result );
	bool Violations = dualRow ? fNewViolatedEdges( Result, dualRow ) : false;
	delete[] dualRow;

	while ( Violations )
	{
		Old = Result;
		++IterationCount;
		dualRow = fSolve( Result );

		if ( abs( Old - Result ) < 1e-6 )
			++SameCount;
		else
			SameCount = 0;

		if ( SameCount == 10 )
			break;

		Violations = dualRow ? fNewViolatedEdges( Result, dualRow ) : false;
		delete[] dualRow;
	}

	return Result;
}

size_t
CDLOSolver::fGetEdgeDOFCount()
{
	std::vector<CEdge*>& Edges = mDomain->mEdges;

	size_t Result = 0;
	for ( size_t i = 0; i < Edges.size(); ++i )
	{
		if ( Edges[i]->Added )
			Result += Edges[i]->DOF();
	}

	return Result;
}

size_t
CDLOSolver::fGetEdgeVarCount()
{
	std::vector<CEdge*>& Edges = mDomain->mEdges;

	size_t Result = 0;
	for ( size_t i = 0; i < Edges.size(); ++i )
	{
		if ( Edges[i]->Added )
		{
			if ( Edges[i]->Type == eEdgeType::FREE ||
				 Edges[i]->Type == eEdgeType::SYMMETRY )
			{
				Result += 3;
			}
			else
				++Result;
		}
	}

	return Result;
}

size_t
CDLOSolver::fGetYieldingEdges()
{
	std::vector<CEdge*>& Edges = mDomain->mEdges;
	size_t Result = 0;
	for ( size_t i = 0; i < Edges.size(); ++i )
	{
		if ( Edges[i]->Added )
		{
			if ( Edges[i]->Type != eEdgeType::FREE &&
				 Edges[i]->Type != eEdgeType::SIMPLE_ANCHORED )
			{
				++Result;
			}
		}
	}

	return Result;
}

size_t
CDLOSolver::fGetEdgeCount()
{
	std::vector<CEdge*>& Edges = mDomain->mEdges;
	size_t Result = 0;
	for ( size_t i = 0; i < Edges.size(); ++i )
	{
		if ( Edges[i]->Added )
			++Result;
	}

	return Result;
}

std::vector<double> 
CDLOSolver::GetEdgeData()
{
	std::vector<CEdge*>& Edges = mDomain->mEdges;
	std::vector<CNode>& Nodes = mDomain->mNodes;

	int YieldCount = 1, RowCount = 1;
	std::vector<double> Result;

	for ( size_t i = 0; i < Edges.size() && mResultArray; ++i )
	{
		if ( Edges[i]->Added )
		{
			if ( Edges[i]->Type == eEdgeType::FREE )
			{
				Result.push_back( mResultArray[RowCount - 1] ); //phin
				++RowCount;
				Result.push_back( mResultArray[RowCount - 1] ); //phit
				++RowCount;
				Result.push_back( mResultArray[RowCount - 1] ); //d
				++RowCount;
				Result.push_back( 0 );

				auto& p1 = Nodes[Edges[i]->N1 - 1].Point;
				auto& p2 = Nodes[Edges[i]->N2 - 1].Point;

				Result.push_back( p1.x );
				Result.push_back( p1.y );
				Result.push_back( p2.x );
				Result.push_back( p2.y );
			}
		    else if (Edges[i]->Type == eEdgeType::SIMPLE_ANCHORED )
			{
				Result.push_back( mResultArray[RowCount - 1] ); //phin
				++RowCount;
				Result.push_back( 0 ); //phit
				Result.push_back( 0 ); //d
				Result.push_back( 0 );

				auto& p1 = Nodes[Edges[i]->N1 - 1].Point;
				auto& p2 = Nodes[Edges[i]->N2 - 1].Point;

				Result.push_back( p1.x );
				Result.push_back( p1.y );
				Result.push_back( p2.x );
				Result.push_back( p2.y );
			}
			else
			{
				double phin = mResultArray[RowCount - 1]; //phin
				++RowCount;

				size_t i1 = mNumDisp + 2 * (YieldCount)-2;
				double pm1 = mResultArray[i1];
				double pm2 = mResultArray[i1 + 1];

				if ( pm1 > 1e-3 || pm2 > 1e-3 )
				{
					Result.push_back( phin ); //phin
					Result.push_back( 0 ); //phit
					Result.push_back( 0 ); //d
					if ( pm1 > 1e-3 )
						Result.push_back( pm1 );
					else
						Result.push_back( pm2 );

					auto& p1 = Nodes[Edges[i]->N1 - 1].Point;
					auto& p2 = Nodes[Edges[i]->N2 - 1].Point;

					Result.push_back( p1.x );
					Result.push_back( p1.y );
					Result.push_back( p2.x );
					Result.push_back( p2.y );
				}

				++YieldCount;
			}
		}
	}

	return Result;
}
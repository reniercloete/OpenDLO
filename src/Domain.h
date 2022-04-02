// Domain.h
// Copyright (c) 2022, Renier Cloete
// This program is released under the BSD license. See the file LICENSE.txt for details.

#pragma once

#include "Point2D.h"
#include "Poly2D.h"
#include "Line2D.h"

#include "Node.h"
#include "Edge.h"

#include "Enums.h"

#include <map>
#include <string>

class CDomain
{
	friend class CDLOSolver;
	friend class CMosekDLOSolver;
	friend class CCoinDLOSolver;
public:
	CDomain();
	virtual ~CDomain();

	void AddBoundaryPoint( const CPoint2D& Point,
						   eEdgeType Type );
	void AddOpeningPoint( size_t Index,
						  const CPoint2D& Point );
	void AddSupport( const CPoint2D& P1,
					 const CPoint2D& P2,
					 eEdgeType Type );

	void Discretize( double Size );
	void BuildEdges();

	void SetLoads( double Live,
				   double Dead )
	{
		mLiveLoad = Live;
		mDeadLoad = Dead;
	}

	void SetYieldMoments( double MpPosx,
						  double MpNegx,
						  double MpPosy,
						  double MpNegy )
	{
		mMpPosx = mMpPosx;
		mMpNegx = mMpNegx;
		mMpPosy = mMpPosy;
		mMpNegy = mMpNegy;
	}

	void Save( const std::string& File );
	void Load( const std::string& File );

	const std::vector<CEdge*>& GetEdges() { return mEdges; };
	const std::vector<CEdge>& GetBoundaryEdges() { return mBoundaryEdges; };
	const std::vector<CPoint2D>& GetBoundaryPoints() { return mPoly.GetPoints(); }

protected:
	struct sSupport
	{
		CLine2D mLine;
		eEdgeType Type;
	};
	double mLiveLoad;
	double mDeadLoad;
	
	CPoly2D mPoly;
	std::vector<CPoly2D> mOpenings;
	std::vector<sSupport> mSupports;

	std::vector<CNode> mNodes;
	std::vector<CEdge> mBoundaryEdges;
	std::vector<CEdge> mMeshEdges;
	std::vector<CEdge> mAdditionalEdges;
	std::vector<CEdge*> mEdges;

	double mMpPosx;
	double mMpNegx;
	double mMpPosy;
	double mMpNegy;

	std::map<size_t, std::vector<size_t>> mNodeMap;

	void fTesselate( double Size );
	void fCreateNodes( double Size );

	size_t fAddNode( const CPoint2D& Point );
	void fAddEdge( size_t N1, size_t N2,
				   eEdgeType Type,
				   std::vector<CEdge>& Edges );
	bool fEdgeAdded( size_t N1,
					 size_t N2 );
	void fMoveToEndOfPartition( std::vector<CEdge*>&Edges,
								size_t& End );
	void fRemoveOverlappedEdges( std::vector<CEdge*>& Edges );
	void fOverlapInternal( std::vector<CEdge*>& Edges,
						   size_t Start,
						   size_t End );
	void fRemoveExteriorEdges( std::vector<CEdge*>& Edges,
							   const CPoly2D& Poly );
	void fCalulculateUDLFactors( std::vector<CEdge*>& Edges,
								 size_t Start,
								 size_t End );
	void fCalculateUDL( std::vector<double>& LoadVector,
					   double UDL );
	
};

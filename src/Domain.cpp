// Domain.cpp
// Copyright (c) 2022, Renier Cloete
// This program is released under the BSD license. See the file LICENSE.txt for details.

#include "Domain.h"

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#ifndef REAL
#define REAL double
#endif
#endif /* not SINGLE */

#include "triangle.h"

#include <assert.h>
#include <thread>
#include <fstream>
#include <algorithm>

CDomain::CDomain():
	mMpPosx( 1.0 ),
	mMpNegx( 1.0 ),
	mMpPosy( 1.0 ),
	mMpNegy( 1.0 )
{
	CEdge::SetNodes( &mNodes );
}

CDomain::~CDomain()
{
	CEdge::SetNodes( nullptr );
}

void 
CDomain::AddBoundaryPoint( const CPoint2D& Point,
						   eEdgeType Type )
{
	size_t Index = mPoly.AddPoint( Point, true );

	if ( Index != -1 )
		mPoly.SetEdgeType( Index, Type );
}

void 
CDomain::AddOpeningPoint( size_t Index,
						  const CPoint2D& Point )
{
	if ( Index >= mOpenings.size() )
	{
		for ( size_t i = 0; i < Index - mOpenings.size() + 1; ++i )
			mOpenings.push_back( CPoly2D() );
	}

	mOpenings[Index].AddPoint( Point, true );
}

void 
CDomain::AddSupport( const CPoint2D& P1,
					 const CPoint2D& P2,
					 eEdgeType Type )
{
	mSupports.push_back( { { P1,P2 },Type } );
}

size_t
CDomain::fAddNode( const CPoint2D& Point )
{
	size_t Result = 0;

	bool Add = true;
	double d;
	for ( size_t i = 0; i < mNodes.size(); ++i )
	{
		d = (mNodes[i].Point - Point).LengthSquared();
		if ( d < 1e-20 )
		{
			Add = false;
			Result = mNodes[i].ID;
			break;
		}
	}

	if ( Add )
	{
		mNodes.push_back( { Point, mNodes.size() + 1 } );
		Result = mNodes.back().ID;
	}

	return Result;
}

void
CDomain::fAddEdge( size_t N1, size_t N2,
				   eEdgeType Type,
				   std::vector<CEdge>& Edges )
{
	double d = (mNodes[N1 - 1].Point - mNodes[N2 - 1].Point).Length();
	Edges.push_back(
		{ 
			N1,
			N2,
			Type,
			d,
			mMpPosx,
			mMpNegx,
			mMpPosy,
			mMpNegy
		}
	);
}

void 
CDomain::fTesselate( double Size )
{
	CVector2D v;
	std::vector<CEdge> NewEdges;

	std::vector<CPoint2D> Points = mPoly.GetPoints();
	Points.push_back( Points.front() );

	for ( size_t i = 1; i < Points.size(); ++i )
	{
		const CPoint2D p1 = Points[i - 1];
		const CPoint2D p2 = Points[i];

		eEdgeType Type = mPoly.GetEdgeType( i - 1 );
				
		v = p2 - p1;
		double Length = v.Length();
		v.Normalize();

		int Number = static_cast<int>(floor( Length / (Size / 2) + 0.5 ));
		double Spacing = Length / Number;

		size_t N1 = fAddNode( p1 );

		for ( int j = 0; j < Number; ++j )
		{
			size_t N2 = fAddNode( p1 + v*(j + 1)*Spacing );
			fAddEdge( N1, N2, Type, mBoundaryEdges );
			N1 = N2;
		}
	}
}

void 
CDomain::fCreateNodes( double Size )
{
	struct sSegmentMarker
	{
		int8_t Padding;
		int8_t Type;
		int16_t ID;
	};
	triangulateio Input;

	Input.numberofpointattributes = 0;
	Input.pointmarkerlist = NULL;
	Input.numberofpointattributes = 0;
	Input.pointattributelist = NULL;

	Input.numberofpoints = (int)mNodes.size();
	Input.pointlist = (double*)malloc( Input.numberofpoints * 2 * sizeof( double ) );

	if ( Input.pointlist == nullptr )
		throw new std::runtime_error( "\"PointList\" allocation unsuccessful" );

	Input.pointmarkerlist = (int*)malloc( Input.numberofpoints * sizeof( int ) );

	if ( Input.pointmarkerlist == nullptr )
		throw new std::runtime_error( "\"PointMarkerList\" allocation unsuccessful" );

	for ( unsigned int i = 0; i<mNodes.size(); ++i )
	{
		Input.pointlist[2 * i] = mNodes[i].Point.x;
		Input.pointlist[2 * i + 1] = mNodes[i].Point.y;
		Input.pointmarkerlist[i] = static_cast<int>(mNodes[i].ID);
	}

	Input.numberofregions = 0;
	Input.numberofholes = 0;

	Input.holelist = 0;

	Input.numberofsegments = (int)(mBoundaryEdges.size());

	Input.segmentlist = new int[Input.numberofsegments * 2];
	Input.segmentmarkerlist = (int*)malloc( Input.numberofsegments * sizeof( int ) );

	for ( size_t i = 0; i<mBoundaryEdges.size(); ++i )
	{
		auto& Edge = mBoundaryEdges[i];

		Input.segmentlist[2 * i] = static_cast<int>(Edge.N1) - 1;
		Input.segmentlist[2 * i + 1] = static_cast<int>(Edge.N2) - 1;
		Input.segmentmarkerlist[i] = static_cast<int>(Edge.ID());
	}

	Input.numberofregions = 0;
	Input.regionlist = nullptr;

	struct triangulateio mid, out, vorout, final;

	mid.pointlist = nullptr;				/* Not needed if -N switch used. */
											/* Not needed if -N switch used or number of point attributes is zero: */
	mid.pointattributelist = nullptr;
	mid.pointmarkerlist = nullptr;			/* Not needed if -N or -B switch used. */
	mid.trianglelist = nullptr;				/* Not needed if -E switch used. */
											/* Not needed if -E switch used or number of triangle attributes is zero: */
	mid.triangleattributelist = nullptr;
	mid.neighborlist = nullptr;				/* Needed only if -n switch used. */
											/* Needed only if segments are output (-p or -c) and -P not used: */
	mid.segmentlist = nullptr;				/* Needed only if segments are output (-p or -c) and -P and -B not used: */
	mid.segmentmarkerlist = nullptr;
	mid.edgelist = nullptr;					/* Needed only if -e switch used. */
	mid.edgemarkerlist = nullptr;			/* Needed if -e used and -B not used. */

	vorout.pointlist = nullptr;				/* Needed only if -v switch used. */
											/* Needed only if -v switch used and number of attributes is not zero: */
	vorout.pointattributelist = nullptr;
	vorout.edgelist = nullptr;				/* Needed only if -v switch used. */
	vorout.normlist = nullptr;				/* Needed only if -v switch used. */

											/* Triangulate the points.  Switches are chosen to read and write a  */
											/*   PSLG (p), preserve the convex hull (c), number everything from  */
											/*   zero (z), assign a regional attribute to each element (A), and  */
											/*   produce an edge list (e), a Voronoi diagram (v), and a triangle */
											/*   neighbor list (n).*/

	mid.triangleattributelist = nullptr;
	mid.pointmarkerlist = nullptr;
	mid.edgemarkerlist = nullptr;
	mid.edgelist = nullptr;
	mid.trianglearealist = nullptr;
	mid.neighborlist = nullptr;
	mid.segmentlist = nullptr;
	mid.segmentmarkerlist = nullptr;
	mid.holelist = nullptr;
	mid.regionlist = Input.regionlist;
	mid.normlist = nullptr;

	//-p : PSLG
	//-z : Number from zero
	//-A : Assign regional attribute
	//-Y : Disable steiner points
	//-e : List of edges
	//-v : Voronoi diagram
	//-n : Outputs neighbours
	//-Q : Quiet

	triangulate( "pzAYYevnV", &Input, &mid, &vorout );

	/* Attach area constraints to the triangles in preparation for */
	/*   refining the triangulation.                               */

	double Sz = sqrt( 3 ) / 8 * Size*Size;
	/* Needed only if -r and -a switches used: */
	mid.trianglearealist = (double *)malloc( mid.numberoftriangles * sizeof( double ) );
	for ( int i = 0; i<mid.numberoftriangles; i++ )
	{
		if ( mid.triangleattributelist )
		{
			if ( mid.triangleattributelist[i] != 1 )
				mid.trianglearealist[i] = Sz;//0.4*(Size*Size)*tan( PI / 3 );
			else
				mid.trianglearealist[i] = -Sz;//-0.4*(Size*Size)*tan( PI / 3 );
		}
		else
			mid.trianglearealist[i] = Sz; //0.4*(Size*Size)*tan( PI / 3 );
	}/**/

	 /* Make necessary initializations so that Triangle can return a */
	 /*   triangulation in `out'.                                    */

	out.pointlist = nullptr;				/* Not needed if -N switch used. */
											/* Not needed if -N switch used or number of attributes is zero: */
	out.pointattributelist = nullptr;
	out.trianglelist = nullptr;				/* Not needed if -E switch used. */
											/* Not needed if -E switch used or number of triangle attributes is zero: */
	out.triangleattributelist = nullptr;
	out.pointmarkerlist = nullptr;
	out.edgemarkerlist = nullptr;
	out.edgelist = nullptr;
	out.trianglearealist = nullptr;
	out.neighborlist = nullptr;
	out.segmentlist = nullptr;
	out.segmentmarkerlist = nullptr;
	out.holelist = nullptr;
	out.regionlist = Input.regionlist;
	out.numberofregions = Input.numberofregions;
	out.normlist = nullptr;

	final.pointlist = nullptr;				/* Not needed if -N switch used. */
											/* Not needed if -N switch used or number of attributes is zero: */
	final.pointattributelist = nullptr;
	final.trianglelist = nullptr;			/* Not needed if -E switch used. */
											/* Not needed if -E switch used or number of triangle attributes is zero: */
	final.triangleattributelist = nullptr;
	final.pointmarkerlist = nullptr;
	final.edgemarkerlist = nullptr;
	final.edgelist = nullptr;
	final.trianglearealist = nullptr;
	final.neighborlist = nullptr;
	final.segmentlist = nullptr;
	final.segmentmarkerlist = nullptr;
	final.holelist = nullptr;
	final.regionlist = Input.regionlist;
	final.numberofregions = Input.numberofregions;
	final.normlist = nullptr;

	/* Refine the triangulation according to the attached */
	/*   triangle area constraints.                       */

	triangulate( "pYYrzaneV", &mid, &out, nullptr );
	triangulate( "pYYzaPneqV", &out, &final, nullptr );

	for ( int i = 0; i<final.numberofedges; i++ )
	{
		size_t N1 = static_cast<size_t>(final.edgelist[2 * i]);
		N1 = fAddNode( { final.pointlist[2 * N1], final.pointlist[2 * N1 + 1] } );

		size_t N2 = static_cast<size_t>(final.edgelist[2 * i + 1]);
		N2 = fAddNode( { final.pointlist[2 * N2], final.pointlist[2 * N2 + 1] } );

		if ( final.edgemarkerlist[i] == 0 )
		{
			if ( !fEdgeAdded( N1, N2 ) )
			{
				mNodeMap[N1].push_back( N2 );
				mNodeMap[N2].push_back( N1 );

				fAddEdge( N1, N2, eEdgeType::INTERNAL, mMeshEdges );
				mMeshEdges.back().Removeable = false;
			}
		}
	}
}

bool 
CDomain::fEdgeAdded( size_t N1,
					 size_t N2 )
{
	bool Result = false;

	auto IVal = mNodeMap.find( N1 );
	if ( IVal != mNodeMap.end() )
	{
		for ( auto k : IVal->second )
		{
			if ( N2 == k )
			{
				Result = true;
				break;
			}
		}
	}

	return Result;
}

void 
CDomain::Discretize( double Size )
{
	mNodes.clear();
	mBoundaryEdges.clear();
	mMeshEdges.clear();
	mAdditionalEdges.clear();
	mEdges.clear();

	fTesselate( Size );
	fCreateNodes( Size );
}

void 
CDomain::BuildEdges()
{
	for ( size_t k = 0; k < mBoundaryEdges.size(); ++k )
	{
		bool Skip = false;

		mBoundaryEdges[k].Removeable = false;

		if ( !fEdgeAdded( mBoundaryEdges[k].N1, mBoundaryEdges[k].N2 ) )
		{
			mNodeMap[mBoundaryEdges[k].N1].push_back( mBoundaryEdges[k].N2 );
			mNodeMap[mBoundaryEdges[k].N2].push_back( mBoundaryEdges[k].N1 );
		}
		else
			assert( 0 );
	}

	for ( int i = 0; i < mNodes.size(); ++i )
	{
		for ( int j = i + 1; j < mNodes.size(); ++j )
		{
			if ( !fEdgeAdded( i + 1, j + 1 ) )
			{
				mNodeMap[i + 1].push_back( j + 1 );
				mNodeMap[j + 1].push_back( i + 1 );

				fAddEdge( mNodes[i].ID, mNodes[j].ID, eEdgeType::INTERNAL, mAdditionalEdges );
				mAdditionalEdges.back().Removeable = true;
			}
		}
	}/**/

	mEdges.reserve( mMeshEdges.size() + mBoundaryEdges.size() + mAdditionalEdges.size() );
	for ( size_t i = 0; i < mBoundaryEdges.size(); ++i )
	{
		mBoundaryEdges[i].Added = true;
		mEdges.push_back( &mBoundaryEdges[i] );
	}
	for ( size_t i = 0; i < mMeshEdges.size(); ++i )
	{
		mMeshEdges[i].Added = true;
		mEdges.push_back( &mMeshEdges[i] );
	}
	for ( size_t i = 0; i < mAdditionalEdges.size(); ++i )
	{
		mEdges.push_back( &mAdditionalEdges[i] );
	}

	fRemoveOverlappedEdges( mEdges );
	fRemoveExteriorEdges( mEdges, mPoly );


	size_t Count = mEdges.size() / 4;

	size_t Start1 = 0;
	size_t End1 = Count;

	size_t Start2 = End1;
	size_t End2 = Start2 + Count;

	size_t Start3 = End2;
	size_t End3 = Start3 + Count;

	size_t Start4 = End3;
	size_t End4 = mEdges.size();

	std::thread T1( [&]() { fCalulculateUDLFactors( mEdges, Start1, End1 ); } );
	std::thread T2( [&]() { fCalulculateUDLFactors( mEdges, Start2, End2 ); } );
	std::thread T3( [&]() { fCalulculateUDLFactors( mEdges, Start3, End3 ); } );
	std::thread T4( [&]() { fCalulculateUDLFactors( mEdges, Start4, End4 ); } );

	T1.join();
	T2.join();
	T3.join();
	T4.join();
}

void
CDomain::fCalulculateUDLFactors( std::vector<CEdge*>& Edges,
								 size_t Start,
								 size_t End )
{
	std::array<double, 3> EdgefL;

	for ( size_t i = Start; i < End; ++i )
		Edges[i]->GetUDLLoadVector( EdgefL, mPoly );
}

void 
CDomain::fMoveToEndOfPartition( std::vector<CEdge*>& Edges,
								size_t& End )
{
	size_t i = End;

	while ( i < Edges.size() &&
			abs( Edges[i - 1]->Line.Slope() - Edges[i]->Line.Slope() ) < 1e-9 )
	{
		++i;
	}

	End = i;
}

void
CDomain::fOverlapInternal( std::vector<CEdge*>& Edges,
						   size_t Start,
						   size_t End )
{
	std::vector<CPoint2D> Intersections;

	for ( size_t i = Start; i < End; ++i )
	{
		CEdge* EdgeI = Edges[i];
		const CLine2D& LineI = EdgeI->Line;
		const CPoint2D& MaxI = LineI.Max();
		const CPoint2D& MinI = LineI.Min();

		if ( !EdgeI->Delete )
		{
			for ( size_t j = i + 1; j < End; ++j )
			{
				CEdge* EdgeJ = Edges[j];
				const CLine2D& LineJ = EdgeJ->Line;
				const CPoint2D& MaxJ = LineJ.Max();
				const CPoint2D& MinJ = LineJ.Min();

				if ( abs( LineI.Slope() - LineJ.Slope() ) > 1e-9 )
				{
					break;
				}

				if ( !EdgeJ->Delete )
				{
					bool skip = MaxJ.x < MinI.x;
					skip |= MinJ.x > MaxI.x;
					skip |= MaxJ.y < MinI.y;
					skip |= MinJ.y > MaxI.y;
					skip |= !Colinear( LineI, LineJ );
					if ( !skip )
					{
						Intersect( LineI, LineJ, Intersections );
						if ( Intersections.size() > 1 )  //Overlap
						{
							if ( EdgeI->Length > EdgeJ->Length && !EdgeI->Delete )
							{
								EdgeI->Delete = true;
								break;
							}
							else
							{
								EdgeJ->Delete = true;
							}
						}
					}
				}
			}
		}
	}
}

void 
CDomain::fRemoveOverlappedEdges( std::vector<CEdge*>& Edges )
{
	std::sort( Edges.begin(), Edges.end(), []( const CEdge* l, const CEdge* r ) -> bool
	{
		return l->Line.Slope() < r->Line.Slope();
	} );

	size_t Count = Edges.size() / 4;

	size_t Start1 = 0;
	size_t End1 = Count;

	fMoveToEndOfPartition( Edges, End1 );

	size_t Start2 = End1;
	size_t End2 = Start2 + Count;

	fMoveToEndOfPartition( Edges, End2 );

	size_t Start3 = End2;
	size_t End3 = Start3 + Count;

	fMoveToEndOfPartition( Edges, End3 );

	size_t Start4 = End3;
	size_t End4 = Edges.size();

	std::thread T1( [&]() { fOverlapInternal( Edges, Start1, End1 ); } );
	std::thread T2( [&]() { fOverlapInternal( Edges, Start2, End2 ); } );
	std::thread T3( [&]() { fOverlapInternal( Edges, Start3, End3 ); } );
	std::thread T4( [&]() { fOverlapInternal( Edges, Start4, End4 ); } );

	T1.join();
	T2.join();
	T3.join();
	T4.join();

	Count = 0;
	for ( size_t i = Edges.size(); i >= 1; --i )
	{
		if ( Edges[i - 1]->Delete )
		{
			++Count;
			Edges.erase( Edges.begin() + i - 1 );
		}
	}
}

void 
CDomain::fRemoveExteriorEdges( std::vector<CEdge*>& Edges,
							   const CPoly2D& Poly )
{
	CPoint2D p;
	std::vector<CPoint2D> Intersections;

	for ( size_t i = 0; i < Edges.size(); ++i )
	{
		if ( Edges[i]->Removeable )
		{

			Intersections.clear();
			Poly.GetOrderedIntersections( Edges[i]->Line, Intersections );

			for ( size_t j = 0; j < Intersections.size() - 1; ++j )
			{
				p = (Intersections[j] + Intersections[j + 1]) / 2;
				if ( Poly.PointOnPoly( p ) == -1 && !Poly.PointInPoly( p ) )
				{
					Edges[i]->Delete = true;
					break;
				}
			}
		}
	}

	size_t Count = 0;
	for ( size_t i = Edges.size(); i >= 1; --i )
	{
		if ( Edges[i - 1]->Delete )
		{
			++Count;
			Edges.erase( Edges.begin() + i - 1 );
		}
	}
}

void 
CDomain::fCalculateUDL( std::vector<double>& LoadVector,
						double UDL )
{
	size_t Index = 0;

	std::array<double, 3> UDLVector;
	for ( size_t i = 0; i < mEdges.size(); ++i )
	{
		if ( mEdges[i]->Added )
		{
			mEdges[i]->GetUDLLoadVector( UDLVector, mPoly );

			for ( size_t j = 0; j < mEdges[i]->DOF(); ++j )
			{
				LoadVector[Index] += UDL*UDLVector[j];
				++Index;
			}
		}
	}
}

void
CDomain::Save( const std::string& File )
{
	std::ofstream output;
	output.open( File, std::ofstream::out );

	output << mNodes.size() << " ";
	for ( size_t i = 0; i < mNodes.size(); ++i )
		output <<
		mNodes[i].ID << " " <<
		mNodes[i].Point.x << " " <<
		mNodes[i].Point.y << " ";

	output << mBoundaryEdges.size() << " ";
	for ( size_t i = 0; i < mBoundaryEdges.size(); ++i )
		output <<
		mBoundaryEdges[i].N1 << " " <<
		mBoundaryEdges[i].N2 << " " <<
		(int)mBoundaryEdges[i].Type << " ";

	output <<
		mMpPosx << " " <<
		mMpNegx << " " <<
		mMpPosy << " " <<
		mMpNegy << " ";

	output << 0.5;

	output.close();
}

void
CDomain::Load( const std::string& File )
{
	std::ifstream input;
	input.open( File, std::ofstream::in );

	size_t sz;
	input >> sz;

	CPoint2D Point;
	size_t Edge;

	for ( size_t i = 0; i < sz; ++i )
	{
		CNode node;
		input >> Point.x;
		input >> Point.y;
		input >> Edge;
		AddBoundaryPoint( Point, (eEdgeType)Edge );
	}

	double MpPosx, MpNegx, MpPosy, MpNegy, NodeDensity;
	
	input >> MpPosx;
	input >> MpNegx;
	input >> MpPosy;
	input >> MpNegy;

	input >> NodeDensity;

	SetYieldMoments( MpPosx, MpNegx, MpPosy, MpNegy );
	Discretize( NodeDensity );
	SetLoads( 1, 0 );

	input.close();
}



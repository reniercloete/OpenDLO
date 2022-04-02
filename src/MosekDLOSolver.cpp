// MosekDLOSolver.cpp
// Copyright (c) 2022, Renier Cloete
// This program is released under the BSD license. See the file LICENSE.txt for details.

#include "MosekDLOSolver.h"

#include "Constants.h"
#include "Domain.h"

void
CMosekDLOSolver::fGetColumnSolution()
{
	if ( mResultArray )
	{
		delete[]mResultArray;
		mResultArray = nullptr;
	}

	mResultArray = new double[mNumDisp + mNumYEdges * 2];

	mSize = mNumDisp + mNumYEdges * 2;

	if ( mResultArray )
	{
		MSK_getxx( mCurrentTask,
				   MSK_SOL_BAS,    // Request the basic solution. 
				   mResultArray );
	}
}

double* 
CMosekDLOSolver::fSolve( double& Objective )
{
	double* Result = nullptr;

	mNumDisp = fGetEdgeVarCount();
	mNumYEdges = fGetYieldingEdges();
	mNumDOF = fGetEdgeDOFCount();

	size_t numvar = mNumDisp + mNumYEdges * 2;
	MSKint32t numcon = static_cast<MSKint32t>(mDomain->mNodes.size() * 3 + mNumYEdges + 1);

	mfL.clear();
	mfD.clear();

	mfD.resize( mNumDOF, 0.0 );
	mDomain->fCalculateUDL( mfD, mDomain->mDeadLoad );

	mfL.resize( mNumDOF, 0.0 );
	mDomain->fCalculateUDL( mfL, mDomain->mLiveLoad );
	
	mPtrb.resize( numvar, -1 );
	mPtre.resize( numvar, -1 );
	mSub.clear();
	mVal.clear();

	MSKenv_t env = NULL;

	if ( MSK_makeenv( &env, NULL ) != MSK_RES_OK )
		return nullptr;

	mCurrentTask = fBuildModel( env );
	MSK_putintparam( mCurrentTask, MSK_IPAR_NUM_THREADS, 4 );

	MSKrescodee trmcode;

	/* Run optimizer */
	MSKrescodee	r = MSK_putobjsense( mCurrentTask, MSK_OBJECTIVE_SENSE_MINIMIZE );
	r = MSK_optimizetrm( mCurrentTask, &trmcode );

	/* Print a summary containing information
	about the solution for debugging purposes. */
	MSK_solutionsummary( mCurrentTask, MSK_STREAM_LOG );

	if ( r == MSK_RES_OK )
	{
		MSKsolstae solsta;

		if ( r == MSK_RES_OK )
			r = MSK_getsolsta( mCurrentTask,
							   MSK_SOL_BAS,
							   &solsta );
		switch ( solsta )
		{
		case MSK_SOL_STA_OPTIMAL:
		case MSK_SOL_STA_NEAR_OPTIMAL:
		{
			fGetColumnSolution();

			MSK_getprimalobj( mCurrentTask,
							  MSK_SOL_BAS,    // Request the basic solution. 
							  &Objective );

			Result = new double[numcon];
			if ( Result )
			{
				MSK_gety( mCurrentTask,
						  MSK_SOL_BAS,    // Request the basic solution. 
						  Result );

				/*printf( "Optimal dual solution\n" );
				for ( MSKint32t j = 0; j < 2 * Nodes.size() + GetEdgeCount() + 1; ++j )
				printf( "x[%d]: %e\n", j, dualRow[j] );

				free( dualRow );*/
			}
			else
				r = MSK_RES_ERR_SPACE;

			break;
		}
		case MSK_SOL_STA_DUAL_INFEAS_CER:
		case MSK_SOL_STA_PRIM_INFEAS_CER:
		case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
		case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
			printf( "Primal or dual infeasibility certificate found.\n" );
			break;
		case MSK_SOL_STA_UNKNOWN:
		{
			char symname[MSK_MAX_STR_LEN];
			char desc[MSK_MAX_STR_LEN];

			/* If the solutions status is unknown, print the termination code
			indicating why the optimizer terminated prematurely. */

			MSK_getcodedesc( trmcode,
							 symname,
							 desc );

			printf( "The solution status is unknown.\n" );
			printf( "The optimizer terminitated with code: %s\n", symname );
			break;
		}
		default:
			printf( "Other solution status.\n" );
			break;
		}
	}

	if ( r != MSK_RES_OK )
	{
		/* In case of an error print error code and description. */
		char symname[MSK_MAX_STR_LEN];
		char desc[MSK_MAX_STR_LEN];

		printf( "An error occurred while optimizing.\n" );
		MSK_getcodedesc( r,
						 symname,
						 desc );
		printf( "Error %s - '%s'\n", symname, desc );
	}

	/* Delete the task and the associated data. */
	MSK_deletetask( &mCurrentTask );

	/* Delete the environment and the associated data. */
	MSK_deleteenv( &env );

	return Result;
}

MSKtask_t
CMosekDLOSolver::fBuildModel( MSKenv_t env )
{
	MSKint32t nummult = static_cast<MSKint32t>(mNumYEdges * 2);
	MSKint32t numvar = static_cast<MSKint32t>(mNumDisp) + nummult;
	MSKint32t numcon = static_cast<MSKint32t>(mDomain->mNodes.size() * 3 + mNumYEdges + 1);
	MSKrescodee	r;
	MSKtask_t task;
	/* Create the optimization task. */
	r = MSK_maketask( env, numcon, numvar, &task );

	/* Directs the log task stream to the 'printstr' function. */
	//if ( r == MSK_RES_OK )
	//r = MSK_linkfunctotaskstream( task, MSK_STREAM_LOG, NULL, printstr );

	/* Bounds on constraints. */
	std::vector<MSKboundkeye> bkc( numcon, MSK_BK_FX );
	std::vector<double>		  blc( numcon, 0.0 );
	std::vector<double>		  buc( numcon, 0.0 );

	blc.back() = buc.back() = 1.0;

	/* Bounds on variables. */
	std::vector<MSKboundkeye> bkx( numvar );
	std::vector<double>       blx( numvar );
	std::vector<double>       bux( numvar );

	for ( size_t i = 0; i < mNumDisp; ++i )
	{
		bkx[i] = MSK_BK_FR;
		blx[i] = -MSK_INFINITY;
		bux[i] = +MSK_INFINITY;
	}

	for ( size_t i = mNumDisp; i < mNumDisp + nummult; ++i )
	{
		bkx[i] = MSK_BK_LO;
		blx[i] = 0;
		bux[i] = +MSK_INFINITY;
	}

	/* Append 'numcon' empty constraints.
	The constraints will initially have no bounds. */
	if ( r == MSK_RES_OK )
		r = MSK_appendcons( task, numcon );

	/* Append 'numvar' variables.
	The variables will initially be fixed at zero (x=0). */
	if ( r == MSK_RES_OK )
		r = MSK_appendvars( task, numvar );

	fCalculateCompatibilityMatrix( task );
	fCalculatePlasticMultiplierTerms( task );

	int Row = 0;

	size_t dsz = mfD.size();
	size_t psz = mObjP.size();

	for ( size_t i = 0; i < dsz; ++i )
	{
		MSK_putcj( task, Row, mfD[i] );
		MSK_putvarbound( task,
						 Row,           /* Index of variable.*/
						 bkx[Row],      /* Bound key.*/
						 blx[Row],      /* Numerical value of lower bound.*/
						 bux[Row] );	/* Numerical value of upper bound.*/
		++Row;
	}

	for ( size_t i = 0; i < psz; ++i )
	{
		MSK_putcj( task, Row, mObjP[i] );
		MSK_putvarbound( task,
						 Row,           /* Index of variable.*/
						 bkx[Row],      /* Bound key.*/
						 blx[Row],      /* Numerical value of lower bound.*/
						 bux[Row] );	/* Numerical value of upper bound.*/
		++Row;
	}

	for (MSKint32t i = 0; i < numcon; ++i )
		r = MSK_putconbound( task,
							 i,           /* Index of constraint.*/
							 bkc[i],      /* Bound key.*/
							 blc[i],      /* Numerical value of lower bound.*/
							 buc[i] );     /* Numerical value of upper bound.*/

	return task;
}

void 
CMosekDLOSolver::fCalculateCompatibilityMatrix( MSKtask_t task )
{
	std::vector<CEdge*>& Edges = mDomain->mEdges;

	std::vector<double> RowVector;
	std::vector<std::vector<double>> Matrix;
	RowVector.resize( 3, 0.0 );
	Matrix.resize( 6, RowVector );

	int n1, n2;
	size_t ColIndex = 0;

	int Count = 0;

	for ( auto& Edge : Edges )
	{
		if ( Edge->Added &&
			 Edge->Type != eEdgeType::FREE &&
			 Edge->Type != eEdgeType::SIMPLE_ANCHORED )
			++Count;
	}

	mObjP.resize( 2 * Count );

	int NodeSize = static_cast<int>(mDomain->mNodes.size());
	int NumYEdges = static_cast<int>(mNumYEdges);

	int col1, col2;
	Count = 1;
	int YieldingCount = 1;

	std::array<size_t, 3> an1 = { 0,0,0 };
	std::array<size_t, 3> an2 = { 0,0,0 };

	for ( int i = 0; i < Edges.size(); ++i )
	{
		if ( Edges[i]->Added )
		{
			n1 = static_cast<int>(Edges[i]->N1);
			n2 = static_cast<int>(Edges[i]->N2);

			an1 = { 0,1,2 };
			an2 = { 3,4,5 };

			if ( n1 > n2 )
			{
				std::swap( n1, n2 );
				an1 = { 3,4,5 };
				an2 = { 0,1,2 };
			}

			Edges[i]->GetCompatibilityMatrix( Matrix, true );

			if ( Edges[i]->Type != eEdgeType::FREE &&
				 Edges[i]->Type != eEdgeType::SIMPLE_ANCHORED )
			{
				mObjP[2 * YieldingCount - 2] = Edges[i]->MpPos*Edges[i]->Length;
				mObjP[2 * YieldingCount - 1] = Edges[i]->MpNeg*Edges[i]->Length;

				++YieldingCount;
			}

			for ( size_t j = 0; j < Edges[i]->DOF(); ++j )
			{
				col1 = static_cast<int>( mVal.size() );

				if ( j == 0 )
				{
					if ( abs( Matrix[0][j] ) > 0 )
					{
						mVal.push_back( Matrix[an1[0]][j] );
						mSub.push_back( 3 * n1 - 3 );
					}

					if ( abs( Matrix[1][j] ) > 0 )
					{
						mVal.push_back( Matrix[an1[1]][j] );
						mSub.push_back( 3 * n1 - 2 );
					}

					if ( abs( Matrix[3][j] ) > 0 )
					{
						mVal.push_back( Matrix[an2[0]][j] );
						mSub.push_back( 3 * n2 - 3 );
					}

					if ( abs( Matrix[4][j] ) > 0 )
					{
						mVal.push_back( Matrix[an2[1]][j] );
						mSub.push_back( 3 * n2 - 2 );
					}
				}
				else if ( j == 1 )
				{
					if ( abs( Matrix[0][j] ) > 0 )
					{
						mVal.push_back( Matrix[an1[0]][j] );
						mSub.push_back( 3 * n1 - 3 );
					}

					if ( abs( Matrix[1][j] ) > 0 )
					{
						mVal.push_back( Matrix[an1[1]][j] );
						mSub.push_back( 3 * n1 - 2 );
					}

					if ( abs( Matrix[2][j] ) > 0 )
					{
						mVal.push_back( Matrix[an1[2]][j] );
						mSub.push_back( 3 * n1 - 1 );
					}

					if ( abs( Matrix[3][j] ) > 0 )
					{
						mVal.push_back( Matrix[an2[0]][j] );
						mSub.push_back( 3 * n2 - 3 );
					}

					if ( abs( Matrix[4][j] ) > 0 )
					{
						mVal.push_back( Matrix[an2[1]][j] );
						mSub.push_back( 3 * n2 - 2 );
					}

					if ( abs( Matrix[5][j] ) > 0 )
					{
						mVal.push_back( Matrix[an2[2]][j] );
						mSub.push_back( 3 * n2 - 1 );
					}
				}
				else if ( j == 2 )
				{
					mVal.push_back( Matrix[an1[2]][j] );
					mSub.push_back( 3 * n1 - 1 );

					mVal.push_back( Matrix[an2[2]][j] );
					mSub.push_back( 3 * n2 - 1 );

				}


				if ( Edges[i]->Type != eEdgeType::FREE &&
					 Edges[i]->Type != eEdgeType::SIMPLE_ANCHORED )
				{
					mVal.push_back( -1 );
					mSub.push_back( Count - 1 + NodeSize * 3 );

					++Count;
				}


				if ( abs( mfL[ColIndex] ) > EPSILON )
				{
					mVal.push_back( mfL[ColIndex] );
					mSub.push_back( NodeSize * 3 + NumYEdges );
				}

				col2 = static_cast<int>(mVal.size());

				mPtrb[ColIndex] = col1;
				mPtre[ColIndex] = col2;

				MSK_putacol( task,
							 static_cast<MSKint32t>(ColIndex),							/* Variable (column) index./**/
							 mPtre[ColIndex] - mPtrb[ColIndex],	/* Number of non-zeros in column j./**/
							 &mSub[0] + mPtrb[ColIndex],		/* Pointer to row indexes of column j./**/
							 &mVal[0] + mPtrb[ColIndex] );		/* Pointer to Values of column j./**/

				++ColIndex;
			}

		}
	}
}

void 
CMosekDLOSolver::fCalculatePlasticMultiplierTerms( MSKtask_t task )
{
	int col1, col2;

	size_t ColumnIndex = mNumDisp;
	int Row = 3 * static_cast<int>(mDomain->mNodes.size());
	std::vector<CEdge*>& Edges = mDomain->mEdges;

	int Count = 1;

	for ( auto& Edge : Edges )
	{
		if ( Edge->Added )
		{
			if ( Edge->Type != eEdgeType::FREE &&
				 Edge->Type != eEdgeType::SIMPLE_ANCHORED )
			{
				col1 = static_cast<int>(mVal.size());

				mVal.push_back( 1 );
				mSub.push_back( Count - 1 + Row );

				col2 = static_cast<int>(mVal.size());

				mPtrb[ColumnIndex] = col1;
				mPtre[ColumnIndex] = col2;

				MSK_putacol( task,
							 static_cast<MSKint32t>(ColumnIndex),		// Variable (column) index.
							 mPtre[ColumnIndex] - mPtrb[ColumnIndex],	// Number of non-zeros in column j.
							 &mSub[0] + mPtrb[ColumnIndex],				// Pointer to row indexes of column j.
							 &mVal[0] + mPtrb[ColumnIndex] );			// Pointer to Values of column j./**/

				++ColumnIndex;

				col1 = static_cast<int>(mVal.size());

				mVal.push_back( -1 );
				mSub.push_back( Count - 1 + Row );

				col2 = static_cast<int>(mVal.size());

				mPtrb[ColumnIndex] = col1;
				mPtre[ColumnIndex] = col2;

				MSK_putacol( task,
							 static_cast<MSKint32t>(ColumnIndex),		// Variable (column) index.
							 mPtre[ColumnIndex] - mPtrb[ColumnIndex],	// Number of non-zeros in column j.
							 &mSub[0] + mPtrb[ColumnIndex],				// Pointer to row indexes of column j.
							 &mVal[0] + mPtrb[ColumnIndex] );			// Pointer to Values of column j./**/

				++ColumnIndex;
				++Count;
			}
		}
	}
}
// CoinDLOSolver.cpp
// Copyright (c) 2022, Renier Cloete
// This program is released under the BSD license. See the file LICENSE.txt for details.

#include "CoinDLOSolver.h"

#include "Constants.h"
#include "Domain.h"
#include "ClpCholeskyPardiso.hpp"

CCoinDLOSolver::~CCoinDLOSolver()
{
	delete mModel;
}

void
CCoinDLOSolver::fGetColumnSolution()
{
	if (mResultArray)
	{
		delete[]mResultArray;
		mResultArray = nullptr;
	}

	mResultArray = new double[mNumDisp + mNumYEdges * 2];

	mSize = mNumDisp + mNumYEdges * 2;

	if (mResultArray)
	{
		memcpy( mResultArray, mModel->getColSolution(), sizeof( double ) * mSize );
	}
}

double*
CCoinDLOSolver::fSolve( double& Objective )
{
	double* Result = nullptr;

	mNumDisp = fGetEdgeVarCount();
	mNumYEdges = fGetYieldingEdges();
	mNumDOF = fGetEdgeDOFCount();

	size_t numvar = mNumDisp + mNumYEdges * 2;
	size_t numcon = mDomain->mNodes.size() * 3 + mNumYEdges + 1;

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

	mModel = new ClpInterior();
	ClpCholeskyPardiso* cholesky = new ClpCholeskyPardiso();
	static_cast<ClpInterior*>(mModel)->setCholesky( cholesky );
	mModel->setPrimalTolerance( 1e-8 );
	mModel->setDualTolerance( 1e-8 );
	mModel->setLogLevel( 0 );
	mModel->resize( static_cast<int>(numcon), 0 );

	int Row = 0;

	fCalculateCompatibilityMatrix();
	fCalculatePlasticMultiplierTerms();

	size_t dsz = mfD.size();
	size_t psz = mObjP.size();

	for (size_t i = 0; i < dsz; ++i)
	{
		mModel->setObjectiveCoefficient( Row, mfD[i] );
		mModel->setColumnLower( Row, -COIN_DBL_MAX );
		mModel->setColumnUpper( Row, COIN_DBL_MAX );
		++Row;
	}

	for (size_t i = 0; i < psz; ++i)
	{
		mModel->setObjectiveCoefficient( Row, mObjP[i] );
		mModel->setColumnLower( Row, 0.0 );
		mModel->setColumnUpper( Row, COIN_DBL_MAX );
		++Row;
	}

	for (size_t i = 0; i < numcon - 1; ++i)
		mModel->setRowBounds( static_cast<int>(i), 0, 0 );

	mModel->setRowBounds( static_cast<int>(numcon) - 1, 1, 1 );

	static_cast<ClpInterior*>(mModel)->primalDual();

	Objective = mModel->objectiveValue();

	Result = new double[numcon];
	if (Result)
	{

		const double* row = mModel->dualRowSolution();
		memcpy( Result, row, sizeof( double ) * (numcon) );
	}

	fGetColumnSolution();

	return Result;
}

void
CCoinDLOSolver::fBuildModel()
{
}

void
CCoinDLOSolver::fCalculateCompatibilityMatrix()
{
	std::vector<CEdge*>& Edges = mDomain->mEdges;

	std::vector<double> RowVector;
	std::vector<std::vector<double>> Matrix;
	RowVector.resize( 3, 0.0 );
	Matrix.resize( 6, RowVector );

	size_t n1, n2;
	size_t ColIndex = 0;

	size_t Count = 0;

	for (auto& Edge : Edges)
	{
		if (Edge->Added &&
			 Edge->Type != eEdgeType::FREE &&
			 Edge->Type != eEdgeType::SIMPLE_ANCHORED)
			++Count;
	}

	mObjP.resize( 2 * Count );

	size_t NodeSize = mDomain->mNodes.size();
	size_t NumYEdges = mNumYEdges;

	size_t col1, col2;
	Count = 1;
	int YieldingCount = 1;

	std::array<size_t, 3> an1 = { 0,0,0 };
	std::array<size_t, 3> an2 = { 0,0,0 };

	for (size_t i = 0; i < Edges.size(); ++i)
	{
		if (Edges[i]->Added)
		{
			n1 = Edges[i]->N1;
			n2 = Edges[i]->N2;

			an1 = { 0,1,2 };
			an2 = { 3,4,5 };

			if (n1 > n2)
			{
				std::swap( n1, n2 );
				an1 = { 3,4,5 };
				an2 = { 0,1,2 };
			}

			Edges[i]->GetCompatibilityMatrix( Matrix, true );

			if (Edges[i]->Type != eEdgeType::FREE &&
				 Edges[i]->Type != eEdgeType::SIMPLE_ANCHORED)
			{
				mObjP[2 * YieldingCount - 2] = Edges[i]->MpPos * Edges[i]->Length;
				mObjP[2 * YieldingCount - 1] = Edges[i]->MpNeg * Edges[i]->Length;

				++YieldingCount;
			}

			for (size_t j = 0; j < Edges[i]->DOF(); ++j)
			{
				col1 = mVal.size();

				if (j == 0)
				{
					if (abs( Matrix[0][j] ) > 0)
					{
						mVal.push_back( Matrix[an1[0]][j] );
						mSub.push_back( static_cast<int>(3 * n1) - 3 );
					}

					if (abs( Matrix[1][j] ) > 0)
					{
						mVal.push_back( Matrix[an1[1]][j] );
						mSub.push_back( static_cast<int>(3 * n1) - 2 );
					}

					if (abs( Matrix[3][j] ) > 0)
					{
						mVal.push_back( Matrix[an2[0]][j] );
						mSub.push_back( static_cast<int>(3 * n2) - 3 );
					}

					if (abs( Matrix[4][j] ) > 0)
					{
						mVal.push_back( Matrix[an2[1]][j] );
						mSub.push_back( static_cast<int>(3 * n2) - 2 );
					}
				}
				else if (j == 1)
				{
					if (abs( Matrix[0][j] ) > 0)
					{
						mVal.push_back( Matrix[an1[0]][j] );
						mSub.push_back( static_cast<int>(3 * n1) - 3 );
					}

					if (abs( Matrix[1][j] ) > 0)
					{
						mVal.push_back( Matrix[an1[1]][j] );
						mSub.push_back( static_cast<int>(3 * n1) - 2 );
					}

					if (abs( Matrix[2][j] ) > 0)
					{
						mVal.push_back( Matrix[an1[2]][j] );
						mSub.push_back( static_cast<int>(3 * n1) - 1 );
					}

					if (abs( Matrix[3][j] ) > 0)
					{
						mVal.push_back( Matrix[an2[0]][j] );
						mSub.push_back( static_cast<int>(3 * n2) - 3 );
					}

					if (abs( Matrix[4][j] ) > 0)
					{
						mVal.push_back( Matrix[an2[1]][j] );
						mSub.push_back( static_cast<int>(3 * n2) - 2 );
					}

					if (abs( Matrix[5][j] ) > 0)
					{
						mVal.push_back( Matrix[an2[2]][j] );
						mSub.push_back( static_cast<int>(3 * n2) - 1 );
					}
				}
				else if (j == 2)
				{
					mVal.push_back( Matrix[an1[2]][j] );
					mSub.push_back( static_cast<int>(3 * n1) - 1 );

					mVal.push_back( Matrix[an2[2]][j] );
					mSub.push_back( static_cast<int>(3 * n2) - 1 );

				}


				if (Edges[i]->Type != eEdgeType::FREE &&
					 Edges[i]->Type != eEdgeType::SIMPLE_ANCHORED)
				{
					mVal.push_back( -1.0 );
					mSub.push_back( static_cast<int>(Count + NodeSize * 3) - 1 );

					++Count;
				}


				if (abs( mfL[ColIndex] ) > EPSILON)
				{
					mVal.push_back( mfL[ColIndex] );
					mSub.push_back( static_cast<int>(NodeSize * 3 + NumYEdges) );
				}

				col2 = mVal.size();

				mPtrb[ColIndex] = static_cast<int>(col1);
				mPtre[ColIndex] = static_cast<int>(col2);

				mModel->addColumn( mPtre[ColIndex] - mPtrb[ColIndex],
								   &mSub[0] + mPtrb[ColIndex],
								   &mVal[0] + mPtrb[ColIndex] );

				++ColIndex;
			}

		}
	}
}

void
CCoinDLOSolver::fCalculatePlasticMultiplierTerms()
{
	size_t col1, col2;

	size_t ColumnIndex = mNumDisp;
	size_t Row = 3 * mDomain->mNodes.size();
	std::vector<CEdge*>& Edges = mDomain->mEdges;

	size_t Count = 1;

	for (auto& Edge : Edges)
	{
		if (Edge->Added)
		{
			if (Edge->Type != eEdgeType::FREE &&
				 Edge->Type != eEdgeType::SIMPLE_ANCHORED)
			{
				col1 = mVal.size();

				mVal.push_back( 1 );
				mSub.push_back( static_cast<int>(Count + Row) - 1 );

				col2 = mVal.size();

				mPtrb[ColumnIndex] = static_cast<int>(col1);
				mPtre[ColumnIndex] = static_cast<int>(col2);

				mModel->addColumn( mPtre[ColumnIndex] - mPtrb[ColumnIndex],
								   &mSub[0] + mPtrb[ColumnIndex],
								   &mVal[0] + mPtrb[ColumnIndex] );

				++ColumnIndex;

				col1 = mVal.size();

				mVal.push_back( -1 );
				mSub.push_back( static_cast<int>(Count + Row) - 1 );

				col2 = mVal.size();

				mPtrb[ColumnIndex] = static_cast<int>(col1);
				mPtre[ColumnIndex] = static_cast<int>(col2);

				mModel->addColumn( mPtre[ColumnIndex] - mPtrb[ColumnIndex],
								   &mSub[0] + mPtrb[ColumnIndex],
								   &mVal[0] + mPtrb[ColumnIndex] );

				++ColumnIndex;
				++Count;
			}
		}
	}
}
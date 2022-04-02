// DLOSolver.h
// Copyright (c) 2022, Renier Cloete
// This program is released under the BSD license. See the file LICENSE.txt for details.

#pragma once

#include <vector>
#include <map>

class CDomain;

class CDLOSolver
{
public:
	CDLOSolver():
		mResultArray(nullptr),
		mSize(0)
	{
	};
	virtual ~CDLOSolver();

	virtual double Solve( CDomain* Domain );
	double* GetResultArray( size_t& Size )
	{
		Size = mSize;
		return mResultArray;
	};

	std::vector<double> GetEdgeData();
	
protected:
	CDomain* mDomain;
	double* mResultArray;
	size_t mSize;
	size_t mNumDisp;
	size_t mNumYEdges;
	size_t mNumDOF;

	std::vector<int>	mPtrb;
	std::vector<int>	mPtre;
	std::vector<int>	mSub;
	std::vector<double> mVal;
	std::vector<double> mObjP;
	std::vector<double> mfL;
	std::vector<double> mfD;

	virtual double* fSolve( double& Objective ) = 0;
	virtual void fGetColumnSolution() = 0;
	virtual bool fNewViolatedEdges( double Lambda,
									double* rowDual );
	void fCalculateNodalForces( std::map<size_t, std::array<double, 3>>& Forces,
								double* rowDual );

	size_t fGetEdgeCount();
	size_t fGetEdgeDOFCount();
	size_t fGetEdgeVarCount();
	size_t fGetYieldingEdges();
};
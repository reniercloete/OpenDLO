// CoinDLOSolver.h
// Copyright (c) 2022, Renier Cloete
// This program is released under the BSD license. See the file LICENSE.txt for details.

#pragma once

#include "DLOSolver.h"
#include "ClpInterior.hpp"

class CCoinDLOSolver : public CDLOSolver
{
public:
	CCoinDLOSolver():
		mModel( nullptr )
	{
	};

	virtual ~CCoinDLOSolver();



protected:
	ClpModel* mModel = nullptr;
	double* fSolve( double& Objective ) override;
	void fGetColumnSolution() override;

	void fBuildModel();
	void fCalculateCompatibilityMatrix();
	void fCalculatePlasticMultiplierTerms();

};
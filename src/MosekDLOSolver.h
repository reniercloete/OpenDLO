// MosekDLOSolver.h
// Copyright (c) 2022, Renier Cloete
// This program is released under the BSD license. See the file LICENSE.txt for details.

#pragma once

#include "DLOSolver.h"
#include "mosek.h"

class CMosekDLOSolver : public CDLOSolver
{
public:
	CMosekDLOSolver():
		mCurrentTask(nullptr)
	{
	};

	virtual ~CMosekDLOSolver() {};



protected:
	MSKtask_t mCurrentTask;
	double* fSolve( double& Objective ) override;
	void fGetColumnSolution() override;

	MSKtask_t fBuildModel( MSKenv_t env );
	void fCalculateCompatibilityMatrix( MSKtask_t task );
	void fCalculatePlasticMultiplierTerms( MSKtask_t task );

};
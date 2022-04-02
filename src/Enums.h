// Enum.h
// Copyright (c) 2022, Renier Cloete
// This program is released under the BSD license. See the file LICENSE.txt for details.

#pragma once

enum class eEdgeType
{
	FREE = 0,
	SYMMETRY,
	FIXED,
	SIMPLE_ANCHORED,
	SIMPLE_NONANCHORED,
	KNIFE_EDGE_ANCHORED,
	KNIFE_EDGE_UNANCHORED,
	INTERNAL
};
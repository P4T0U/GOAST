// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef OPTIMIZATION_HH
#define OPTIMIZATION_HH


#include "Optimization/optInterface.h"
#include "Optimization/optParameters.h"
#include "Optimization/optUtils.h"
#include "Optimization/Constraints.h"
#include "Optimization/Objectives.h"
#include "Optimization/LinearOperator.h"
#include "Optimization/Functionals.h"

#include "Core/LinearSolver.h"

#include "Optimization/stepsizeControl.h"

#include "Optimization/gradientDescent.h"
#include "Optimization/quasiNewton.h"
#include "Optimization/GaussNewton.h"
#include "Optimization/simpleNewton.h"
#include "Optimization/TrustRegionNewton.h"
#include "Optimization/LineSearchNewton.h"

#include "Optimization/augmentedLagrange.h"
#include "Optimization/ByrdOmojokunSQP.h"

//=============================================================================
#endif // OPTIMIZATION_HH defined
//=============================================================================

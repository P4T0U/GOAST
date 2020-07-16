// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef NRIC_H
#define NRIC_H

// Basic NRIC stuff, admissibility + LTheta op = parametrization
#include <goast/NRIC/Admissibility.h>
#include <goast/NRIC/NRICMap.h>

// Reconstruction
#include <goast/NRIC/LeastSquaresReconstruction.h>
#include <goast/NRIC/DirectReconstruction.h>

// Deformation energies
#include <goast/NRIC/NonlinearDeformations.h>
#include <goast/NRIC/QuadraticDeformations.h>

// Other
#include <goast/NRIC/GaussCurvature.h>
#include <goast/NRIC/TrigonometryOperators.h>
#include <goast/NRIC/LinearBlending.h>

#endif //NRIC_H

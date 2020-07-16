// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Header including Ipopt and providing helper functions
 * \author Simon
 */
#ifndef __IPOPTINCLUDES_H
#define __IPOPTINCLUDES_H

#ifdef GOAST_WITH_IPOPT

#ifdef __GNUC__
#pragma GCC system_header
#endif

#include <cstddef>


#define HAVE_CSTDDEF

#include <IpTNLP.hpp>

#undef HAVE_CSTDDEF

#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include "IpTNLPAdapter.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpIpoptData.hpp"
#include "IpOrigIpoptNLP.hpp"

/**
 * \brief Yield a string describing an given Ipopt status
 * \author Simon
 * \param ipoptStatus the Ipopt status code
 * \return the describing string
 */
std::string getIpoptStatus( Ipopt::ApplicationReturnStatus ipoptStatus ) {
  std::string status = "";

  switch ( ipoptStatus ) {
    case 0 :
      status = "Solve Succeeded";
      break;

    case 1 :
      status = "Solved To Acceptable Level";
      break;

    case 2 :
      status = "Infeasible Problem Detected";
      break;

    case 3 :
      status = "Search Direction Becomes Too Small";
      break;

    case 4 :
      status = "Diverging Iterates";
      break;

    case 5 :
      status = "User Requested Stop";
      break;

    case 6 :
      status = "Feasible Point Found";
      break;

    case -1 :
      status = "Maximum Iterations Exceeded";
      break;

    case -2 :
      status = "Restoration Failed";
      break;

    case -3 :
      status = "Error In Step Computation";
      break;

    case -10 :
      status = "Not Enough Degrees Of Freedom";
      break;

    case -11 :
      status = "Invalid Problem Definition";
      break;

    case -12 :
      status = "Invalid Option";
      break;

    case -13 :
      status = "Invalid Number Detected";
      break;

    case -100 :
      status = "Unrecoverable Exception";
      break;

    case -101 :
      status = "NonIpopt Exception Thrown";
      break;

    case -102 :
      status = "Insufficient Memory";
      break;

    case -199 :
      status = "Internal Error";
      break;
  }

  return status;

}

/**
 * \brief Print Ipopt status in a legible way
 * \param ipoptStatus the Ipopt status code
 * \param outputOnlyIfFailed whether the output should only be printed if it describes an unsuccessful scenario
 */
void outputIpoptStatus( Ipopt::ApplicationReturnStatus ipoptStatus, const bool outputOnlyIfFailed = false ) {
  std::string status = getIpoptStatus( ipoptStatus );
  // status: 0 - Solve_Succeeded, 1 - Solved_To_Acceptable_Level
  if ( !outputOnlyIfFailed && (ipoptStatus == 0 || ipoptStatus == 1))
    std::cout << "Ipopt finished with status: " << status.c_str() << std::endl;
}


#endif //GOAST_WITH_IPOPT

#endif //__IPOPTINCLUDES_H
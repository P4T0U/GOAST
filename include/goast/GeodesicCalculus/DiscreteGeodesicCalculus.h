// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Collection of all important functions to perform geodesic interpolation and extrapolation.
 * \author Heeren
 * 
 */
#ifndef GEODESICCALCULUS_HH
#define GEODESICCALCULUS_HH

#include <goast/Core.h>
#include <goast/Optimization.h>

#include "ConvexGeodesicCombination.h"
#include "GeodesicInterpolation.h"
#include "GeodesicExtrapolation.h"

#include <iostream>


//======================================================================================================================
// RIEMANNIAN CONVEX COMBINATIONS
//======================================================================================================================

/**
 * \brief Compute weighted 3-point geodesic interpolation
 * \author Heeren
 *
 * Given shapes \f$ S_A \f$ and \f$ S_B \f$ and some real value \f$ 0 < \lambda < 1 \f$,
 * the convex combination \f$ S \f$ is defined as the minimizer of \f$ (1-\lambda) W[S_A, S] +\lambda W[S, S_B] \f$
 * If the mask vector size is larger than zero, a boundary mask with Dirichletnodes is used in the optimization.
 * Otherwise, the optimization is done by means of a Lagrange ansatz fixing the first two modes prescribed by shape \f$ S_A \f$.
 */
template<typename ConfiguratorType>
void convexInterpolation( const MeshTopologySaver& Topology, 
                          const typename ConfiguratorType::VectorType& StartGeom, 
                          const typename ConfiguratorType::VectorType& EndGeom,
                          const std::vector<int>& Mask,
                          const DeformationBase<ConfiguratorType>& W,                          
                          const OptimizationParameters<ConfiguratorType>& optPars,
                          double lambda,
                          typename ConfiguratorType::VectorType& Result,
                          bool quiet = true ){
    ConvexInterpolationOp<ConfiguratorType> convexOp( Topology, StartGeom, EndGeom, W, lambda, optPars, quiet );
    if( Mask.size() > 0 )
      convexOp.setBoundaryMask( Mask );
    convexOp.execute( Result );    
}

/**
 * \brief Compute weighted 3-point geodesic extrapolation
 * \author Heeren
 *
 * Given shapes \f$ S_A \f$ and \f$ S_v \f$ and some real value \f$ 0 < \lambda < 1 \f$,
 * the convex extrapolation is given as the root of the functional \f$ F[S] = (1-\lambda) W_{,2}[S_A, S_v] +\lambda W_{,1}[S_v, S] \f$
 * Note that \f$ F[S] = 0 \f$ is a necessary condition for \f$ (S_A, S_v, S) \f$ being a Riemannian convex combination.
 */
template<typename ConfiguratorType>
void convexExtrapolation( const MeshTopologySaver& Topology, 
                          const typename ConfiguratorType::VectorType& StartGeom, 
                          const typename ConfiguratorType::VectorType& varGeom,
                          const std::vector<int>& Mask,
                          const DeformationBase<ConfiguratorType>& W,                          
                          const OptimizationParameters<ConfiguratorType>& optPars,
                          double lambda,
                          typename ConfiguratorType::VectorType& Result,
                          bool quiet = true ){
    ConvexExtrapolationOp<ConfiguratorType> convexOp( Topology, StartGeom, varGeom, W, lambda, optPars, quiet );
    if( Mask.size() > 0 )
      convexOp.setBoundaryMask( Mask );
    convexOp.execute( Result );    
}

//======================================================================================================================
// GEODESIC INTERPOLATON
//======================================================================================================================

/**
 * \brief Compute full discrete geodesic path
 * \author Heeren
 *
 * In detail, compute discrete geodesic \f$ (s_0, ..., s_K) \f$, where \f$ K+1 \f$ is the total length.
 * Here \f$ s_0 \f$ and \f$ s_K \f$ are given as arguments and the result is given as Path in \f$(s_1, ..., s_{K-1})\f$
 *
 * \note If Path.size() > 0 it is assumed that it contains some reasonable initialization with \f$ K-1 \f$ shapes.
 * If the mask vector size is larger than zero, a boundary mask with Dirichlet nodes is used in the optimization.
 * Otherwise, the optimization is done by means of a Lagrange ansatz fixing the first two modes prescribed by shape \f$ s_0 \f$.
 */
template<typename ConfiguratorType>
void integerInterpolation( const MeshTopologySaver& Topology, 
                           const typename ConfiguratorType::VectorType& StartGeom, 
                           const typename ConfiguratorType::VectorType& EndGeom,
                           const std::vector<int>& Mask,
                           const DeformationBase<ConfiguratorType>& W,
                           const OptimizationParameters<ConfiguratorType>& optPars,
                           int totalLength,
                           typename ConfiguratorType::VectorType& Path,
                           bool quiet = true ){
  GeodesicInterpolation<ConfiguratorType> interpolationOp( Topology, StartGeom, EndGeom, W, totalLength, optPars, quiet );  
  if( Mask.size() > 0 )
    interpolationOp.setBoundaryMask( Mask );
  int initScheme = Path.size() > 0 ? -1 : optPars.getInitializationScheme(); 
  interpolationOp.execute( Path, initScheme );    
}

/**
 * \brief Compute single shape in discrete geodesic path
 * \author Heeren
 *
 * Compute single shape \f$ s_k \f$ in disrete geodesic \f$ (s_0, ..., s_K)\f$, where \f$ K+1 \f$ is the total length and \f$ k \f$ the evaluation
 *
 * \note If initPath.size() > 0 it is assumed that it contains some reasonable initialization with K-1 shapes.
 */
template<typename ConfiguratorType>
void integerInterpolation( const MeshTopologySaver& Topology, 
                           const typename ConfiguratorType::VectorType& StartGeom, 
                           const typename ConfiguratorType::VectorType& EndGeom,
                           const std::vector<int>& Mask,
                           const DeformationBase<ConfiguratorType>& W,
                           const OptimizationParameters<ConfiguratorType>& optPars,
                           int totalLength,
                           int evaluation,
                           const typename ConfiguratorType::VectorType& initPath,
                           typename ConfiguratorType::VectorType& Shape,
                           bool quiet = true ){
  // initial shape
  if( evaluation == 0 ){
    Shape = StartGeom;
    return;    
  }
    
  // totalLength = K+1
  if( evaluation == totalLength - 1 ){
      Shape = EndGeom;
      return;
  }
    
  // compute geodesic  
  typename ConfiguratorType::VectorType path( initPath );
  integerInterpolation<ConfiguratorType>( Topology, StartGeom, EndGeom, Mask, W, optPars, totalLength, path, quiet );  
  
  // write to shape
  Shape.resize( 3 * Topology.getNumVertices() );
  for( int i = 0; i < 3 * Topology.getNumVertices(); i++ )
      Shape[i] = path[3 * (evaluation-1) * Topology.getNumVertices() + i ];    
}

/**
 * \brief Compute single shape in between shapes of discrete geodesic path
 * \author Heeren
 *
 * Compute single shape  \f$ s \approx s(t) \f$ along disrete geodesic at some evaluation time \f$ 0 \leq t \leq 1 \f$.
 * First, compute discrete geodesic \f$ (s_0, ..., s_K) \f$ by using the function integerInterpolation above.
 * Second, compute convex combination between \f$ s_k \f$ and \f$ s_{k+1} \f$ where  \f$ k = \max\{ n \in \N : n \leq tK \} \f$ with an appropiate value \f$ 0 < \lambda< 1 \f$.
 *
 * \note If initPath.size() > 0 it is assumed that it contains some reasonable initialization with K-1 shapes.
 */
template<typename ConfiguratorType>
void doubleInterpolation( const MeshTopologySaver& Topology, 
                          const typename ConfiguratorType::VectorType& StartGeom, 
                          const typename ConfiguratorType::VectorType& EndGeom,
                          const std::vector<int>& Mask,
                          const DeformationBase<ConfiguratorType>& W,
                          const OptimizationParameters<ConfiguratorType>& optPars,
                          int totalLength,
                          double evaluation,
                          typename ConfiguratorType::VectorType& initPath,
                          typename ConfiguratorType::VectorType& Shape,
                          bool quiet = true ){

  integerInterpolation<ConfiguratorType>( Topology, StartGeom, EndGeom, Mask, W, optPars, totalLength, initPath, quiet );  
  
  // get left and right shape
  int K = totalLength - 1;
  int leftIdx = std::floor( evaluation * K );
  int rightIdx = leftIdx + 1;
  
  typename ConfiguratorType::VectorType leftShape( StartGeom ), rightShape( EndGeom );
  if( leftIdx > 0 )
    for( int i = 0; i < 3 * Topology.getNumVertices(); i++ )
      leftShape[i]  = initPath[3 * (leftIdx-1) * Topology.getNumVertices() + i ]; 
  if( rightIdx < K )
    for( int i = 0; i < 3 * Topology.getNumVertices(); i++ )
      rightShape[i] = initPath[3 * (rightIdx-1) * Topology.getNumVertices() + i ];
  
  // compute convex combination
  double lambda = evaluation*K - leftIdx;
  if(!quiet) std::cerr << "Compute convex combination between l = " << leftIdx << " and r = " << rightIdx << " with lambda = " << lambda << std::endl;
  convexInterpolation<ConfiguratorType>( Topology, leftShape, rightShape, Mask, W, optPars, lambda, Shape, quiet );

}

//======================================================================================================================
// GEODESIC EXTRAPOLATON
//======================================================================================================================

/**
 * \brief Compute sequence of discrete geodesic extrapolation.
 * \author Heeren
 *
 * Given \f$ S_A \f$ and \f$ S_v \f$, compute a discrete geodesic \f$ (S_0, S_1, S_2, ..., S_K) \f$ where \f$ K-1 \f$ is the number of extrapolation steps.
 *  a) If forwardShooting = true: \f$ S_0 = S_A \f$ and  \f$ S_1 = S_v \f$
 *     For illustration with K=2:  Given start shape \f$ S_0 \f$ and variational shape \f$ S_1 \f$, compute \f$ S_2 \f$ s.t. \f$ (S_0, S_1, S_2) \f$ is geodesic.
 *  b) If forwardShooting = false: \f$ S_{K-1} = S_A \f$ and \f$ S_K = S_v \f$
 *     For illustration with K=2: given start shape \f$ S_1 \f$ and variational shape \f$ S_2 \f$, compute \f$ S_0 \f$ s.t. \f$(S_0, S_1, S_2)\f$ is geodesic.
 *
 * \note Result is stored as std::vector where each element is a single shape \f$ S_k \f$ (see also next function!)
 */
template<typename ConfiguratorType>
void integerExtrapolation( const MeshTopologySaver& Topology, 
                           const typename ConfiguratorType::VectorType& StartGeom, 
                           const typename ConfiguratorType::VectorType& varGeom,
                           const std::vector<int>& Mask,
                           const DeformationBase<ConfiguratorType>& W,
                           const OptimizationParameters<ConfiguratorType>& optPars,
                           int steps,
                           std::vector<typename ConfiguratorType::VectorType>& Path,
                           bool forwardShooting = true,
                           bool quiet = true ){
  GeodesicExtrapolation<ConfiguratorType> extrapolationOp( Topology, StartGeom, varGeom, W, optPars, quiet );  
  if( Mask.size() > 0 )
    extrapolationOp.setBoundaryMask( Mask );
  extrapolationOp.execute( steps, Path, forwardShooting );     
}

/**
 * \brief Compute sequence of discrete geodesic extrapolation.
 * \author Heeren
 *
 * Same function as integerExtrapolation above but here the result is stored as concatenation in single vector
 */
template<typename ConfiguratorType>
void integerExtrapolation( const MeshTopologySaver& Topology, 
                           const typename ConfiguratorType::VectorType& StartGeom, 
                           const typename ConfiguratorType::VectorType& varGeom,
                           const std::vector<int>& Mask,
                           const DeformationBase<ConfiguratorType>& W,
                           const OptimizationParameters<ConfiguratorType>& optPars,
                           int steps,
                           typename ConfiguratorType::VectorType& Path,
                           bool forwardShooting = true,
                           bool quiet = true ){
  GeodesicExtrapolation<ConfiguratorType> extrapolationOp( Topology, StartGeom, varGeom, W, optPars, quiet );  
  if( Mask.size() > 0 )
    extrapolationOp.setBoundaryMask( Mask );
  extrapolationOp.execute( steps, Path, forwardShooting );     
}
#endif

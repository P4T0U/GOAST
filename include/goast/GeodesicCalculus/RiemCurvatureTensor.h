// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Collection of functions to compute covariant derivatives, Riemann curvature tensor and sectional curvatures
 * \author Heeren
 * 
 */
#ifndef RIEMCURVATURETENSOR_HH
#define RIEMCURVATURETENSOR_HH

#include <goast/Core.h>
#include <goast/Optimization.h>
#include <goast/GeodesicCalculus.h>
#include <goast/DiscreteShells/DeformationTransfer.h>

#include "DiscreteGeodesicCalculus.h"

#include <iostream>

//======================================================================================================================
// VECTOR FIELDS (NOT USED SO FAR!!!)
//======================================================================================================================
//
template<typename ConfiguratorType>
class AbstractVectorField{
  typedef typename ConfiguratorType::VectorType VectorType;
public:
AbstractVectorField( ) {}
 // Destroy polymorphic class correctly, important!
virtual ~AbstractVectorField () {}
virtual void evaluate( const VectorType& point, VectorType& eval ) const = 0;
};

//
template<typename ConfiguratorType>
class ConstantVectorField : public AbstractVectorField<ConfiguratorType>{
    
typedef typename ConfiguratorType::VectorType VectorType;
const VectorType& _constField;
    
ConstantVectorField( const VectorType& constField ) : _constField(constField) {}

void evaluate( const VectorType& /*point*/, VectorType& eval ) const {
  eval = _constField;     
}
    
};

//
template<typename ConfiguratorType>
class MeanCurvatureVectorField : public AbstractVectorField<ConfiguratorType>{
    
typedef typename ConfiguratorType::VectorType VectorType;

const MeshTopologySaver& _topology;
    
MeanCurvatureVectorField( const MeshTopologySaver& Topology ) : _topology(Topology){}

void evaluate( const VectorType& point, VectorType& eval ) const {
  //
  if( point.size() != 3 * _topology.getNumVertices() )
      throw BasicException("MeanCurvatureVectorField::evaluate(): sizes don't match!");
  
  //
}
    
};


//======================================================================================================================
// DISCRETE PARALLEL TRANSPORT
//======================================================================================================================

/**
 * \brief Discrete version of parallel transport (via one step of Schild's ladder)
 * \author Heeren
 *
 * Let \f$ y_0 \f$ and \f$ y_1 \f$ a start and end geometry (i.e shapes), respectively, \f$ \xi_0 \f$ a tangent vector at \f$ y_0 \f$ representing a variation of \f$ y_0 \f$.
 *
 * Let \f$ y_0^p = y_0 + \xi_0 \f$ the corresponding "variational start geometry", then we are looking for a "variational end geometry" \f$ y_1^p = y_1 + \xi_1 \f$,
 * such that \f$ \xi_1 = P_{y_0\to y_1}(\xi_0) \f$ is the vector \f$ \xi_0 \f$ transported parallely from \f$ y_0\f$  to \f$ y_1 \f$ via Schild's ladder construction.
 *
 * One iteration of Schild's ladder consists of two steps:
 *
 * i) Compute midpoint \f$ y^c \f$ of three-point geodesic \f$ (y_0^p, y^c, y_1) \f$.
 *
 * ii) Compute \f$ y_1^p \f$ via discrete geodesic extrapolation, i.e. a three-point geodesic \f$ (y_0, y^c, y_1^p) \f$.
 */
template<typename ConfiguratorType>
void parTranspSingleStepShort( const MeshTopologySaver& Topology, 
                                     const typename ConfiguratorType::VectorType& StartGeom,
                                     const typename ConfiguratorType::VectorType& EndGeom,
                                     const typename ConfiguratorType::VectorType& VarStartGeom,
                                     const std::vector<int>& Mask,
                                     const DeformationBase<ConfiguratorType>& W,
                                     const OptimizationParameters<ConfiguratorType>& optPars,
                                     typename ConfiguratorType::VectorType& VarEndGeom,
                                     bool quiet = true ){
  typename ConfiguratorType::VectorType SchildLadderMidpoint;
  if( !quiet ) std::cerr << "1) Compute discrete 3-geodesic to get the midpoint."  << std::endl;
  integerInterpolation( Topology, VarStartGeom, EndGeom, Mask, W, optPars, 3, SchildLadderMidpoint, quiet );
  if( !quiet ) std::cerr << "2) Compute discrete EXP2 to get the transported point." << std::endl;
  integerExtrapolation( Topology, StartGeom, SchildLadderMidpoint, Mask, W, optPars, 1, VarEndGeom, true, quiet );
}

/**
 * \brief Discrete version of inverse parallel transport (via one step of Schild's ladder)
 * \author Heeren
 *
 * Let \f$ y_0 \f$ and \f$ y_1 \f$ a start and end geometry (i.e shapes), respectively, \f$ \xi_1 \f$ a tangent vector at \f$ y_1 \f$ representing a variation of \f$ y_1 \f$.
 *
 * Let \f$ y_1^p = y_1 + \xi_1 \f$ the corresponding "variational end geometry", then we are looking for a "variational start geometry" \f$ y_0^p = y_0 + \xi_0 \f$,
 * such that \f$ \xi_1 = P_{y_0\to y_1}(\xi_0) \f$ is the vector \f$ \xi_0 \f$ transported parallely from \f$ y_0 \f$ to \f$ y_1 \f$ via Schild's ladder construction.
 *
 * One iteration of INVERSE Schild's ladder consists of two steps:
 *
 * i) Compute midpoint \f$ y^c \f$ of three-point geodesic \f$(y_0, y^c, y_1^p)\f$.
 *
 * ii) Compute \f$ y_0^p \f$ via (BACKWARD!) discrete geodesic extrapolation, i.e. a three-point geodesic \f$(y_0^p, y^c, y_1)\f$.
 *
 * \note Since W is not symmetric, \f$ y_0^p \f$ is NOT obtained via forward discrete extrapolation, i.e. computing \f$ y \f$ such that \f$(y_1, y^c, y)\f$ is a three-point geodesic!
 */
template<typename ConfiguratorType>
void invParTranspSingleStepShort( const MeshTopologySaver& Topology, 
                                  const typename ConfiguratorType::VectorType& StartGeom,
                                  const typename ConfiguratorType::VectorType& EndGeom,
                                  const typename ConfiguratorType::VectorType& VarEndGeom,
                                  const std::vector<int>& Mask,
                                  const DeformationBase<ConfiguratorType>& W,
                                  const OptimizationParameters<ConfiguratorType>& optPars,
                                  typename ConfiguratorType::VectorType& VarStartGeom,
                                  bool quiet = true,
                                  std::string savenameStem = "" ){  
  if( !quiet ) std::cerr << "Start inverse parallel transport with single step (i.e. K=1)."  << std::endl;
  typename ConfiguratorType::VectorType SchildLadderMidpoint;
  if( !quiet ) std::cerr << "1) Compute discrete 3-geodesic to get the midpoint."  << std::endl;
  integerInterpolation( Topology, StartGeom, VarEndGeom, Mask, W, optPars, 3, SchildLadderMidpoint, quiet );
  if( !quiet ) std::cerr << "2) Compute discrete EXP2 to get the inverse transported point." << std::endl;
  integerExtrapolation( Topology, SchildLadderMidpoint, EndGeom,  Mask, W, optPars, 1, VarStartGeom, false, quiet );
  
  // save parallelogram 
  if( savenameStem.size() > 0 ){
      TriMesh mesh( Topology.getGrid() );
      
      std::ostringstream basePoint;
      basePoint<< savenameStem << "_basePoint.ply";
      setGeometry( mesh, StartGeom );
      OpenMesh::IO::write_mesh( mesh, basePoint.str() );
      
      std::ostringstream endPoint;
      endPoint<< savenameStem << "_endPoint.ply";
      setGeometry( mesh, EndGeom );
      OpenMesh::IO::write_mesh( mesh, endPoint.str() );
      
      std::ostringstream varEnd;
      varEnd<< savenameStem << "_varEnd.ply";
      setGeometry( mesh, VarEndGeom );
      OpenMesh::IO::write_mesh( mesh, varEnd.str() );
      
      std::ostringstream varStart;
      varStart<< savenameStem << "_varStart.ply";
      setGeometry( mesh, VarStartGeom );
      OpenMesh::IO::write_mesh( mesh, varStart.str() );
      
      std::ostringstream saveNameMidpoint;
      saveNameMidpoint<< savenameStem << "_midPoint.ply";
      setGeometry( mesh, SchildLadderMidpoint );
      OpenMesh::IO::write_mesh( mesh, saveNameMidpoint.str() );
  }
}


//======================================================================================================================
// DIFFERENT VERSIONS OF DISCRETE FIRST COVARIANT DERIVATIVE
//======================================================================================================================

/**
 * \brief First-order consistent approximation of first covariant derivative via FORWARD difference quotient
 * \author Heeren
 *
 * Computes at StartGeom the first covariant derivative of (constant) tangent vector field Eta in tangent direction U via forward difference quotient of step size Tau.
 *
 * The discrete inverse parallel transport is realized by means of a scaling factor Epsilon (where usually Tau = Epsilon)
 *
 * \note Eta, U, FirCovDer are shape differences (i.e. tangent vectors/displacements), not shapes themselves!
 */
template<typename ConfiguratorType>
void firstCovariantDerivativeForward( const MeshTopologySaver& Topology, 
                                    const typename ConfiguratorType::VectorType& StartGeom,
                                    const typename ConfiguratorType::VectorType& Eta,
                                    const typename ConfiguratorType::VectorType& U,
                                    const double Tau,
                                    const double Epsilon,
                                    const std::vector<int>& Mask,
                                    const DeformationBase<ConfiguratorType>& W,
                                    const OptimizationParameters<ConfiguratorType>& optPars,
                                    typename ConfiguratorType::VectorType& FirCovDer,
                                    bool quiet = true ){  
  typename ConfiguratorType::VectorType yURight( StartGeom );
  yURight += Tau * U;
    
  // Here we actually use that Eta is constant. Otherwise we need evaluations at this point!   
  typename ConfiguratorType::VectorType yURightEta( yURight );                                   
  yURightEta += Epsilon * Eta;
                                   
  // perform two inverse parallel transports on RIGHT hand side
  typename ConfiguratorType::VectorType xiLeft;  // xiRight is stored in FirCovDer
  invParTranspSingleStepShort( Topology, StartGeom, yURight, yURightEta, Mask, W, optPars, FirCovDer, quiet );  
  
  // subtract position (i.e. reference shape) to get tangent vectors (i.e. displacements)
  // However, we have to register before subtraction!
  // NOTE Since StartGeom is supposed to be registered wrt. the global reference shape (stored in Topology) we also register to this reference shape here.
  if( Mask.size() == 0 ){
    if( !quiet ) std::cerr << "Perform registration with respect to reference mesh (prescribed by topology class)..." << std::endl;
    typename ConfiguratorType::VectorType referenceGeom, OptimalRBMdofs;
    getGeometry( Topology.getGrid(), referenceGeom );
    MomentumRegistrationOperator<ConfiguratorType>( Topology, referenceGeom ).execute( FirCovDer, OptimalRBMdofs );    
  }  
  FirCovDer -= StartGeom;
  
  FirCovDer -= Epsilon * Eta;
  FirCovDer /= Tau * Epsilon;
}

/**
 * \brief First-order consistent approximation of first covariant derivative via CENTRAL difference quotient
 * \author Heeren
 *
 * Computes at StartGeom the first covariant derivative of (constant) tangent vector field Eta in tangent direction U via symmetric differences of step size Tau.
 *
 * The discrete inverse parallel transport is realized by means of a scaling factor Epsilon (where usually Tau = Epsilon)
 *
 * Although a central difference quotient is used, the order is still linear in tau due to missing reflection! (cf. to function firstCovariantDerivative below)
 *
 * \note Eta, U, FirCovDer are shape differences (i.e. tangent vectors/displacements), not shapes themselves!
 */
template<typename ConfiguratorType>
void firstCovariantDerivativeSymm( const MeshTopologySaver& Topology, 
                                const typename ConfiguratorType::VectorType& StartGeom,
                                const typename ConfiguratorType::VectorType& Eta,
                                const typename ConfiguratorType::VectorType& U,
                                const double Tau,
                                const double Epsilon,
                                const std::vector<int>& Mask,
                                const DeformationBase<ConfiguratorType>& W,
                                const OptimizationParameters<ConfiguratorType>& optPars,
                                typename ConfiguratorType::VectorType& FirCovDer,
                                bool quiet = true ){
  if( !quiet ) std::cerr << "--------------------------------------------------------------" << std::endl;
  if( !quiet ) std::cerr << "COMPUTE FIRST COVARIANT DERIVATIVE" << std::endl;
  if( !quiet ) std::cerr << "--------------------------------------------------------------" << std::endl;
  
  typename ConfiguratorType::VectorType yURight( StartGeom ), yULeft( StartGeom );
  yURight += Tau * U;
  yULeft  -= Tau * U;
    
  // Here we actually use that Eta is constant. Otherwise we need evaluations at this point!   
  typename ConfiguratorType::VectorType yURightEta( yURight ), yULeftEta( yULeft );                                   
  yURightEta += Epsilon * Eta;
  yULeftEta  += Epsilon * Eta;
                                   
  // perform two inverse parallel transports on RIGHT hand side
  typename ConfiguratorType::VectorType xiLeft;  // xiRight is stored in FirCovDer
  invParTranspSingleStepShort( Topology, StartGeom, yURight, yURightEta, Mask, W, optPars, FirCovDer, quiet );
  invParTranspSingleStepShort( Topology, StartGeom, yULeft,  yULeftEta, Mask, W, optPars, xiLeft, quiet );
  
  FirCovDer -= xiLeft;
  FirCovDer /= 2 * Tau * Epsilon;
}

/**
 * \brief Second-order consistent approximation of first covariant derivative via CENTRAL and REFLECTED difference quotient
 * \author Heeren
 *
 * Computes at StartGeom the first covariant derivative of (constant) tangent vector field Eta in tangent direction U via symmetric differences of step size Tau.
 *
 *  The discrete inverse parallel transport is realized by means of a scaling factor Epsilon (where usually Tau = Epsilon)
 *
 *  Important fact for proper reflection: in discrete world we have \f$ P^{-1}(-\eta) \f$ IS NOT EQUAL \f$ -P^{-1}(\eta)\f$.
 *  Hence central difference quotient given by \f[ \frac{P^{-1}(\epsilon \eta^r) - (-P^{-1}(-\epsilon \eta^l))}{2\tau\epsilon}\, , \f]
 *  where \f$ \eta^r = y + \tau u \f$ and  \f$ \eta^l = y - \tau u \f$.
 *
 *  \note Eta, U, FirCovDer are shape differences (i.e. tangent vectors/displacements), not shapes themselves!
 */
template<typename ConfiguratorType>
void firstCovariantDerivative( const MeshTopologySaver& Topology, 
                               const typename ConfiguratorType::VectorType& StartGeom,
                               const typename ConfiguratorType::VectorType& Eta,
                               const typename ConfiguratorType::VectorType& U,
                               const double Tau,
                               const double Epsilon,
                               const std::vector<int>& Mask,
                               const DeformationBase<ConfiguratorType>& W,
                               const OptimizationParameters<ConfiguratorType>& optPars,
                               typename ConfiguratorType::VectorType& FirCovDer,
                               bool quiet = true,
                               std::string savenameStem = "" ){  
  
  if( !quiet ) std::cerr << "\n----------------------------------------------------------------------------------" << std::endl;
  if( !quiet ) std::cerr << "START FIRST COVARIANT DERIVATIVE" << std::endl;
  typename ConfiguratorType::VectorType yURight( StartGeom ), yULeft( StartGeom );
  yURight += Tau * U;
  yULeft  -= Tau * U;
    
  // Here we actually use that Eta is constant. Otherwise we need evaluations at this point!   
  typename ConfiguratorType::VectorType yURightEta( yURight ), yULeftEta( yULeft );                                   
  yURightEta += Epsilon * Eta;
  yULeftEta  -= Epsilon * Eta; //NOTE the sign here, this is where the proper reflection happens!
  
  // save covariant derivative? 
  std::ostringstream rightName, leftName;
  if( savenameStem.size() > 0 ){
      rightName << savenameStem << "_right";
      leftName<< savenameStem << "_left";
  }
  
  // perform two inverse parallel transports
  typename ConfiguratorType::VectorType xiLeft;  // xiRight is stored in FirCovDer
  // NOTE inverse parallel transport returns points (i.e. shapes) rather than tangent vectors (i.e displacements)
  if( !quiet ) std::cerr << "\n----------------------------------------------------------------------------------" << std::endl;
  if( !quiet ) std::cerr << "Inverse parallel transport on RIGHT hand side." << std::endl;
  invParTranspSingleStepShort( Topology, StartGeom, yURight, yURightEta, Mask, W, optPars, FirCovDer, quiet, rightName.str() );
  if( !quiet ) std::cerr << "\n----------------------------------------------------------------------------------" << std::endl;
  if( !quiet ) std::cerr << "Inverse parallel transport on LEFT hand side." << std::endl;
  invParTranspSingleStepShort( Topology, StartGeom, yULeft,  yULeftEta, Mask, W, optPars, xiLeft, quiet, leftName.str() );  
  
  // subtract position (i.e. reference shape) to get tangent vectors (i.e. displacements)
  // However, we have to register before subtraction!
  // NOTE Since StartGeom is supposed to be registered wrt. the global reference shape (stored in Topology) we also register to this reference shape here.
  if( Mask.size() == 0 ){
    if( !quiet ) std::cerr << "Perform registration with respect to reference mesh (prescribed by topology class)..." << std::endl;
    typename ConfiguratorType::VectorType referenceGeom;
    getGeometry( Topology.getGrid(), referenceGeom );
    MomentumRegistrationOperator<ConfiguratorType>( Topology, referenceGeom ).execute( FirCovDer );    
    MomentumRegistrationOperator<ConfiguratorType>( Topology, referenceGeom ).execute( xiLeft ); 
  }
  FirCovDer -= StartGeom;
  xiLeft -= StartGeom;
  
  // compute sum (!!) of displacements due to change in sign in transported vector (see above)
  FirCovDer += xiLeft;
  FirCovDer /= 2 * Tau * Epsilon;
  if( !quiet ) std::cerr << "----------------------------------------------------------------------------------\n" << std::endl;
}



//======================================================================================================================
// DISCRETE SECOND COVARIANT DERIVATIVE, RIEMANN CURVATURE TENSOR AND SECTIONAL CURVATURE
//======================================================================================================================

/**
 * \brief Approximation of second covariant derivative by means of nested first covariant derivatives (quadratic convergence in tau!)
 * \author Heeren
 *
 * Computes at StartGeom the second covariant derivative \f$\nabla_V \nabla_U \eta\f$ of (constant) tangent vector field Eta in tangent directions U and V.
 *
 * In detail, the first covariant derivative is applied twice, i.e. we distinguish between outer derivative (here, V) and inner derivative (here, U).
 * Both derivatives are approximated with (symmetrized,reflected) central difference quotient,
 * where we use different step sizes tau for the outer and inner derivative, respectively.
 *
 * For quadratic convergence, however, we need that \f$innerTau = outerTau^\beta\f$ with \f$\beta \geq 3/2\f$.
 *
 * \note Eta, U, V, SecCovDer are shape differences, not shapes themselves!
 */
template<typename ConfiguratorType>
void secondCovariantDerivative( const MeshTopologySaver& Topology, 
                                const typename ConfiguratorType::VectorType& StartGeom,
                                const typename ConfiguratorType::VectorType& Eta,
                                const typename ConfiguratorType::VectorType& outerDir,
                                const typename ConfiguratorType::VectorType& innerDir,                                
                                const double OuterTau,
                                const double InnerTau,
                                const std::vector<int>& Mask,
                                const DeformationBase<ConfiguratorType>& W,
                                const OptimizationParameters<ConfiguratorType>& optPars,
                                typename ConfiguratorType::VectorType& SecCovDer,
                                bool quiet = true,
                                std::string savenameStem = "" ){
  if( !quiet ) std::cerr << "\n\n----------------------------------------------------------------------------------" << std::endl;
  if( !quiet ) std::cerr << "START OF SECOND COVARIANT DERIVATIVE" << std::endl;
  if( !quiet ) std::cerr << "----------------------------------------------------------------------------------" << std::endl;
  
  TriMesh auxMesh( Topology.getGrid() );
  typename ConfiguratorType::RealType outerEps = OuterTau;
  typename ConfiguratorType::RealType innerEps = InnerTau;
  
  // let V be the outer direction
  typename ConfiguratorType::VectorType yVRight( StartGeom ), yVLeft( StartGeom );
  yVRight += OuterTau * outerDir;
  yVLeft  -= OuterTau * outerDir;
  
  // save meshes?
  if( savenameStem.size() > 0 ){
      if( !quiet ) std::cerr << "Save meshes..." << std::endl;
      std::ostringstream savenameRight, savenameLeft;
      savenameRight << savenameStem << "_yVRight.ply";
      setGeometry( auxMesh, yVRight );
      OpenMesh::IO::write_mesh(auxMesh, savenameRight.str() );  
      savenameLeft << savenameStem << "_yVLeft.ply";
      setGeometry( auxMesh, yVLeft );
      OpenMesh::IO::write_mesh(auxMesh, savenameLeft.str() ); 
  }
  
  // save inner covariant derivatives? 
  std::ostringstream innerRightName, innerLeftName;
  if( savenameStem.size() > 0 ){
      innerRightName << savenameStem << "_innerRight";
      innerLeftName<< savenameStem << "_innerLeft";
  }
  
  // compute first covariant derivative at y_r = y + tau * V and y_l = y - tau * V, denoted by \xi_r and \xi_l, respectively
  // let \xi denote the vector field representing the inner covariant derivative.
  typename ConfiguratorType::VectorType xiRight, xiLeft;
  if( !quiet ) std::cerr << "\n\n----------------------------------------------------------------------------------" << std::endl;
  if( !quiet ) std::cerr << "Start two INNER first covariant derivatives (RIGHT and LEFT)." << std::endl;
  if( !quiet ) std::cerr << "----------------------------------------------------------------------------------" << std::endl;
  if( !quiet ) std::cerr << "Compute RIGHT INNER first covariant derivative with tau = " << InnerTau  << std::endl;
  firstCovariantDerivative( Topology, yVRight, Eta, innerDir, InnerTau, innerEps, Mask, W, optPars, xiRight, quiet, innerRightName.str() );
  if( !quiet ) std::cerr << "\n\n----------------------------------------------------------------------------------" << std::endl;
  if( !quiet ) std::cerr << "----------------------------------------------------------------------------------" << std::endl;
  if( !quiet ) std::cerr << "Compute LEFT INNER first covariant derivative with tau = " << InnerTau  << std::endl;
  firstCovariantDerivative( Topology, yVLeft, Eta, innerDir, InnerTau, innerEps, Mask, W, optPars, xiLeft, quiet, innerLeftName.str() );

  // get variational shapes from shape differences \xi at y_r and y_l
  typename ConfiguratorType::VectorType yVRightXi( yVRight ), yVLeftXi( yVLeft ), temp;
  yVRightXi += outerEps * xiRight;
  yVLeftXi  -= outerEps * xiLeft; //NOTE the sign here!
  
  // save meshes?
  if( savenameStem.size() > 0 ){
      if( !quiet ) std::cerr << "Save meshes..." << std::endl;
      std::ostringstream savenameRight, savenameLeft;
      savenameRight << savenameStem << "_yVRightXi.ply";
      setGeometry( auxMesh, yVRightXi );
      OpenMesh::IO::write_mesh(auxMesh, savenameRight.str() );  
      savenameLeft << savenameStem << "_yVLeftXi.ply";
      setGeometry( auxMesh, yVLeftXi );
      OpenMesh::IO::write_mesh(auxMesh, savenameLeft.str() ); 
  }
  
  // save outer covariant derivatives? 
  std::ostringstream outerRightName, outerLeftName;
  if( savenameStem.size() > 0 ){
      outerRightName << savenameStem << "_outerRight";
      outerLeftName<< savenameStem << "_outerLeft";
  }
  
  if( !quiet ) std::cerr << "\n\n----------------------------------------------------------------------------------" << std::endl;
  if( !quiet ) std::cerr << "----------------------------------------------------------------------------------" << std::endl;
  if( !quiet ) std::cerr << "Start OUTER first covariant derivative with tau = " << OuterTau << std::endl;
  // NOTE inverse parallel transport returns points (i.e. shapes) rather than tangent vectors (i.e displacements)
  if( !quiet ) std::cerr << "\n----------------------------------------------------------------------------------" << std::endl;
  if( !quiet ) std::cerr << "Inverse parallel transport on RIGHT hand side." << std::endl;
  invParTranspSingleStepShort( Topology, StartGeom, yVRight, yVRightXi, Mask, W, optPars, SecCovDer, quiet, outerRightName.str() );
  if( !quiet ) std::cerr << "\n----------------------------------------------------------------------------------" << std::endl;
  if( !quiet ) std::cerr << "Inverse parallel transport on LEFT hand side." << std::endl;
  invParTranspSingleStepShort( Topology, StartGeom, yVLeft, yVLeftXi, Mask, W, optPars, temp, quiet, outerLeftName.str() );    
  // subtract position (i.e. reference shape) to get tangent vectors (i.e. displacements)
  if( Mask.size() == 0 ){
    if( !quiet ) std::cerr << "Perform registration with respect to reference mesh (prescribed by topology class)..." << std::endl;
    typename ConfiguratorType::VectorType referenceGeom;
    getGeometry( Topology.getGrid(), referenceGeom );
    MomentumRegistrationOperator<ConfiguratorType>( Topology, referenceGeom ).execute( SecCovDer );    
    MomentumRegistrationOperator<ConfiguratorType>( Topology, referenceGeom ).execute( temp ); 
  }
  SecCovDer -= StartGeom;
  temp -= StartGeom;
  
  //throw BasicException("STOP HERE!");
  
  // compute sum (!!) of displacements due to change in sign in transported vector (see above)
  SecCovDer += temp;
  SecCovDer /= 2. * OuterTau * outerEps;
  
  // save meshes?
  if( savenameStem.size() > 0 ){
      if( !quiet ) std::cerr << "Save final mesh..." << std::endl;
      std::ostringstream savename;
      savename << savenameStem << "_final.ply";
      setGeometry( auxMesh, StartGeom + SecCovDer );
      OpenMesh::IO::write_mesh(auxMesh, savename.str() ); 
  }
  
  if( !quiet ) std::cerr << "----------------------------------------------------------------------------------" << std::endl;
  if( !quiet ) std::cerr << "END OF SECOND COVARIANT DERIVATIVE" << std::endl;
  if( !quiet ) std::cerr << "----------------------------------------------------------------------------------\n\n" << std::endl;
}

/**
 * \brief Approximation of Riemannian curvature tensor R(U,V)Eta at StartGeom for (constant) vector fields U, V, Eta.
 * \author Heeren
 *
 * Since \f$\eta\f$ is constant, we have \f$ \nabla_{[U,V]} \eta = 0 \f$ and hence \f[ R(U,V)\eta = \nabla_U \nabla_V \eta - \nabla_V \nabla_U \eta \, .\f]
 * We use the second covariant derivative with outer step size \f$ \tau \f$ and inner step size \f$ \tau^{3/2} \f$.
 *
 * \note Eta, U, V, RUVEta are shape differences, not shapes themselves!
 */
template<typename ConfiguratorType>
void curvatureTensor( const MeshTopologySaver& Topology, 
                      const typename ConfiguratorType::VectorType& StartGeom,
                      const typename ConfiguratorType::VectorType& Eta,
                      const typename ConfiguratorType::VectorType& U,
                      const typename ConfiguratorType::VectorType& V,
                      const double Tau,
                      const std::vector<int>& Mask,
                      const DeformationBase<ConfiguratorType>& W,
                      const OptimizationParameters<ConfiguratorType>& optPars,
                      typename ConfiguratorType::VectorType& RUVEta,
                      bool quiet = true,
                      std::string savenameStem = "" ){
    
  if( !quiet ) std::cerr << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;  
  if( !quiet ) std::cerr << "COMPUTE CURVATURE TENSOR" << std::endl;
  if( !quiet ) std::cerr << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;  
  typename ConfiguratorType::RealType innerTau = Tau * std::sqrt( Tau );
    
  std::ostringstream firstSecCovDerivFilename, secondSecCovDerivFilename;
  if( savenameStem.size() > 0 ){
      firstSecCovDerivFilename << savenameStem << "_secCovDeriv_firstInnerTerm";
      secondSecCovDerivFilename << savenameStem << "_secCovDeriv_secondInnerTerm";
  }
  
  typename ConfiguratorType::VectorType secondCovarDerVUEta; // secCovDerUV is stored in RUVEta
  if( !quiet ) std::cerr << "i) Compute Nabla_U Nabla_V Eta." << std::endl;
  secondCovariantDerivative( Topology, StartGeom, Eta, U, V, Tau, innerTau, Mask, W, optPars, RUVEta, quiet, firstSecCovDerivFilename.str() );
  if( !quiet ) std::cerr << "ii) Compute Nabla_V Nabla_U Eta." << std::endl;
  secondCovariantDerivative( Topology, StartGeom, Eta, V, U, Tau, innerTau, Mask, W, optPars, secondCovarDerVUEta, quiet, secondSecCovDerivFilename.str() );
  
  // saving
  if( savenameStem.size() > 0 ){
    if( !quiet ) std::cerr << "Save meshes..." << std::endl;
    TriMesh output( Topology.getGrid() );
    std::ostringstream savenameRight, savenameLeft, savename;
    savenameRight << savenameStem << "_secCovDeriv_RyUVRight.ply";    
    setGeometry( output, StartGeom + RUVEta );
    OpenMesh::IO::write_mesh(output, savenameRight.str() );
    savenameLeft << savenameStem << "_secCovDeriv_RyUVLeft.ply";    
    setGeometry( output, StartGeom + secondCovarDerVUEta );
    OpenMesh::IO::write_mesh(output, savenameLeft.str() );
    savename << savenameStem << "_secCovDeriv_RyUVFinal.ply";    
    setGeometry( output, StartGeom + RUVEta - secondCovarDerVUEta );
    OpenMesh::IO::write_mesh(output, savename.str() );
  }
  
  RUVEta -= secondCovarDerVUEta;
}
 
/**
 * \brief Approximation of sectional curvature
 * \author Heeren
 *
 * Sectional curvature \f$ \kappa_{U,V} \f$ at StartGeom in directions U, V, i.e. \f[ \kappa_{U,V} = \frac{g(U,R(U,V)V)}{g(U,U)g(V,V)-g(U,V)^2} \, .\f]
 *
 * \note U, V are shape differences, not shapes themselves!
 */
template<typename ConfiguratorType>
double sectionalCurvature( const MeshTopologySaver& Topology, 
                           const typename ConfiguratorType::VectorType& StartGeom,
                           const typename ConfiguratorType::VectorType& U,
                           const typename ConfiguratorType::VectorType& V,
                           const double Tau,
                           const std::vector<int>& Mask,
                           const DeformationBase<ConfiguratorType>& W,
                           const OptimizationParameters<ConfiguratorType>& optPars,
                           bool quiet = true, 
                           std::string savenameStem = "" ){

  if( !quiet ) std::cerr << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;  
  if( !quiet ) std::cerr << "COMPUTE SECTIONAL CURVATURE" << std::endl;
  if( !quiet ) std::cerr << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;  
  typename ConfiguratorType::VectorType rUVV, metricEval;
  if(!quiet) std::cerr << "Compute R_p(U,V)V..." << std::endl; 
  curvatureTensor( Topology, StartGeom, V, U, V, Tau, Mask, W, optPars, rUVV, quiet, savenameStem );
  
  // g_y = 0.5 D_2^2 W[y,y]
  if(!quiet) std::cerr << "Compute metric g_p = 0.5 * W_{,22}[p,p]... " << std::endl; 
  typename ConfiguratorType::SparseMatrixType Hessian( 3 * Topology.getNumVertices(), 3 * Topology.getNumVertices() );
  W.applyDefHessian( StartGeom, StartGeom, Hessian, 0.5 );
  //if ( Mask.size() > 0 )
  //      applyMaskToSymmetricMatrix( Mask, Hessian );  
    
  //if(!quiet) std::cerr << "Euclidean norm of R(U,V)V is " << rUVV.norm() << std::endl; 
  //metricEval = Hessian * rUVV;
  //if(!quiet) std::cerr << "Shell metric of R(U,V)V is " << rUVV.dot( metricEval ) << std::endl; 
  
  metricEval = Hessian * U;
  typename ConfiguratorType::RealType gUU = U.dot( metricEval );
  metricEval = Hessian * V;  
  typename ConfiguratorType::RealType gVV = V.dot( metricEval );
  typename ConfiguratorType::RealType gUV = U.dot( metricEval );  
  metricEval = Hessian * rUVV;
  
  typename ConfiguratorType::RealType denom = gUU * gVV - gUV * gUV;
  typename ConfiguratorType::RealType secCurv = U.dot( metricEval ) / denom;
  //std::cerr << "Sectional curvature kappa_p(U,V) = g(U, R(U,V)V) / [g(U,U)g(V,V)-g(U,V)^2] = " << U.dot( metricEval ) << " / " << denom << std::endl; 
  if( !quiet ) std::cerr << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;  
  if( !quiet ) std::cerr << "Sectional curvature kappa_p(U,V) = g(U, R(U,V)V) / [g(U,U)g(V,V)-g(U,V)^2] = " << U.dot( metricEval ) << " / ( " << gUU << " * " << gVV << " - (" << gUV << ")^2 ) = "; 
  if( !quiet ) std::cerr << U.dot( metricEval ) << " /  " << denom << " = " << secCurv << std::endl; 
  if( !quiet ) std::cerr << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl << std::endl;  
  
  return  secCurv;
}

//======================================================================================================================
// DISCRETE GEODESIC TRIANGLES
//======================================================================================================================

/**
 * \brief Compute approximation of (sign of) sectional curvature by means of geodesic triangle
 * \author Heeren
 *
 * Consider a geodesic triangle with nodes \f$p\f$, \f$c_1\f$ and \f$c_2\f$.
 * For given EVEN integer \f$K\f$, we first compute a discrete K-geodesic \f$(q_0, ..., q_K)\f$ from \f$c_1\f$ to \f$c_2\f$.
 * The length \f$l_e\f$ of the (curved) edge \f$e = (c_1, c_2)\f$ is given by the length of this geodesic.
 *
 * Next, we compute \f$K+1\f$ discrete K-geodesics from \f$p\f$ to \f$q_k\f$ for \f$k = 0, ..., K\f$ with corresponding lengths \f$l_k\f$.
 * We assume that the geodesic triangle is equal-sided, i.e \f$l_0 = l_K\f$.
 * Since the triangle is equal-sided, the height h_e in the corresponding FLAT triangle \f$(p, c_1, c_2)\f$ wrt. edge e satisfies \f[(l_0)^2 = (h_e)^2 + (0.5 \cdot l_e)^2 \, .\f]
 *
 * The (curved) height d_e of the geodesic triangle wrt. edge is given by \f$l_{K/2}\f$.
 * Now if \f$d_e \approx h_e\f$, we conclude that the sectional curvature \f$ \kappa_p(c_1-p, c_2-p) \f$ is approximately zero.
 * If \f$d_e > h_e\f$, we conclude that \f$ \kappa_p(c_1-p, c_2-p) > 0 \f$ and if \f$d_e < h_e\f$, we conclude that \f$ \kappa_p(c_1-p, c_2-p) < 0 \f$.
 */
template<typename ConfiguratorType>
void investigateGeodesicTriangle( const MeshTopologySaver& Topology, 
                                   const typename ConfiguratorType::VectorType& point,
                                   const typename ConfiguratorType::VectorType& corner1,
                                   const typename ConfiguratorType::VectorType& corner2,
                                   const int K,
                                   const std::vector<int>& Mask,
                                   const DeformationBase<ConfiguratorType>& W,
                                   const OptimizationParameters<ConfiguratorType>& optPars,
                                   bool quiet = true ){
    typedef typename ConfiguratorType::VectorType VectorType;
    typename ConfiguratorType::RealType energy;
    
    if(!quiet) std::cerr << "===========================================================" << std::endl;
    if(!quiet) std::cerr << "Start to investigate geodesic triangle." << std::endl;
    if(!quiet) std::cerr << "===========================================================" << std::endl;
                                                                                                                      
    // compute edge between two corners
    VectorType edgePath;
    if(!quiet) std::cerr << "Compute geodesic along edge..." << std::endl;
    integerInterpolation( Topology, corner1, corner2, Mask, W, optPars, K+1, edgePath, true );
    DiscretePathEnergy<ConfiguratorType>( W, K, corner1, corner2 ).apply( edgePath, energy );
    double lengthEdge = std::sqrt( energy[0] );
    
    // bring into more convenient form
    std::vector< Eigen::Ref<const VectorType> > shapes;
    shapes.reserve(K+1);
    const int numLocalDofs = 3 * Topology.getNumVertices();    
    shapes.push_back( corner1 );
    for( int k = 0; k < K-1; k++ )
       shapes.push_back( edgePath.segment(k*numLocalDofs, numLocalDofs) );
    shapes.push_back( corner2 );
    
    // compute distances from "point" to midpoint on edge    
    VectorType distances(3);
    if(!quiet) std::cerr << "Compute geodesics from edge end- and midpoint to point..." << std::endl;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for( int k = 0; k < 3; k++ ){
        VectorType Path;
        integerInterpolation( Topology, shapes[k*K/2], point, Mask, W, optPars, K+1, Path, true );          
        DiscretePathEnergy<ConfiguratorType>( W, K, shapes[k*K/2], point ).apply( Path, energy );
        distances[k] = std::sqrt( energy[0] );
    }    
    
    if( std::abs( distances[0] - distances[2] ) > 1e-05 )
        std::cerr << "WARNING in investigateGeodesicTriangle(): geodesic triangle is not equal-sided!" << std::endl;
    
    // We consider an equal-sided triangle with edges (e,f,f).
    // Then the height h_e with respect to edge e is given by  |f|^2 = (1/2 * |e|)^2 + |h_e|^2,
    // hence |h_e| = \sqrt{ |f|^2 - 0.25 \cdot |e|^2 }
    double height =  std::sqrt( distances[0] * distances[0] - 0.25 * lengthEdge * lengthEdge );  
    
    // print results to console
    if(!quiet) {
      std::cerr << "\n===================================" << std::endl;
      std::cerr << "Results:" << std::endl;
      std::cerr << "Basis edge length = " << lengthEdge << std::endl;
      std::cerr << "Basis edge height of flat triangle = " << height << std::endl;
      std::cerr << "Basis edge height of curved triangle = " << distances[1] << std::endl;
      if( height > distances[1] )
        std::cerr << "Since flat height is larger, we have a negative curvature.\n\n" << std::endl;
      else
        std::cerr << "Since curved height is larger, we have a positive curvature.\n\n" << std::endl;
    }
}

/**
 * \brief Compute area of flat triangle defined by geodesic lengths between three points in shell space
 * \author Heeren
 *
 * For three points \f$p_0, p_1, p_2\f$ compute the three pairwise geodesic lengths (based on discrete K-geodesics).
 * Return the area of the flat Euclidean triangle induced by these three lengths (by means of Heron's formula).
 */
template<typename ConfiguratorType>
double getAreaOfFlatTriangle( const MeshTopologySaver& Topology, 
                                const typename ConfiguratorType::VectorType& point0,
                                const typename ConfiguratorType::VectorType& point1,
                                const typename ConfiguratorType::VectorType& point2,
                                const int K,
                                const std::vector<int>& Mask,
                                const DeformationBase<ConfiguratorType>& W,
                                const OptimizationParameters<ConfiguratorType>& optPars,
                                bool quiet = true ){
    typedef typename ConfiguratorType::VectorType VectorType;
    typename ConfiguratorType::RealType energy;
    VectorType edgePath;
    
    if(!quiet) std::cerr << "===========================================================" << std::endl;
    if(!quiet) std::cerr << "Start to compute area of FLAT triangle." << std::endl;
    if(!quiet) std::cerr << "===========================================================" << std::endl;
    
    // compute edges and edge lengths    
    if(!quiet) std::cerr << "Compute geodesic of length " << K+1 << " along first edge..." << std::endl;
    integerInterpolation( Topology, point0, point1, Mask, W, optPars, K+1, edgePath, quiet );
    DiscretePathEnergy<ConfiguratorType>( W, K, point0, point1 ).apply( edgePath, energy );
    double l1 = std::sqrt( energy[0] );
    edgePath.resize(0);
    
    if(!quiet) std::cerr << "Compute geodesic of length " << K+1 << " along second edge..." << std::endl;
    integerInterpolation( Topology, point1, point2, Mask, W, optPars, K+1, edgePath, quiet );
    DiscretePathEnergy<ConfiguratorType>( W, K, point1, point2 ).apply( edgePath, energy );
    double l2 = std::sqrt( energy[0] );
    edgePath.resize(0);
    
    if(!quiet) std::cerr << "Compute geodesic of length " << K+1 << " along third edge..." << std::endl;
    integerInterpolation( Topology, point2, point0, Mask, W, optPars, K+1, edgePath, quiet );
    DiscretePathEnergy<ConfiguratorType>( W, K, point2, point0 ).apply( edgePath, energy );
    double l3 = std::sqrt( energy[0] );
    
    // use Heron's formula to compute area
    double s = (l1 + l2 + l3) / 2.;
    if( (s < l1) || (s < l2) || (s < l3) ){
      std::cerr << "getAreaOfFlatTriangle(): computed lengths are " << l1 << " " << l2 << " " << l3 << std::endl;
      std::cerr << "Heron: " << (s-l1) << " " << (s-l2) << " " << (s-l3) << std::endl;
      throw BasicException("getAreaOfFlatTriangle(): degenerated triangle! triangle inequality does not hold!");
    }
    return std::sqrt( s * (s-l1) * (s-l2) * (s-l3) );
    
}

/**
 * \brief Compute area of geodesic triangle defined by three points in shell space
 * \author Heeren
 *
 * Subdivide the geodesic triangle \f$T\f$ into a set all smaller triangles \f$t_1, ..., t_n\f$ and
 * set \f$ \mathrm{area}(T) = \sum_i |t_i|\f$, where \f$|t_i|\f$ denotes the area of the corresponding flat triangle.
 */
template<typename ConfiguratorType>
double getAreaOfGeodesicTriangle( const MeshTopologySaver& Topology, 
                                const typename ConfiguratorType::VectorType& point0,
                                const typename ConfiguratorType::VectorType& point1,
                                const typename ConfiguratorType::VectorType& point2,
                                const int outerK,
                                const int innerK,
                                const std::vector<int>& Mask,
                                const DeformationBase<ConfiguratorType>& W,
                                const OptimizationParameters<ConfiguratorType>& optPars,
                                bool quiet = true ){
    if( outerK < 2 )
        throw BasicException("getAreaOfGeodesicTriangle(): outer K too smal, should be at least 2!");
    
    typedef typename ConfiguratorType::VectorType VectorType;
    
    if(!quiet) std::cerr << "===========================================================" << std::endl;
    if(!quiet) std::cerr << "Start to compute area of GEODESIC triangle." << std::endl;
    if(!quiet) std::cerr << "===========================================================" << std::endl;
    
    // compute geodesics along edges
    VectorType edgePath1, edgePath2;
    if(!quiet) std::cerr << "Compute geodesic along first edge..." << std::endl;
    integerInterpolation( Topology, point0, point1, Mask, W, optPars, outerK+1, edgePath1, quiet );
    if(!quiet) std::cerr << "Compute geodesic along second edge..." << std::endl;
    integerInterpolation( Topology, point0, point2, Mask, W, optPars, outerK+1, edgePath2, quiet );
    
    // bring into more convenient form
    std::vector< Eigen::Ref<const VectorType> > edgePoints1, edgePoints2;
    const int numLocalDofs = 3 * Topology.getNumVertices();  
    edgePoints1.reserve(outerK);
    edgePoints2.reserve(outerK); 
    for( int k = 0; k < outerK-1; k++ ){
       edgePoints1.push_back( edgePath1.segment(k*numLocalDofs, numLocalDofs) );
       edgePoints2.push_back( edgePath2.segment(k*numLocalDofs, numLocalDofs) );
    }
    edgePoints1.push_back( point1 );
    edgePoints2.push_back( point2 );
    
    // prepare subdivision of triangle
    std::vector<VectorType> upperPoints, lowerPoints;
    upperPoints.push_back( point0 );
    lowerPoints.push_back( edgePoints1[0] );
    lowerPoints.push_back( edgePoints2[0] );
    double fullArea = getAreaOfFlatTriangle<ConfiguratorType>( Topology, point0, edgePoints1[0], edgePoints2[0], innerK, Mask, W, optPars, true );
    
    // start to subdivide triangle
    if(!quiet) std::cerr << "Start to subdivide triangle." << std::endl;
    int numSubTriangs = 1;
    for( int k = 1; k < outerK; k++ ){
        
        // update upper points
        upperPoints.resize( lowerPoints.size() );
        upperPoints = lowerPoints;  
        lowerPoints.resize(0);
        
        // compute k-geodesic from left to right
        VectorType tempPath;
        if(!quiet) std::cerr << "Compute geodesic of length " << k + 2 << std::endl;
        integerInterpolation( Topology, edgePoints1[k], edgePoints2[k], Mask, W, optPars, k+2, tempPath, true );
        
        // update lower points
        lowerPoints.push_back( edgePoints1[k] );
        for( int i = 0; i < k; i++ )
          lowerPoints.push_back( tempPath.segment(i*numLocalDofs, numLocalDofs) );
        lowerPoints.push_back( edgePoints2[k] );
        
        // add area of 2k+1 subtriangles
        std::cerr << "Add area of " << 2*k+1 << " subtriangles. " << std::endl;
        fullArea += getAreaOfFlatTriangle<ConfiguratorType>( Topology, lowerPoints[0], lowerPoints[1], upperPoints[0],  innerK, Mask, W, optPars, true );
        numSubTriangs += 2*k+1;
        for( int i = 0; i < k; i++ ){
            fullArea += getAreaOfFlatTriangle<ConfiguratorType>( Topology, lowerPoints[i+1], upperPoints[i], upperPoints[i+1], innerK, Mask, W, optPars, true );
            fullArea += getAreaOfFlatTriangle<ConfiguratorType>( Topology, lowerPoints[i+1], lowerPoints[i+2], upperPoints[i+1], innerK, Mask, W, optPars, true );
        }        
    }    
    return fullArea;    
    
}

/**
 * \brief Shoot two discrete geodesic rays from foot point and compute distances between points on rays
 * \author Heeren
 *
 * For foot point \f$y\f$ and two points \f$p_K\f$ and \f$q_K\f$, compute discrete geodesics \f$(p_0, p_1, ..., p_K)\f$ and \f$(q_0, q_1, ..., q_K)\f$
 * where \f$p_0 = q_0 = y\f$. Then compute lengths of geodesics connecting \f$p_k\f$ and \f$q_k\f$ for \f$k = 1, ..., K\f$.
 */
template<typename ConfiguratorType>
void shootTwoGeodesicsFromPoint( const MeshTopologySaver& Topology,
                                   const typename ConfiguratorType::VectorType& point0,
                                   const typename ConfiguratorType::VectorType& point1,
                                   const typename ConfiguratorType::VectorType& point2,
                                   const int K,
                                   const std::vector<int>& Mask,
                                   const DeformationBase<ConfiguratorType>& W,
                                   const OptimizationParameters<ConfiguratorType>& optPars,
                                   bool quiet = true ){
    typedef typename ConfiguratorType::VectorType VectorType;
    typename ConfiguratorType::RealType energy;
    VectorType edgePath1, edgePath2;
    
    // compute edges and edge lengths    
    if(!quiet) std::cerr << "Compute geodesic of length " << K+1 << " from p0 to p1..." << std::endl;
    integerInterpolation( Topology, point0, point1, Mask, W, optPars, K+1, edgePath1, true );
    if(!quiet) std::cerr << "Compute geodesic of length " << K+1 << " from p0 to p2..." << std::endl;
    integerInterpolation( Topology, point0, point2, Mask, W, optPars, K+1, edgePath2, true );
    
     // bring into more convenient form
    std::vector< Eigen::Ref<const VectorType> > shapes1, shapes2;
    shapes1.reserve(K);
    shapes2.reserve(K);
    const int numLocalDofs = 3 * Topology.getNumVertices();    
    for( int k = 0; k < K-1; k++ ){
       shapes1.push_back( edgePath1.segment(k*numLocalDofs, numLocalDofs) );
       shapes2.push_back( edgePath2.segment(k*numLocalDofs, numLocalDofs) );
    }
    shapes1.push_back( point1 );
    shapes2.push_back( point2 );
    
    // now compute connections between two outgoing geodesics
    VectorType lengths(K);
    for( int i = 0; i < K; i++ ){
        VectorType tempPath;
        if(!quiet) std::cerr << "Compute geodesic of length " << 2*i+3 << "..." << std::endl;
        integerInterpolation( Topology, shapes1[i], shapes2[i], Mask, W, optPars, 2*i+3, tempPath, true );
        DiscretePathEnergy<ConfiguratorType>( W, 2*i+2, shapes1[i], shapes2[i] ).apply( tempPath, energy );
        lengths[i] = std::sqrt( energy[0] );
        std::cerr << "Length of " << i+1 << "th connection is " << lengths[i] << std::endl;
        if( i > 0 )
            std::cerr << "slope = " << lengths[i] - lengths[i-1] << std::endl;
    }
    
    std::cerr << "Print all slopes: ";
    for( int i = 1; i < K; i++ )
        std::cerr << lengths[i] - lengths[i-1] << " ";
    std::cerr << std::endl;
    
}

#endif

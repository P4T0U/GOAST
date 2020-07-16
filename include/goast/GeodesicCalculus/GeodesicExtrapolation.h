// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Energies and optimization operator for discrete geodesic extrapolation.
 * \author Heeren
 * 
 */

#ifndef GEODESICEXTRAPOLATION_HH
#define GEODESICEXTRAPOLATION_HH

#include <iostream>

#include <goast/Core.h>
#include <goast/Optimization.h>

#include <goast/NRIC/LeastSquaresReconstruction.h>
#include "ConvexGeodesicCombination.h"

//===============================================================================================================================
//===============================================================================================================================

/**
 * \brief Vector-valued discrete Exp2 operator.
 * \author Heeren
 *
 * Given two shapes \f$ S_0 \f$ and \f$ S_1 \f$, this class realizes
 * the vector-valued functional \f$ F[S] = W_{,2}[S_0, S_1] + W_{,1}[S_1, S] \f$,
 * where \f$ W \f$ is a deformation energy passed in the constructor
 * Here \f$ W_{,i}\f$ denotes the derivative w.r.t. the ith argument of \f$ W[.,.] \f$.
 *
 * Note that \f$ (S_0, S_1, S_2) \f$ is a 3-point-geodesic iff. \f$ F[S_2] = 0 \f$.
 */
template <typename ConfiguratorType >
class Exp2Energy : public BaseOp< typename ConfiguratorType::VectorType > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;

  const DeformationBase<ConfiguratorType>& _W;
  const VectorType& _shape0;
  const VectorType& _shape1;
  VectorType _constPart;
  int _numDofs;

public:
    Exp2Energy ( const DeformationBase<ConfiguratorType>& W,
                 const VectorType& shape0,
                 const VectorType& shape1 ) :
      _W(W),
      _shape0(shape0),
      _shape1(shape1),
      _numDofs(shape0.size()){
        // initialization
        calcConstPartOfEnergy();
      }


  //! The vertex positions of S_2 are given as argument.
  void apply( const VectorType& shape2, VectorType& Dest ) const {

    if( shape2.size() != _numDofs )
      throw BasicException ( "Exp2Energy::apply(): arg has wrong size!" );
    if( Dest.size() != _numDofs )
      Dest.resize( _numDofs );
    
    // add constant partial
    Dest = _constPart;
    // add other part
    _W.applyAddUndefGradient( _shape1, shape2, Dest );

  }

protected:
  // pre-cimpute constant part of energy W[S_0, S_1]
  void calcConstPartOfEnergy() {
    _constPart.resize( _numDofs );
    _W.applyDefGradient( _shape0, _shape1, _constPart );
  }

};

//! \brief Matrix-valued derivative of discrete Exp2 operator.
//! \author Heeren
template <typename ConfiguratorType >
class Exp2Gradient : public BaseOp< typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType       TripletType;
  typedef std::vector<TripletType> TripletListType;

  const DeformationBase<ConfiguratorType>& _W;
  const VectorType& _shape1;
  int _numDofs;

public:
    Exp2Gradient ( const DeformationBase<ConfiguratorType>& W,
                 const VectorType& shape0,
                 const VectorType& shape1 ) :
      _W(W),
      _shape1(shape1),
      _numDofs(shape0.size()){  }
      
  //! The vertex positions of S_2 are given as argument.
  void apply( const VectorType& shape2, MatrixType& Dest ) const {

    if( shape2.size() != _numDofs )
      throw BasicException ( "Exp2Gradient::apply(): arg has wrong size!" );
    if( (Dest.rows() != _numDofs) || (Dest.cols() != _numDofs) )
      Dest.resize( _numDofs, _numDofs );

    // fill triplet lists
    TripletListType tripletList;
    tripletList.reserve( _W.numOfNonZeroHessianEntries() );      
         
    pushTriplets( shape2, tripletList );
    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
  
  void pushTriplets( const VectorType& Arg, TripletListType& tripletList ) const {
    if( Arg.size() != _numDofs )
      throw BasicException ( "Exp2Gradient::pushTriplets(): arg has wrong size!" );
    _W.pushTripletsMixedHessian( _shape1, Arg, tripletList, 0, 0, false, 1. );
  }


};

//==============================================================================================================================
//==============================================================================================================================

/**
 * \brief Vector-valued INVERSE/BACKWARD discrete Exp2 operator.
 * \author Heeren
 *
 * Given two shapes \f$ S_1 \f$ and \f$ S_2 \f$ this class realizes
 * the vector-valued functional \f$ F[S] = W_{,2}[S, S_1] +  W_{,1}[S_1, S_2] \f$
 * needed for the inverse/backward Exp2 operator.
 *
 * That means \f$ F[S] = 0 \f$ iff.  \f$ (S, S_1, S_2) \f$ is a discrete geodesic
 */
template <typename ConfiguratorType>
class Exp2InvEnergy : public BaseOp<typename ConfiguratorType::VectorType>{

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;

  const DeformationBase<ConfiguratorType>& _W;
  const VectorType& _shape1, _shape2;
  VectorType _constPart;
  int _numDofs;

public:
    Exp2InvEnergy ( const DeformationBase<ConfiguratorType>& W,
                 const VectorType& shape1,
                 const VectorType& shape2 ) :
      _W(W),
      _shape1(shape1),
      _shape2(shape2),
      _numDofs(shape1.size()){
        // initialization
        calcConstPartOfEnergy();
      }


  //! The vertex positions of S_0 are given as argument.
  void apply( const VectorType& shape0, VectorType& Dest ) const {

    if( shape0.size() != _numDofs )
      throw BasicException ( "Exp2InvEnergy::apply(): arg has wrong size!" );
    if( Dest.size() != _numDofs )
      Dest.resize( _numDofs );
    
    // add constant partial
    Dest = _constPart;
    // add other part
    _W.applyAddDefGradient( shape0, _shape1, Dest );

  }

protected:
  // pre-cimpute constant part of energy W[S_0, S_1]
  void calcConstPartOfEnergy() {
    _constPart.resize( _numDofs );
    _W.applyUndefGradient( _shape1, _shape2, _constPart );
  }
};


//! \brief Matrix-valued derivative of INVERSE/BACKWARD discrete Exp2 operator.
//! \author Heeren
template <typename ConfiguratorType >
class Exp2InvGradient : public BaseOp< typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType       TripletType;
  typedef std::vector<TripletType> TripletListType;

  const DeformationBase<ConfiguratorType>& _W;
  const VectorType& _shape1;
  int _numDofs;

public:
    Exp2InvGradient ( const DeformationBase<ConfiguratorType>& W,
                 const VectorType& shape1,
                 const VectorType& shape2 ) :
      _W(W),
      _shape1(shape1),
      _numDofs(shape1.size()){  }
      
  //! The vertex positions of S_2 are given as argument.
  void apply( const VectorType& shape0, MatrixType& Dest ) const {

    if( shape0.size() != _numDofs )
      throw BasicException ( "Exp2Gradient::apply(): arg has wrong size!" );
    if( (Dest.rows() != _numDofs) || (Dest.cols() != _numDofs) )
      Dest.resize( _numDofs, _numDofs );
      
    // fill triplet lists
    TripletListType tripletList;
    tripletList.reserve( _W.numOfNonZeroHessianEntries() );      
    pushTriplets( shape0, tripletList );
    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
  
  void pushTriplets( const VectorType& shape0, TripletListType& tripletList ) const {
    if( shape0.size() != _numDofs )
      throw BasicException ( "Exp2Gradient::pushTriplets(): arg has wrong size!" );
    _W.pushTripletsMixedHessian( shape0, _shape1, tripletList, 0, 0, true, 1. );
  }

};

//===============================================================================================================================
//===============================================================================================================================

/**
 * \brief Discrete geodesic extrapolation / shooting operator
 * \author Heeren
 *
 * Performs (inverse) geodesic shooting by iteratively computing roots of the (inverse) discrete Exp2 operator functional.
 */
template<typename ConfiguratorType>
class GeodesicExtrapolation{
    
protected:    
  typedef typename ConfiguratorType::RealType          RealType;

  typedef typename ConfiguratorType::VectorType        VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;
  
  
  const MeshTopologySaver& _topology;  
  const VectorType& _positionGeom, _variationGeom;
  const DeformationBase<ConfiguratorType>& _W;  
  
  const OptimizationParameters<ConfiguratorType>& _optPars;
  bool _quiet;
  std::vector<int> const* _bdryMask;  
  
public:    
    GeodesicExtrapolation( const MeshTopologySaver& Topology,
                                const VectorType& PositionGeom, 
                                const VectorType& VariationGeom, 
                                const DeformationBase<ConfiguratorType>& W,
                                const OptimizationParameters<ConfiguratorType>& optPars,
                                bool quiet = true ) 
    :  _topology(Topology), _positionGeom(PositionGeom), _variationGeom(VariationGeom), _W(W), _optPars(optPars), _quiet(quiet), _bdryMask(NULL){}
    
  // set boundary mask    
  void setBoundaryMask( const std::vector<int>& Mask ){
    _bdryMask = &Mask;  
  }      
    
  //! Given p as "position" and q as "variational shape" compute EXP^K_p( v ), 
  //! where v := (q - p)  and K is the number of shooting steps
  void execute( int shootingSteps, VectorType& extrapolatedPath, bool forwardShooting = true ) const { 
      std::vector<VectorType> extrapolatedPathAlias;
      execute( shootingSteps, extrapolatedPathAlias, forwardShooting );      
      // write back 
      int numLocalDofs = 3 * _topology.getNumVertices();
      extrapolatedPath.resize( shootingSteps * numLocalDofs );
      for( int i = 0; i < shootingSteps; i++ )
          extrapolatedPath.segment( i * numLocalDofs, numLocalDofs ) = extrapolatedPathAlias[2+i];
  }
    
  //! Given p as "position" and q as "variational shape" compute EXP^K_p( v ), 
  //! where v := (q - p)  and K is the number of shooting steps
  void execute( int shootingSteps, std::vector<VectorType>& extrapolatedPath, bool forwardShooting = true ) const {    

    // resize path
    extrapolatedPath.resize(shootingSteps+2);
    extrapolatedPath[0] = forwardShooting ? _positionGeom : _variationGeom;
    extrapolatedPath[1] = forwardShooting ? _variationGeom : _positionGeom;
    
    if(!_quiet){
        std::cerr << "=================================================================" << std::endl;
        if(forwardShooting)
            std::cerr << "Start forward shooting with " << shootingSteps << " steps." << std::endl;
        else
            std::cerr << "Start backward shooting with " << shootingSteps << " steps." << std::endl;
        if(_bdryMask)
            std::cerr << "Use " << _bdryMask->size()/3 << " Dirichlet boundary nodes" << std::endl;
        std::cerr << "=================================================================" << std::endl;
    }
    
    // shooting 
    for( int k = 0; k < shootingSteps; k++ ){

      if( !_quiet ) std::cerr << "Step " << k+1 << " of " << shootingSteps << " steps..." << std::endl;
      
      // resize
      extrapolatedPath[k+2].resize( _positionGeom.size() );
      
      if( forwardShooting ){
        typedef Exp2Energy<ConfiguratorType>   GradType;
        typedef Exp2Gradient<ConfiguratorType> HessType;
        GradType DF( _W, extrapolatedPath[k], extrapolatedPath[k+1] );
        HessType D2F( _W, extrapolatedPath[k], extrapolatedPath[k+1] );   
        //solveLeastSquaresLineSearchNewtonSingleShape<ConfiguratorType, GradType, HessType >( _topology, DF, D2F, extrapolatedPath[k+1], extrapolatedPath[k+2], _optPars, _bdryMask );
        //solveLeastSquaresGradientDescentSingleShape<ConfiguratorType, GradType, HessType >( _topology, DF, D2F, extrapolatedPath[k+1], extrapolatedPath[k+2], _optPars, _bdryMask );
        solveNewtonSingleShape<ConfiguratorType, GradType, HessType >( _topology, DF, D2F, extrapolatedPath[k+1], extrapolatedPath[k+2], _optPars, _bdryMask );        
      }
      else{
        typedef Exp2InvEnergy<ConfiguratorType>   GradType;
        typedef Exp2InvGradient<ConfiguratorType> HessType;
        GradType DF( _W, extrapolatedPath[k+1], extrapolatedPath[k] );
        HessType D2F( _W, extrapolatedPath[k+1], extrapolatedPath[k] );   
        //solveLeastSquaresLineSearchNewtonSingleShape<ConfiguratorType, GradType, HessType >( _topology, DF, D2F, extrapolatedPath[k+1], extrapolatedPath[k+2], _optPars, _bdryMask );
        //solveLeastSquaresGradientDescentSingleShape<ConfiguratorType, GradType, HessType >( _topology, DF, D2F, extrapolatedPath[k+1], extrapolatedPath[k+2], _optPars, _bdryMask );
        solveNewtonSingleShape<ConfiguratorType, GradType, HessType >( _topology, DF, D2F, extrapolatedPath[k+1], extrapolatedPath[k+2], _optPars, _bdryMask );
      }
    }    

  }
  
  void execute( int shootingSteps, bool forwardShooting, std::string filename ) const {
      // shoot
      std::vector<VectorType> extrapolatedPath;
      execute( shootingSteps, extrapolatedPath, forwardShooting );
      //save
      saveSolution( extrapolatedPath, filename );      
  }
  
  //
  void saveSolution( const std::vector<VectorType>& extrapolatedPath, std::string saveNameStem ) const {      
      TriMesh auxMesh( _topology.getGrid() );
      for( uint k = 0; k < extrapolatedPath.size(); k++ ){
        std::ostringstream saveName;
        saveName << saveNameStem << "_" << k << ".ply";
        setGeometry( auxMesh, extrapolatedPath[k] );
        if ( !OpenMesh::IO::write_mesh( auxMesh, saveName.str() ) )
          throw BasicException("GeodesicExtrapolation::saveSolution: could not write mesh!");
      }
  }

};

#endif

// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Nonlinear least-squares reconstruction of nodal positions from lengths and angles using quadratic energy
 * \author Heeren
 *
 * Based on Fröhlich, S., & Botsch, M. (2011). Example-Driven Deformations Based on Discrete Shells. Computer Graphics
 * Forum, 30(8), 2246–2257.
 */


#ifndef NRIC_LEASTSQUARESRECONSTRUCTION_H
#define NRIC_LEASTSQUARESRECONSTRUCTION_H

#include <chrono>
#include <ctime>
#include <goast/Core.h>
#include <goast/Optimization.h>

//======================================================================================================================================
//======================================================================================================================================
/**
 * \brief Compute the weights of the quadratic discrete shells energy
 * \author Heeren
 * \param Topology topology of the mesh
 * \param Geometry nodal positions
 * \param Weights computed weights
 * \param squaredValues Whether the weights should be squared or not
 *
 * Computed as in Froehlich and Botsch, 2011:
 * edge length weight is for edge e is 1/l_e
 * dihedral angle weight is for edge e is l_e/sqrt(A_e)
 */
template<typename ConfiguratorType>
void getReconstructionWeights( const MeshTopologySaver &Topology,
                               const typename ConfiguratorType::VectorType &Geometry,
                               typename ConfiguratorType::VectorType &Weights,
                               bool squaredValues = false ) {

  if ( Geometry.size() != 3 * Topology.getNumVertices())
    throw BasicException( "getReconstructionWeights(): input geometry has wrong size!!!" );

  Weights.resize( 2 * Topology.getNumEdges());
  Weights.setZero();

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::RealType RealType;

  for ( int edgeIdx = 0; edgeIdx < Topology.getNumEdges(); ++edgeIdx ) {

    int pi( Topology.getAdjacentNodeOfEdge( edgeIdx, 0 )),
            pj( Topology.getAdjacentNodeOfEdge( edgeIdx, 1 )),
            pk( Topology.getOppositeNodeOfEdge( edgeIdx, 0 )),
            pl( Topology.getOppositeNodeOfEdge( edgeIdx, 1 ));

    // set up vertices
    VecType Pi, Pj, P, temp;
    getXYZCoord<VectorType, VecType>( Geometry, Pi, pi );
    getXYZCoord<VectorType, VecType>( Geometry, Pj, pj );
    RealType edgeLengthSqr = dotProduct( Pj - Pi, Pj - Pi );
    Weights[edgeIdx] = squaredValues ? 1. / edgeLengthSqr : 1. / std::sqrt( edgeLengthSqr );

    if ( std::min( pl, pk ) < 0 )
      continue;

    // compute volume A = (|T_1| + |T_2|) / 3,
    // where T_1 and T_2 share the edge e
    RealType vol = 0.;
    getXYZCoord<VectorType, VecType>( Geometry, P, pk );
    temp.makeCrossProduct( P - Pj, Pi - P );
    vol += temp.norm() / 6.;
    getXYZCoord<VectorType, VecType>( Geometry, P, pl );
    temp.makeCrossProduct( P - Pi, Pj - P );
    vol += temp.norm() / 6.;

    Weights[Topology.getNumEdges() + edgeIdx] = squaredValues ? edgeLengthSqr / vol : std::sqrt( edgeLengthSqr / vol );
  }

}


//======================================================================================================================================
//======================================================================================================================================


/**
 * \brief Sets all angle weights below a user specified threshold to zero
 * \author Heeren
 * \param Angles considered dihedral angles
 * \param AngleWeights considered weights for the different angles
 * \param threshold threshold on the size of the angle
 * \param quiet enable or disable debug output
 * \return Number of changed weights
 *
 * \todo Is this correct?!
 */
template<typename ConfiguratorType>
int correctWeights( const typename ConfiguratorType::VectorType &Angles,
                    typename ConfiguratorType::VectorType &AngleWeights, double threshold, bool quiet = true ) {
  int numOfCorrAngles = 0;
  if ( Angles.size() != AngleWeights.size())
    throw BasicException( "correctWeights(): sizes do not match!!!" );
  for ( int i = 0; i < Angles.size(); i++ ) {
    if ( std::abs( Angles[i] ) > threshold ) {
      AngleWeights[i] = 0.;
      if ( !quiet ) std::cerr << i << "th dihedral a. is " << Angles[i] << std::endl;
      numOfCorrAngles++;
    }
  }
  return numOfCorrAngles;
}

//======================================================================================================================================
//======================================================================================================================================
//======================================================================================================================================

/**
 * \brief Weighted residual between given target lengths and angles and ones given by nodal positions
 * \author Heeren
 *
 * Given a mesh with n vertices and m edges. For prescribed edge lengths L = (L_e)_e and dihedral angles D = (D_e)_e
 * the function maps nodal positions (x,y,z) \in \R^{3n} to the two m-vectors ( sqrt(alpha_e) * (l_e[x,y,z] - L_e ) )_e
 * and ( sqrt(beta_e) * (d_e(x,y,z) - D_e) )_e, i.e. f: R^{3n} -> R^{2m}.
 */
template<typename ConfiguratorType>
class LinearReconstructionFunctional : public BaseOp<typename ConfiguratorType::VectorType> {

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::RealType RealType;

  const MeshTopologySaver _topology;
  const VectorType &_targetLengthsAngles, _weights;
  RealType _edgeLengthWeight, _bendingWeight;
  bool _useBending;

public:
  /**
   * \param Topology topology of the underlying mesh
   * \param prescribedLengthsAngles target lengths and angles to reconstruct
   * \param Weights weights of the quadratic energy in lengths and angles, i.e. per entry weights for the output
   * \param edgeLengthWeight additional weight for the edge lengths (also called membrane weight)
   * \param bendingWeight  additional weight for the dihedral angles (also called bending weight)
   */
  LinearReconstructionFunctional( const MeshTopologySaver &Topology,
                                  const VectorType &prescribedLengthsAngles,
                                  const VectorType &Weights,
                                  RealType edgeLengthWeight,
                                  RealType bendingWeight )
          : _topology( Topology ),
            _targetLengthsAngles( prescribedLengthsAngles ),
            _weights( Weights ),
            _edgeLengthWeight( std::sqrt( edgeLengthWeight )),
            _bendingWeight( std::sqrt( bendingWeight )),
            _useBending( _bendingWeight > 1.e-12 ) {
    if ( prescribedLengthsAngles.size() != 2 * Topology.getNumEdges())
      throw BasicException(
              "LinearReconstructionFunctional::LinearReconstructionFunctional(): wrong size of prescribed arguments!!!" );
    if ( Weights.size() != 2 * Topology.getNumEdges())
      throw BasicException(
              "LinearReconstructionFunctional::LinearReconstructionFunctional(): wrong size of weights!!!" );
  }

  /**
   * \brief Update bending weight used in the energy
   * \param mu new bending weight
   */
  void setBendingWeight( RealType mu ) {
    _bendingWeight = mu;
    _useBending = _bendingWeight > 1.e-12;
  }

  /**
   * \param Arg nodal positions
   * \param Dest weighted lengths and angles
   */
  void apply( const VectorType &Arg, VectorType &Dest ) const {

    if ( Arg.size() != 3 * _topology.getNumVertices())
      throw BasicException( "LinearReconstructionFunctional::apply(): wrong size of arguments in Arg!!!" );

    if ( Dest.size() != 2 * _topology.getNumEdges())
      Dest.resize( 2 * _topology.getNumEdges());
    Dest.setZero();

    // run over all edges
    for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ) {

      int pi( _topology.getAdjacentNodeOfEdge( edgeIdx, 0 )),
              pj( _topology.getAdjacentNodeOfEdge( edgeIdx, 1 )),
              pk( _topology.getOppositeNodeOfEdge( edgeIdx, 0 )),
              pl( _topology.getOppositeNodeOfEdge( edgeIdx, 1 ));

      // set up vertices
      VecType Pi, Pj;
      getXYZCoord<VectorType, VecType>( Arg, Pi, pi );
      getXYZCoord<VectorType, VecType>( Arg, Pj, pj );

      // compute length of edge
      RealType edgeLength = std::sqrt( dotProduct( Pj - Pi, Pj - Pi ));
      Dest[edgeIdx] += _edgeLengthWeight * _weights[edgeIdx] * ( edgeLength - _targetLengthsAngles[edgeIdx] );

      // no bending at boundary edges
      if ( std::min( pl, pk ) < 0 )
        continue;

      // add bending contribution at all?
      if ( !_useBending )
        continue;

      VecType Pk, Pl;
      getXYZCoord<VectorType, VecType>( Arg, Pk, pk );
      getXYZCoord<VectorType, VecType>( Arg, Pl, pl );

      // compute dihedral angle
      RealType dihedralAngle = getDihedralAngle( Pi, Pj, Pk, Pl );
      Dest[_topology.getNumEdges() + edgeIdx] += _bendingWeight * _weights[_topology.getNumEdges() + edgeIdx] *
                                                 ( dihedralAngle -
                                                   _targetLengthsAngles[_topology.getNumEdges() + edgeIdx] );
    }

  }

};

/**
 * \brief Derivative of the weighted residual between given target lengths and angles and ones given by nodal positions
 * \author Heeren
 *
 * \sa LinearReconstructionFunctional
 */
template<typename ConfiguratorType>
class LinearReconstructionDerivative
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::VecType VecType;

  typedef std::vector<TripletType> TripletListType;


  const MeshTopologySaver _topology;
  const VectorType &_targetLengthsAngles, _weights;
  RealType _edgeLengthWeight, _bendingWeight;
  bool _useBending;

public:
  /**
   * \param Topology topology of the underlying mesh
   * \param prescribedLengthsAngles target lengths and angles to reconstruct
   * \param Weights weights of the quadratic energy in lengths and angles, i.e. per entry weights for the output
   * \param edgeLengthWeight additional weight for the edge lengths (also called membrane weight)
   * \param bendingWeight  additional weight for the dihedral angles (also called bending weight)
   */
  LinearReconstructionDerivative( const MeshTopologySaver &Topology,
                                  const VectorType &prescribedLengthsAngles,
                                  const VectorType &Weights,
                                  RealType edgeLengthWeight,
                                  RealType bendingWeight )
          : _topology( Topology ),
            _targetLengthsAngles( prescribedLengthsAngles ),
            _weights( Weights ),
            _edgeLengthWeight( std::sqrt( edgeLengthWeight )),
            _bendingWeight( std::sqrt( bendingWeight )),
            _useBending( _bendingWeight > 1.e-12 ) {
    if ( prescribedLengthsAngles.size() != 2 * Topology.getNumEdges())
      throw BasicException(
              "LinearReconstructionDerivative::LinearReconstructionDerivative(): wrong size of prescribed arguments!!!" );
    if ( Weights.size() != 2 * Topology.getNumEdges())
      throw BasicException(
              "LinearReconstructionDerivative::LinearReconstructionDerivative(): wrong size of weights!!!" );
  }

  /**
   * \brief Update bending weight used in the energy
   * \param mu new bending weight
   */
  void setBendingWeight( RealType mu ) {
    _bendingWeight = mu;
    _useBending = _bendingWeight > 1.e-12;
  }

  void apply( const VectorType &Arg, MatrixType &Dest ) const {

    if ( Arg.size() != 3 * _topology.getNumVertices())
      throw BasicException( "LinearReconstructionDerivative::apply(): wrong size of arguments in Arg!!!" );

    if (( Dest.rows() != 2 * _topology.getNumEdges()) || ( Dest.cols() != 3 * _topology.getNumVertices()))
      Dest.resize( 2 * _topology.getNumEdges(), 3 * _topology.getNumVertices());
    Dest.setZero();

    // set up triplet list
    TripletListType tripletList;
    pushTriplets( Arg, tripletList, false );
    // fill matrix from triplets
    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend());

  }

  void applyTransposed( const VectorType &Arg, MatrixType &Dest ) const {

    if ( Arg.size() != 3 * _topology.getNumVertices())
      throw BasicException( "LinearReconstructionDerivative::applyTransposed(): wrong size of arguments in Arg!!!" );

    if (( Dest.rows() != 3 * _topology.getNumVertices()) || ( Dest.cols() != 2 * _topology.getNumEdges()))
      Dest.resize( 3 * _topology.getNumVertices(), 2 * _topology.getNumEdges());
    Dest.setZero();

    // set up triplet list
    TripletListType tripletList;
    pushTriplets( Arg, tripletList, true );
    // fill matrix from triplets
    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend());
  }

  void pushTriplets( const VectorType &Arg, TripletListType &tripletList, bool transposed = false ) const {

    tripletList.clear();
    tripletList.reserve( 3 * ( 2 + 4 ) * _topology.getNumEdges());
    int colOffset = _topology.getNumVertices();
    int rowOffset = _topology.getNumEdges();

    // run over all edges
    for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ) {

      int pi( _topology.getAdjacentNodeOfEdge( edgeIdx, 0 )),
              pj( _topology.getAdjacentNodeOfEdge( edgeIdx, 1 )),
              pk( _topology.getOppositeNodeOfEdge( edgeIdx, 0 )),
              pl( _topology.getOppositeNodeOfEdge( edgeIdx, 1 ));

      // set up vertices
      VecType Pi, Pj;
      getXYZCoord<VectorType, VecType>( Arg, Pi, pi );
      getXYZCoord<VectorType, VecType>( Arg, Pj, pj );
      VecType edge( Pi - Pj );
      RealType edgeLength = edge.norm();

      // assemble in global matrix
      for ( int i = 0; i < 3; i++ ) {
        pushTriplet( tripletList, edgeIdx, i * colOffset + pi,
                     _edgeLengthWeight * _weights[edgeIdx] * edge[i] / edgeLength, transposed );
        pushTriplet( tripletList, edgeIdx, i * colOffset + pj,
                     -1. * _edgeLengthWeight * _weights[edgeIdx] * edge[i] / edgeLength, transposed );
      }

      // no bending at boundary edges
      if ( std::min( pl, pk ) < 0 )
        continue;

      // add bending contribution at all?
      if ( !_useBending )
        continue;

      VecType Pk, Pl;
      getXYZCoord<VectorType, VecType>( Arg, Pk, pk );
      getXYZCoord<VectorType, VecType>( Arg, Pl, pl );

      // compute first derivatives of dihedral angle
      VecType thetak, thetal, thetai, thetaj;
      getThetaGradK( Pi, Pj, Pk, thetak );
      getThetaGradK( Pj, Pi, Pl, thetal );
      getThetaGradI( Pi, Pj, Pk, Pl, thetai );
      getThetaGradJ( Pi, Pj, Pk, Pl, thetaj );

      // assemble in global matrix
      for ( int i = 0; i < 3; i++ ) {
        pushTriplet( tripletList, rowOffset + edgeIdx, i * colOffset + pi,
                     _bendingWeight * _weights[rowOffset + edgeIdx] * thetai[i], transposed );
        pushTriplet( tripletList, rowOffset + edgeIdx, i * colOffset + pj,
                     _bendingWeight * _weights[rowOffset + edgeIdx] * thetaj[i], transposed );
        pushTriplet( tripletList, rowOffset + edgeIdx, i * colOffset + pk,
                     _bendingWeight * _weights[rowOffset + edgeIdx] * thetak[i], transposed );
        pushTriplet( tripletList, rowOffset + edgeIdx, i * colOffset + pl,
                     _bendingWeight * _weights[rowOffset + edgeIdx] * thetal[i], transposed );
      }

    }

  }

protected:
  void pushTriplet( TripletListType &tripletList, int row, int col, RealType value, bool transposed ) const {
    if ( transposed )
      tripletList.emplace_back( col, row, value );
    else
      tripletList.emplace_back( row, col, value );
  }

};

/**
 * \brief Second derivative of the squared(!) weighted residual between target lengths and angles and ones given by nodal positions
 * \author Heeren
 *
 * This implements the full hessian of the squared residual F(x)^T*F(x), where F is the weighted residual,
 * i.e. LinearReconstructionFunctional.
 *
 * \sa LinearReconstructionFunctional
 */
template<typename ConfiguratorType>
class LinearReconstructionHessian
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;

  typedef std::vector<TripletType> TripletListType;


  const MeshTopologySaver _topology;
  const VectorType &_targetLengthsAngles, _weights;
  RealType _edgeLengthWeight, _bendingWeight;
  bool _useBending;

public:
  /**
   * \param Topology topology of the underlying mesh
   * \param prescribedLengthsAngles target lengths and angles to reconstruct
   * \param Weights weights of the quadratic energy in lengths and angles, i.e. per entry weights for the output
   * \param edgeLengthWeight additional weight for the edge lengths (also called membrane weight)
   * \param bendingWeight  additional weight for the dihedral angles (also called bending weight)
   */
  LinearReconstructionHessian( const MeshTopologySaver &Topology,
                               const VectorType &prescribedLengthsAngles,
                               const VectorType &Weights,
                               RealType edgeLengthWeight,
                               RealType bendingWeight )
          : _topology( Topology ),
            _targetLengthsAngles( prescribedLengthsAngles ),
            _weights( Weights ),
            _edgeLengthWeight( std::sqrt( edgeLengthWeight )),
            _bendingWeight( std::sqrt( bendingWeight )),
            _useBending( _bendingWeight > 1.e-12 ) {
    if ( prescribedLengthsAngles.size() != 2 * Topology.getNumEdges())
      throw BasicException( "LinearReconstructionDerivative::LinearReconstructionDerivative(): "
                            "wrong size of prescribed arguments!!!" );
    if ( Weights.size() != 2 * Topology.getNumEdges())
      throw BasicException( "LinearReconstructionDerivative::LinearReconstructionDerivative(): "
                            "wrong size of weights!!!" );
  }

  void setBendingWeight( RealType mu ) {
    _bendingWeight = mu;
    _useBending = _bendingWeight > 1.e-12;
  }

  // Arg contains nodal positions
  void apply( const VectorType &Arg, MatrixType &Dest ) const {

    if ( Arg.size() != 3 * _topology.getNumVertices())
      throw BasicException( "LinearReconstructionDerivative::apply(): wrong size of arguments in Arg!!!" );

    if (( Dest.rows() != 3 * _topology.getNumVertices()) || ( Dest.cols() != 3 * _topology.getNumVertices()))
      Dest.resize( 3 * _topology.getNumVertices(), 3 * _topology.getNumVertices());
    Dest.setZero();

    // set up triplet list
    TripletListType tripletList;
    pushTriplets( Arg, tripletList, false );
    // fill matrix from triplets
    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend());

  }

  void applyTransposed( const VectorType &Arg, MatrixType &Dest ) const {

    if ( Arg.size() != 3 * _topology.getNumVertices())
      throw BasicException( "LinearReconstructionDerivative::applyTransposed(): wrong size of arguments in Arg!!!" );

    if (( Dest.rows() != 3 * _topology.getNumVertices()) || ( Dest.cols() != 3 * _topology.getNumVertices()))
      Dest.resize( 3 * _topology.getNumVertices(), 3 * _topology.getNumVertices());
    Dest.setZero();

    // set up triplet list
    TripletListType tripletList;
    pushTriplets( Arg, tripletList, true );
    // fill matrix from triplets
    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend());
  }

  void pushTriplets( const VectorType &Arg, TripletListType &tripletList, bool transposed = false ) const {
    tripletList.reserve( 10 * 18 * _topology.getNumEdges());

    const int numVertices = _topology.getNumVertices(); // colOffset
    const int numEdges = _topology.getNumEdges(); // rowOffset

    for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ) {
      int pi( _topology.getAdjacentNodeOfEdge( edgeIdx, 0 )),
              pj( _topology.getAdjacentNodeOfEdge( edgeIdx, 1 )),
              pk( _topology.getOppositeNodeOfEdge( edgeIdx, 0 )),
              pl( _topology.getOppositeNodeOfEdge( edgeIdx, 1 ));

      // set up vertices
      VecType Pi, Pj;
      getXYZCoord<VectorType, VecType>( Arg, Pi, pi );
      getXYZCoord<VectorType, VecType>( Arg, Pj, pj );
      VecType edge( Pi - Pj );
      RealType edgeLength = edge.norm();

      RealType factor = 1. * _edgeLengthWeight * _edgeLengthWeight * _weights[edgeIdx] * _weights[edgeIdx];

      MatType tensorProduct;
      tensorProduct.makeTensorProduct( edge, edge );
      tensorProduct *= factor * _targetLengthsAngles[edgeIdx] / ( edgeLength * edgeLength * edgeLength );
      tensorProduct.addToDiagonal( factor * ( edgeLength - _targetLengthsAngles[edgeIdx] ) / edgeLength );

      localToGlobal( tripletList, pi, pi, tensorProduct );
      localToGlobal( tripletList, pj, pj, tensorProduct );
      tensorProduct *= -1.;
      localToGlobal( tripletList, pi, pj, tensorProduct );

      // no bending at boundary edges
      if ( std::min( pl, pk ) == -1 )
        continue;

      VecType Pk, Pl;
      getXYZCoord<VectorType, VecType>( Arg, Pk, pk );
      getXYZCoord<VectorType, VecType>( Arg, Pl, pl );

      // compute first derivatives of dihedral angle
      VecType thetak, thetal, thetai, thetaj;
      getThetaGradK( Pi, Pj, Pk, thetak );
      getThetaGradK( Pj, Pi, Pl, thetal );
      getThetaGradI( Pi, Pj, Pk, Pl, thetai );
      getThetaGradJ( Pi, Pj, Pk, Pl, thetaj );

      factor = 1. * _bendingWeight * _bendingWeight * _weights[numEdges + edgeIdx] * _weights[numEdges + edgeIdx];

      RealType delThetaDouble =
              factor * ( getDihedralAngle( Pi, Pj, Pk, Pl ) - _targetLengthsAngles[numEdges + edgeIdx] );

      // now compute second derivatives of dihedral angle
      MatType H, aux;


      //kk
      getHessThetaKK( Pi, Pj, Pk, aux );
      tensorProduct.makeTensorProduct( thetak, thetak );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H );
      H *= 1;
      localToGlobal( tripletList, pk, pk, H );

      //ik & ki (Hki = Hik)
      getHessThetaIK( Pi, Pj, Pk, aux );
      tensorProduct.makeTensorProduct( thetai, thetak );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H );
      localToGlobal( tripletList, pi, pk, H );

      //jk & kj (Hkj = Hjk)
      getHessThetaJK( Pi, Pj, Pk, aux );
      tensorProduct.makeTensorProduct( thetaj, thetak );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H );
      localToGlobal( tripletList, pj, pk, H );

      //ll
      getHessThetaKK( Pj, Pi, Pl, aux );
      tensorProduct.makeTensorProduct( thetal, thetal );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H );
      H *= 1;
      localToGlobal( tripletList, pl, pl, H );

      //il & li (Hli = Hil)
      getHessThetaJK( Pj, Pi, Pl, aux );
      tensorProduct.makeTensorProduct( thetai, thetal );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H );
      localToGlobal( tripletList, pi, pl, H );

      //jl & lj (Hlj = Hjl)
      getHessThetaIK( Pj, Pi, Pl, aux );
      tensorProduct.makeTensorProduct( thetaj, thetal );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H );
      localToGlobal( tripletList, pj, pl, H );

      //kl/lk: Hkl = 0 and Hlk = 0
      tensorProduct.makeTensorProduct( thetak, thetal );
      tensorProduct *= factor;
      localToGlobal( tripletList, pk, pl, tensorProduct );

      //ii
      getHessThetaII( Pi, Pj, Pk, Pl, aux );
      tensorProduct.makeTensorProduct( thetai, thetai );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H );
      H *= 1;
      localToGlobal( tripletList, pi, pi, H );

      //jj
      getHessThetaII( Pj, Pi, Pl, Pk, aux );
      tensorProduct.makeTensorProduct( thetaj, thetaj );
      getWeightedMatrixSum( factor, tensorProduct, delThetaDouble, aux, H );
      H *= 1;
      localToGlobal( tripletList, pj, pj, H );

      //ij & ji (Hij = Hji)
      getHessThetaJI( Pi, Pj, Pk, Pl, H );
      H *= delThetaDouble;
      tensorProduct.makeTensorProduct( thetai, thetaj );
      H.addMultiple( tensorProduct, factor );
      localToGlobal( tripletList, pi, pj, H );
    }
  }

protected:
  void pushTriplet( TripletListType &tripletList, int row, int col, RealType value, bool transposed = false ) const {
    if ( transposed )
      tripletList.push_back( TripletType( col, row, value ));
    else
      tripletList.push_back( TripletType( row, col, value ));
  }

  void localToGlobal( TripletListType &tripletList, int k, int l, const MatType &localMatrix ) const {
    int numV = _topology.getNumVertices();
    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ )
        tripletList.push_back( TripletType( i * numV + k, j * numV + l, localMatrix( i, j )));

    if ( k != l ) {
      for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
          tripletList.push_back( TripletType( i * numV + l, j * numV + k, localMatrix( j, i )));
    }
  }

};



//======================================================================================================================================
//======================================================================================================================================
//======================================================================================================================================




#endif
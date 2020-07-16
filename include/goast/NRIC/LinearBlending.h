// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Linear blending of lengths and angles for deformation synthesis, e.g. handle editing
 * \author Heeren
 *
 * Based on Fröhlich, S., & Botsch, M. (2011). Example-Driven Deformations Based on Discrete Shells. Computer Graphics
 * Forum, 30(8), 2246–2257.
 */

#ifndef NRIC_LINEARBLENDING_H
#define NRIC_LINEARBLENDING_H

/**
 * \brief Residual energy for deformation synthesis via linear blending of lengths and angles
 * \author Heeren
 *
 * Full residual energy  F: \R^n x \R^n x \R^n x \R^k --> \R^m x \R^m from \cite{FrBo11},
 * including lengths/angles from potentially multiple input shapes
 * Unknowns are given by nodal positions of optimal mesh (i.e. in \R^{3n}) plus k blending weights (in \R^k)
 * Output are the vectors of weighted edge length and dihedral angle residuals, respectively, both living in \R^{m}
 * a) For interpolation or local averaging, the blending weights are supposed to be fixed
 * b) For reconstruction/projecting onto example space, these can be optimized as well
 */
template<typename ConfiguratorType>
class LinearBlendingResidual : public BaseOp<typename ConfiguratorType::VectorType> {
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;

  const MeshTopologySaver _topology;
  const VectorType &_targetLengthsAngles, _weights;
  const RealType _edgeLengthWeight, _bendingWeight;
  int _numData, _numLocalDOFs;
  const std::vector<VectorType> &_inputDataLengthsAngles;

public:
  /**
   * \param Topology Topology of the mesh
   * \param targetLengthsAngles Base point of the affine example space in lengths and angles
   * \param IntegrationWeights Weights of the quadratic energy in lengths and angles
   * \param edgeLengthWeight Weight of the membran, i.e. edge lengths, part of the quadratic energy
   * \param bendingWeight Weight of the bending, i.e. dihedral angle, part of the quadratic energy
   * \param inputDataLengthsAngles Examples shapes in lengths and angles
   */
  LinearBlendingResidual( const MeshTopologySaver &Topology,
                          const VectorType &targetLengthsAngles,
                          const VectorType &IntegrationWeights,
                          RealType edgeLengthWeight,
                          RealType bendingWeight,
                          const std::vector<VectorType> &inputDataLengthsAngles )
          : _topology( Topology ),
            _targetLengthsAngles( targetLengthsAngles ),
            _weights( IntegrationWeights ),
            _edgeLengthWeight( std::sqrt( edgeLengthWeight )),
            _bendingWeight( std::sqrt( bendingWeight )),
            _numData( inputDataLengthsAngles.size()),
            _numLocalDOFs( 3 * Topology.getNumVertices()),
            _inputDataLengthsAngles( inputDataLengthsAngles ) {}

  /**
   * \param Arg Nodal positions and blending weights stacked in one vector
   * \param Dest Lengths and angles residual weighted according to the given weights
   */
  void apply( const VectorType &Arg, VectorType &Dest ) const {

    if ( Arg.size() != _numData + _numLocalDOFs )
      throw BasicException( "LinearBlendingResidual::apply: wrong number of components in Arg!" );

    if ( Dest.size() != 2 * _topology.getNumEdges())
      Dest.resize( 2 * _topology.getNumEdges());
    Dest.setZero();


    int numEdges = _topology.getNumEdges();

    VectorType lengthAngleDifferenceSum = _targetLengthsAngles;
    for ( int i = 0; i < _numData; i++ )
      lengthAngleDifferenceSum += Arg[_numLocalDOFs + i] * ( _inputDataLengthsAngles[i] - _targetLengthsAngles );

    // run over all edges
    for ( int edgeIdx = 0; edgeIdx < numEdges; ++edgeIdx ) {
      int pi( _topology.getAdjacentNodeOfEdge( edgeIdx, 0 )),
              pj( _topology.getAdjacentNodeOfEdge( edgeIdx, 1 )),
              pk( _topology.getOppositeNodeOfEdge( edgeIdx, 0 )),
              pl( _topology.getOppositeNodeOfEdge( edgeIdx, 1 ));

      // set up vertices
      VecType Pi, Pj;
      getXYZCoord<typename VectorType::ConstSegmentReturnType, VecType>( Arg.segment( 0, _numLocalDOFs ), Pi, pi );
      getXYZCoord<typename VectorType::ConstSegmentReturnType, VecType>( Arg.segment( 0, _numLocalDOFs ), Pj, pj );

      // compute squared length of edge
      RealType edgeLength = std::sqrt( dotProduct( Pj - Pi, Pj - Pi ));
      Dest[edgeIdx] = _edgeLengthWeight * _weights[edgeIdx] * ( edgeLength - lengthAngleDifferenceSum[edgeIdx] );

      // no bending at boundary edges
      if ( std::min( pl, pk ) < 0 )
        continue;

      VecType Pk, Pl;
      getXYZCoord<typename VectorType::ConstSegmentReturnType, VecType>( Arg.segment( 0, _numLocalDOFs ), Pk, pk );
      getXYZCoord<typename VectorType::ConstSegmentReturnType, VecType>( Arg.segment( 0, _numLocalDOFs ), Pl, pl );

      // compute dihedral angle
      RealType dihedralAngle = getDihedralAngle( Pi, Pj, Pk, Pl );
      Dest[numEdges + edgeIdx] = _bendingWeight * _weights[numEdges + edgeIdx] *
                                 ( dihedralAngle - lengthAngleDifferenceSum[numEdges + edgeIdx] );
    }
  }

};


/**
* \brief Residual energy for deformation synthesis via linear blending of lengths and angles
* \author Heeren
*
* \sa LinearBlendingResidual
*/
template<typename ConfiguratorType>
class LinearBlendingResidualGradient
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {


  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;
  typedef typename ConfiguratorType::VecType VecType;


  const MeshTopologySaver _topology;
  const VectorType &_targetLengthsAngles, _weights;
  RealType _edgeLengthWeight, _bendingWeight;
  int _numData, _numLocalDOFs;
  const std::vector<VectorType> &_inputDataLengthsAngles;

public:
  /**
   * \param Topology Topology of the mesh
   * \param targetLengthsAngles Base point of the affine example space in lengths and angles
   * \param IntegrationWeights Weights of the quadratic energy in lengths and angles
   * \param edgeLengthWeight Weight of the membran, i.e. edge lengths, part of the quadratic energy
   * \param bendingWeight Weight of the bending, i.e. dihedral angle, part of the quadratic energy
   * \param inputDataLengthsAngles Examples shapes in lengths and angles
   * \param FixBlendingWeights Whether gradient of blending weights should be zero, i.e. they are fixed
   */
  LinearBlendingResidualGradient( const MeshTopologySaver &Topology,
                                  const VectorType &targetLengthsAngles,
                                  const VectorType &IntegrationWeights,
                                  RealType edgeLengthWeight,
                                  RealType bendingWeight,
                                  const std::vector<VectorType> &inputDataLengthsAngles )
          : _topology( Topology ),
            _targetLengthsAngles( targetLengthsAngles ),
            _weights( IntegrationWeights ),
            _edgeLengthWeight( std::sqrt( edgeLengthWeight )),
            _bendingWeight( std::sqrt( bendingWeight )),
            _numData( inputDataLengthsAngles.size()),
            _numLocalDOFs( 3 * Topology.getNumVertices()),
            _inputDataLengthsAngles( inputDataLengthsAngles ) {}

  /**
   * \param Arg Nodal positions and blending weights stacked in one vector
   * \param Dest Triplet list for the gradient matrix
   */
  void pushTriplets( const VectorType &Arg, TripletListType &tripletList ) const {

    if ( Arg.size() != _numData + _numLocalDOFs )
      throw BasicException( "LinearBlendingResidual::pushTriplets(): wrong number of components in Arg!" );

    int numEdges = _topology.getNumEdges();
    int numNodes = _topology.getNumVertices();

    Eigen::Ref<const VectorType> Geometry = Arg.segment( 0, _numLocalDOFs );
    Eigen::Ref<const VectorType> BlendingWeights = Arg.segment( _numLocalDOFs, _numData );

    for ( int edgeIdx = 0; edgeIdx < numEdges; ++edgeIdx ) {

      int pi( _topology.getAdjacentNodeOfEdge( edgeIdx, 0 )),
              pj( _topology.getAdjacentNodeOfEdge( edgeIdx, 1 )),
              pk( _topology.getOppositeNodeOfEdge( edgeIdx, 0 )),
              pl( _topology.getOppositeNodeOfEdge( edgeIdx, 1 ));

      // set up vertices
      VecType Pi, Pj;
      getXYZCoord<Eigen::Ref<const VectorType>, VecType>( Geometry, Pi, pi );
      getXYZCoord<Eigen::Ref<const VectorType>, VecType>( Geometry, Pj, pj );
      VecType edge( Pi - Pj );
      RealType edgeLength = edge.norm();

      // assemble in global matrix
      for ( int i = 0; i < 3; i++ ) {
        tripletList.push_back( TripletType( edgeIdx, pi + i * numNodes,
                                            _edgeLengthWeight * _weights[edgeIdx] * edge[i] / edgeLength ));
        tripletList.push_back( TripletType( edgeIdx, pj + i * numNodes,
                                            -1. * _edgeLengthWeight * _weights[edgeIdx] * edge[i] / edgeLength ));
      }

      // no bending at boundary edges
      if ( std::min( pl, pk ) < 0 )
        continue;

      VecType Pk, Pl;
      getXYZCoord<Eigen::Ref<const VectorType>, VecType>( Geometry, Pk, pk );
      getXYZCoord<Eigen::Ref<const VectorType>, VecType>( Geometry, Pl, pl );

      // compute first derivatives of dihedral angle
      VecType thetak, thetal, thetai, thetaj;
      getThetaGradK( Pi, Pj, Pk, thetak );
      getThetaGradK( Pj, Pi, Pl, thetal );
      getThetaGradI( Pi, Pj, Pk, Pl, thetai );
      getThetaGradJ( Pi, Pj, Pk, Pl, thetaj );

      int offSetEdgeIdx = numEdges + edgeIdx;
      // assemble in global matrix
      for ( int i = 0; i < 3; i++ ) {
        tripletList.push_back(
                TripletType( offSetEdgeIdx, pi + i * numNodes, _bendingWeight * _weights[offSetEdgeIdx] * thetai[i] ));
        tripletList.push_back(
                TripletType( offSetEdgeIdx, pj + i * numNodes, _bendingWeight * _weights[offSetEdgeIdx] * thetaj[i] ));
        tripletList.push_back(
                TripletType( offSetEdgeIdx, pk + i * numNodes, _bendingWeight * _weights[offSetEdgeIdx] * thetak[i] ));
        tripletList.push_back(
                TripletType( offSetEdgeIdx, pl + i * numNodes, _bendingWeight * _weights[offSetEdgeIdx] * thetal[i] ));
      }

      // derivatives wrt. the blending weights
      for ( int i = 0; i < _numData; i++ ) {
        tripletList.push_back( TripletType( edgeIdx, 3 * numNodes + i, -1. * _edgeLengthWeight * _weights[edgeIdx] *
                                                                       ( _inputDataLengthsAngles[i][edgeIdx] -
                                                                         _targetLengthsAngles[edgeIdx] )));
        tripletList.push_back( TripletType( offSetEdgeIdx, 3 * numNodes + i,
                                            -1. * _bendingWeight * _weights[offSetEdgeIdx] *
                                            ( _inputDataLengthsAngles[i][offSetEdgeIdx] -
                                              _targetLengthsAngles[offSetEdgeIdx] )));
      }

    }

  }

  /**
  * \param Arg Nodal positions and blending weights stacked in one vector
  * \param Dest Gradient as sparse matrix
  */
  void apply( const VectorType &Arg, MatrixType &Dest ) const {

    if ( Arg.size() != _numData + _numLocalDOFs )
      throw BasicException( "LinearBlendingResidual::apply: wrong number of components in Arg!" );

    int numEdges = _topology.getNumEdges();
    if ( Dest.rows() != 2 * numEdges || Dest.cols() != _numData + _numLocalDOFs )
      Dest.resize( 2 * numEdges, _numData + _numLocalDOFs );
    Dest.setZero();

    // set up triplet list
    TripletListType tripletList;
    pushTriplets( Arg, tripletList );
    // fill matrix from triplets
    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend());
  }

  /**
   * \param Arg Nodal positions and blending weights stacked in one vector
   * \param Dest Transposed gradient as sparse matrix
   *
   * \todo Add proper transposed assembly via pushTriplets
   */
  void applyTransposed( const VectorType &Arg, MatrixType &Dest ) const {

    if ( Arg.size() != _numData + _numLocalDOFs )
      throw BasicException( "LinearBlendingResidual::applyTransposed: wrong number of components in Arg!" );

    int numEdges = _topology.getNumEdges();
    if ( Dest.rows() != _numData + _numLocalDOFs || Dest.cols() != 2 * numEdges )
      Dest.resize( _numData + _numLocalDOFs, 2 * numEdges );
    Dest.setZero();


    MatrixType temp;
    temp.resize( 2 * numEdges, _numData + _numLocalDOFs );
    temp.setZero();
    // set up triplet list
    TripletListType tripletList;
    pushTriplets( Arg, tripletList );
    // fill matrix from triplets
    temp.setFromTriplets( tripletList.cbegin(), tripletList.cend());

    Dest = temp.transpose();
  }


};

#endif //NRIC_LINEARBLENDING_H

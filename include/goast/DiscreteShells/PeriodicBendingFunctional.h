// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Heuristic bending energy on triangle meshes (for fixed reference domain) and derivatives.
 * \author Heeren
 * \cite EzHeAzRuBe19
 */
 #ifndef PERIODICBENDINGFUNCTIONALS_HH
#define PERIODICBENDINGFUNCTIONALS_HH


//== INCLUDES =================================================================
#include <goast/Core/Auxiliary.h>
#include <goast/Core/LocalMeshGeometry.h>
#include <goast/Core/Topology.h>

/**
 * \brief Periodic version of Discrete Shells bending energy where dihedral angles are replaced by cosine of dihedral angles.
 * \author Heeren
 *
 * Defined in eq. (9) in \cite EzHeAzRuBe19
 *
 * Let \f$ x \f$ be the geometry of a triangle mesh with edge set \f$ E \f$, then this class implements the functional
 * \f[ F[x] = \sum_{e \in E} \alpha_e (t_e - \langle nk_e[x], nl_e[x] \rangle)^2 \f]
 * where \f$ nk_e[x] \f$ and \f$ nl_e[x] \f$ are the two adjacent face normals (to edge e),
 * \f$ \alpha_e =  (l_e)^2 / d_e \f$ are fixed integration weights (the same as in the original Discrete Shell's bending energy),
 * and \f$ t_e = \langle nk_e, nl_e \rangle  \f$ are prescribed angles.
 *
 * All quantities that do not depend on \f$ x \f$ have been precomputed on a reference geometry.
 *
 * Note that \f$ \langle nk_e, nl_e \rangle = \cos(\theta_e) \f$, where \f$ \theta_e \f$ is the dihedral angle at edge \f$e\f$.
 *
 * The functional is refered to as 'periodic', since the energy density is periodic wrt. a rotation of a triangle
 * about the edge to one of its neighboring triangles (due to the cosine!), whereas the classical Discrete Shell's bending energy SimpleBendingFunctional is not.
 */
template<typename ConfiguratorType>
class PeriodicBendingFunctional : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::MatType    MatType;

  const MeshTopologySaver& _topology;
  VectorType  _refCosDihedralAngles, _refSqrEdgeLengths, _refEdgeAreas;
  const RealType _weight;

public:

  PeriodicBendingFunctional( const MeshTopologySaver& topology,
                       const VectorType& ReferenceGeometry,
                       RealType Weight = 1. ) 
  : _topology( topology), 
    _weight(Weight){
        computeReferenceQuantitiesBending<ConfiguratorType>( topology, ReferenceGeometry, _refCosDihedralAngles, _refEdgeAreas, _refSqrEdgeLengths, true );
    }
    
    PeriodicBendingFunctional( const MeshTopologySaver& topology,
                         const VectorType& RefCosDihedralAngles, 
                         const VectorType& RefEdgeAreas, 
                         const VectorType& RefSqrEdgeLengths,
                         RealType Weight = 1. ) : _topology( topology), _refCosDihedralAngles(RefCosDihedralAngles), _refSqrEdgeLengths(RefEdgeAreas), _refEdgeAreas(RefSqrEdgeLengths), _weight(Weight) {}

  // energy evaluation
  void apply( const VectorType& Argument, RealType & Dest ) const {

    if( Argument.size() != 3 * _topology.getNumVertices() ){
      std::cerr << "size of active = " << Argument.size() << " vs. num of dofs = " << 3 * _topology.getNumVertices() << std::endl;
      throw BasicException( "PeriodicBendingFunctional::apply(): wrong input size!");
    }

    Dest = 0.;
    
    for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){
      
      if( !(_topology.isEdgeValid(edgeIdx)) )
	continue;

      int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
          pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
          pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ),
          pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );

      // no bending at boundary edges
      if( std::min( pl, pk) < 0 )
        continue;
      
      // get deformed geometry
      VecType Pi, Pj, Pk, Pl;
      getXYZCoord<VectorType, VecType>( Argument, Pi, pi);
      getXYZCoord<VectorType, VecType>( Argument, Pj, pj);
      getXYZCoord<VectorType, VecType>( Argument, Pk, pk);
      getXYZCoord<VectorType, VecType>( Argument, Pl, pl);

      // compute deformed dihedral angle
      RealType delTheta = getCosOfDihedralAngle( Pi, Pj, Pk, Pl ) - _refCosDihedralAngles[edgeIdx];
      Dest += _weight * delTheta * delTheta * _refSqrEdgeLengths[edgeIdx] / _refEdgeAreas[edgeIdx];

#ifdef DEBUGMODE
      if( std::isnan( Dest ) ){
          std::cerr << "NaN in periodic bending energy in edge " << edgeIdx << "! " << std::endl;
          if( hasNanEntries(Argument) )
            std::cerr << "Argument has NaN entries! " << std::endl;
          else
            std::cerr << "deformed cos(theta) = " << getCosOfDihedralAngle( Pi, Pj, Pk, Pl ) << std::endl;            
          throw BasicException("PeriodicBendingFunctional::apply(): NaN Error!");
      }
#endif
    }
  }

 };

//==========================================================================================================
//! \brief Gradient of PeriodicBendingFunctional
//! \author Heeren
template<typename ConfiguratorType>
class PeriodicBendingGradient : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {

  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::MatType    MatType;

  const MeshTopologySaver& _topology;
  VectorType  _refCosDihedralAngles, _refSqrEdgeLengths, _refEdgeAreas;
  const RealType _weight;
  mutable VecType _projNormal[2];
  mutable RealType _area[2];

public:
PeriodicBendingGradient( const MeshTopologySaver& topology,
                         const VectorType& ReferenceGeometry,
                         RealType Weight = 1. ) : _topology( topology), _weight(Weight) {
                              computeReferenceQuantitiesBending<ConfiguratorType>( topology, ReferenceGeometry, _refCosDihedralAngles, _refEdgeAreas, _refSqrEdgeLengths, true );
                        }
                        
PeriodicBendingGradient( const MeshTopologySaver& topology,
                         const VectorType& RefCosDihedralAngles, 
                         const VectorType& RefEdgeAreas, 
                         const VectorType& RefSqrEdgeLengths,
                         RealType Weight = 1. ) : _topology( topology), _refCosDihedralAngles(RefCosDihedralAngles), _refSqrEdgeLengths(RefEdgeAreas), _refEdgeAreas(RefSqrEdgeLengths), _weight(Weight) {}

void apply( const VectorType& Argument, VectorType& Dest ) const {

  int numDofs = 3 * _topology.getNumVertices();
  if( Argument.size() != numDofs ){
      std::cerr << "size of undef = " << Argument.size() << " vs. num of dofs = " << numDofs << std::endl;
      throw BasicException( "PeriodicBendingGradient::apply(): sizes dont match!");
  }
      
  if( Dest.size() != numDofs )
    Dest.resize( numDofs );  
  Dest.setZero();

  //
  for ( int edgeIdx = 0; edgeIdx < _topology.getNumEdges(); ++edgeIdx ){
    
    if( !(_topology.isEdgeValid(edgeIdx)) )
      continue;

        int pi( _topology.getAdjacentNodeOfEdge(edgeIdx,0) ),
            pj( _topology.getAdjacentNodeOfEdge(edgeIdx,1) ),
            pk( _topology.getOppositeNodeOfEdge(edgeIdx,0) ),
            pl( _topology.getOppositeNodeOfEdge(edgeIdx,1) );

      // no bending at boundary edges
      if( std::min( pl, pk) < 0 )
        continue;

      //! first get defomed quantities
      VecType Pi, Pj, Pk, Pl, temp;
      getXYZCoord<VectorType, VecType>( Argument, Pi, pi);
      getXYZCoord<VectorType, VecType>( Argument, Pj, pj);
      getXYZCoord<VectorType, VecType>( Argument, Pk, pk);
      getXYZCoord<VectorType, VecType>( Argument, Pl, pl);      

      // get normals
      VecType nk, nl, projNk, projNl;
      MatType proj;
      _area[0] = 0.5 * getWeightedNormal( Pi, Pj, Pk, nk, true );
      _area[1] = 0.5 * getWeightedNormal( Pj, Pi, Pl, nl, true );
      getProjection( nk, proj );
      proj.mult( nl, _projNormal[0] );
      getProjection( nl, proj );
      proj.mult( nk, _projNormal[1] );
      
      // compute weighted differnce of dihedral angles
      RealType delTheta = _refCosDihedralAngles[edgeIdx] - (nk*nl);
      delTheta *= -2. * _weight * _refSqrEdgeLengths[edgeIdx] / _refEdgeAreas[edgeIdx];
      
      // D_kN = (2A)^{-1} (Id - NN^T) (Ex)^T, with Ex \in \R^{3,3} s.t. Ex(V) is the cross prodcut between E = Pi - Pj and V, A is area of triangle
      VecType grad;

      // p = k 
      multScaledCrossDiffOp( Pi, Pj, _projNormal[0], 0.5 / _area[0], grad );
      for( int i = 0; i < 3; i++ )
        Dest[i*_topology.getNumVertices() + pk] += delTheta * grad[i];

      // p = l    
      multScaledCrossDiffOp( Pj, Pi, _projNormal[1], 0.5 / _area[1], grad );
      for( int i = 0; i < 3; i++ )
        Dest[i*_topology.getNumVertices() + pl] += delTheta * grad[i];
      
      // p = i 
      multScaledCrossDiffOp( Pj, Pk, _projNormal[0], 0.5 / _area[0], grad );
      for( int i = 0; i < 3; i++ )
        Dest[i*_topology.getNumVertices() + pi] += delTheta * grad[i]; 
      multScaledCrossDiffOp( Pl, Pj, _projNormal[1], 0.5 / _area[1], grad );
      for( int i = 0; i < 3; i++ )
        Dest[i*_topology.getNumVertices() + pi] += delTheta * grad[i];
      
      // p = j
      multScaledCrossDiffOp( Pk, Pi, _projNormal[0], 0.5 / _area[0], grad );
      for( int i = 0; i < 3; i++ )
        Dest[i*_topology.getNumVertices() + pj] += delTheta * grad[i];
      multScaledCrossDiffOp( Pi, Pl, _projNormal[1], 0.5 / _area[1], grad );
      for( int i = 0; i < 3; i++ )
        Dest[i*_topology.getNumVertices() + pj] += delTheta * grad[i];
      
            
#ifdef DEBUGMODE
      if( hasNanEntries( Dest ) ){
          std::cerr << "NaN in simple bending gradient deformed in edge " << edgeIdx << "! " << std::endl;
          if( hasNanEntries(Argument) )
            std::cerr << "Argument has NaN entries! " << std::endl;
          throw BasicException("PeriodicBendingGradient::apply(): NaN Error!");
      }
#endif

  }
}
protected:
  // res = a * C * arg, where C*z = (P0-P1) \times z
  void multScaledCrossDiffOp( const SmallVec3<RealType>& P0, const SmallVec3<RealType>& P1, const SmallVec3<RealType>& arg, RealType Factor, SmallVec3<RealType>& res ) const {
      SmallVec3<RealType> diff(P0 - P1);
      for( int i = 0; i < 3; i++ )
        res[i] = Factor * ( arg[(i+2)%3]*diff[(i+1)%3] - arg[(i+1)%3]*diff[(i+2)%3] );
  }
};

//==========================================================================================================
//==========================================================================================================
//! \brief Hessian of PeriodicBendingFunctional
//! \author Heeren
template<typename ConfiguratorType>
class PeriodicBendingHessian : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {

protected:
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::VecType     VecType;
  typedef typename ConfiguratorType::MatType     MatType;
  
  typedef std::vector<TripletType> TripletListType;
  
  const MeshTopologySaver& _topology;
  const VectorType& _undefShell;
  RealType _factor;
  mutable int _rowOffset, _colOffset;
  
public:
  PeriodicBendingHessian( const MeshTopologySaver& topology,
                             const VectorType& undefShell,
			     const RealType Weight = 1.,
                             int rowOffset = 0, 
                             int colOffset = 0 ) : _topology( topology), _undefShell(undefShell), _factor( Weight ), _rowOffset(rowOffset), _colOffset(colOffset){}
    
  void setRowOffset( int rowOffset ) const {
        _rowOffset = rowOffset;
  }
    
  void setColOffset( int colOffset ) const {
        _colOffset = colOffset;
  }

  //
  void apply( const VectorType& defShell, MatrixType& Dest ) const {    
    assembleHessian( defShell, Dest );
  }

  // assmeble Hessian matrix via triplet list
  void assembleHessian( const VectorType& defShell, MatrixType& Hessian ) const {
    int dofs = 3*_topology.getNumVertices();
    if( (Hessian.rows() != dofs) || (Hessian.cols() != dofs) )
        Hessian.resize( dofs, dofs );
    Hessian.setZero();
    
    // set up triplet list
    TripletListType tripletList;
    // per edge we have 4 active vertices, i.e. 16 combinations each producing a 3x3-matrix
    tripletList.reserve( 16 * 9 * _topology.getNumEdges() );   
    // fill matrix from triplets
    pushTriplets( defShell, tripletList );
    Hessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
  
  // fill triplets
  void pushTriplets( const VectorType& defShell, TripletListType& tripletList ) const {
      throw BasicException("PeriodicBendingHessian has not been implemented yet!");
  }
  
};

#endif
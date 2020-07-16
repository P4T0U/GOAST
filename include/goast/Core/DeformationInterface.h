// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DEFORMATIONINTERFACE_HH
#define DEFORMATIONINTERFACE_HH

//== INCLUDES =================================================================
#include "Auxiliary.h"
#include "LocalMeshGeometry.h"
#include "Topology.h"


//==========================================================================================================
// DEFORMATION ENERGIES
//==========================================================================================================

/**
 * \brief Abstract base class to prescribe structure of deformation energies.
 * \author Heeren
 *
 * Deformation energy  \f$ E \f$  is considered as some functional  \f$ (S_1, S_2) \mapsto E[S_1, S_2] \f$ ,
 * where  \f$ S_1 \f$  and  \f$ S_2 \f$  are refered to as undeformed and deformed configuration, respectively.
 *
 * Several functions are supposed to be provided by derived classes:
 * - applyEnergy()
 * - applyUndefGradient (), i.e. D_1 E[S_1, S_2]
 * - applyDefGradient (), i.e. D_2 E[S_1, S_2]
 * - pushTripletsDefHessian(), i.e. to assemble D_2^2 E[S_1, S_2]
 * - pushTripletsUndefHessian(), i.e. to assemble  D_1^2 E[S_1, S_2]
 * - pushTripletsMixedHessian(), i.e. to assemble D_1 D_2 E[S_1, S_2]
 * where D_i denotes the derivative w.r.t. the ith argument.
 */
template <typename ConfiguratorType>
class DeformationBase {

protected:    
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;  
  typedef std::vector<TripletType> TripletListType;

public:
  DeformationBase( ){}
  
  virtual ~DeformationBase () {}
    
  virtual void applyEnergy ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, RealType& Dest ) const = 0;

  RealType operator() (const VectorType& UndeformedGeom, const VectorType& DeformedGeom) const {
    RealType Energy;
    applyEnergy(UndeformedGeom, DeformedGeom, Energy);
    return Energy;
  }
  
  void applyAddEnergy ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, RealType& Dest, RealType factor = 1.0 ) const {
    RealType Temp;
    applyEnergy( UndeformedGeom, DeformedGeom, Temp );
//    Dest.addMultiple( Temp, factor );
    Dest += Temp * factor;
  }
  
  virtual void applyUndefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest ) const = 0;
  
  void applyAddUndefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, Eigen::Ref<VectorType> Dest, RealType factor = 1.0 ) const {
    VectorType Temp( UndeformedGeom.size() );
    applyUndefGradient( UndeformedGeom, DeformedGeom, Temp );
    Dest += factor * Temp;
  }
  
  virtual void applyDefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest ) const = 0;
  
  void applyAddDefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, Eigen::Ref<VectorType> Dest, RealType factor = 1.0 ) const {
    VectorType Temp( UndeformedGeom.size() );
    applyDefGradient( UndeformedGeom, DeformedGeom, Temp );
    Dest += factor * Temp;
  }
   
  void applyDefHessian ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, MatrixType& Hessian, RealType factor = 1.0 ) const {
    int dofs = UndeformedGeom.size();
    if( (Hessian.rows() != dofs) || (Hessian.cols() != dofs) )
        Hessian.resize( dofs, dofs );
    Hessian.setZero();
    
    // set up triplet list
    TripletListType tripletList;
    tripletList.reserve( numOfNonZeroHessianEntries() );    
    pushTripletsDefHessian ( UndeformedGeom, DeformedGeom, tripletList, 0, 0, factor );
    
    // fill matrix from triplets
    Hessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  } 
  
  void applyUndefHessian ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, MatrixType& Hessian, RealType factor = 1.0 ) const {
    int dofs = UndeformedGeom.size();
    if( (Hessian.rows() != dofs) || (Hessian.cols() != dofs) )
        Hessian.resize( dofs, dofs );
    Hessian.setZero();
    
    // set up triplet list
    TripletListType tripletList;
    tripletList.reserve( numOfNonZeroHessianEntries() );    
    pushTripletsUndefHessian ( UndeformedGeom, DeformedGeom, tripletList, 0, 0, factor );
    
    // fill matrix from triplets
    Hessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  } 
  
  void applyMixedHessian ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, MatrixType& Hessian, const bool FirstDerivWRTDef, RealType factor = 1.0 ) const {
    int dofs = UndeformedGeom.size();
    if( (Hessian.rows() != dofs) || (Hessian.cols() != dofs) )
        Hessian.resize( dofs, dofs );
    Hessian.setZero();
    
    // set up triplet list
    TripletListType tripletList;
    tripletList.reserve( numOfNonZeroHessianEntries() );    
    pushTripletsMixedHessian ( UndeformedGeom, DeformedGeom, tripletList, 0, 0, FirstDerivWRTDef, factor );
    
    // fill matrix from triplets
    Hessian.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  } 
   
  virtual void pushTripletsDefHessian ( const VectorType& UndefGeom, const VectorType& DefGeom, TripletListType& triplets, int rowOffset, int colOffset, RealType factor = 1.0 ) const = 0; 
   
  virtual void pushTripletsUndefHessian ( const VectorType& UndefGeom, const VectorType& DefGeom, TripletListType& triplets, int rowOffset, int colOffset, RealType factor = 1.0 ) const = 0; 
  
  // mixed second derivative of deformation energy E[S_1, S_2], i.e. if "FirstDerivWRTDef" we have D_1 D_2 E[.,.], otherwise D_2 D_1 E[.,.]
  virtual void pushTripletsMixedHessian ( const VectorType& UndefGeom, const VectorType& DefGeom, TripletListType& triplets, int rowOffset, int colOffset, const bool FirstDerivWRTDef, RealType factor = 1.0 ) const = 0; 
  
  // should be overwritenn to imporve performance
  virtual int numOfNonZeroHessianEntries () const {
    return 0;    
  }
  
  virtual RealType getBendingWeight() const {
    throw BasicException("DeformationBase::getBendingWeight() has to be implemented in derived function!");    
  }
};


/**
 * \brief Weighted sum of membrane and bending deformation
 * \author Heeren
 */
template <typename ConfiguratorType, typename MembraneDeformationType, typename BendingDeformationType>
class ShellDeformation : public DeformationBase<ConfiguratorType>{
    
protected:    
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;
  
  const MeshTopologySaver& _topology;
  RealType _bendWeight, _memweight;

public:
  ShellDeformation( const MeshTopologySaver& Topology, RealType bendWeight ) : _topology( Topology ), _bendWeight( bendWeight ), _memweight(1.) {}
  
  ShellDeformation( const MeshTopologySaver& Topology, RealType memWeight, RealType bendWeight ) : _topology( Topology ), _bendWeight( bendWeight ), _memweight(memWeight) {}
  
  void applyEnergy ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, RealType& Dest ) const {
      MembraneDeformationType( _topology, _memweight ).applyEnergy( UndeformedGeom, DeformedGeom, Dest );
      BendingDeformationType( _topology, _bendWeight ).applyAddEnergy( UndeformedGeom, DeformedGeom, Dest );
  }
  
  void applyUndefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest ) const {
      MembraneDeformationType( _topology, _memweight ).applyUndefGradient( UndeformedGeom, DeformedGeom, Dest );
      BendingDeformationType( _topology, _bendWeight ).applyAddUndefGradient( UndeformedGeom, DeformedGeom, Dest );
  }
  
  void applyDefGradient ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, VectorType& Dest ) const {
      MembraneDeformationType( _topology, _memweight ).applyDefGradient( UndeformedGeom, DeformedGeom, Dest );
      BendingDeformationType( _topology, _bendWeight ).applyAddDefGradient( UndeformedGeom, DeformedGeom, Dest );
  }
  
  void pushTripletsDefHessian ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, TripletListType& triplets, int rowOffset, int colOffset, RealType factor = 1.0 ) const {
      //TODO parallelize?!
      MembraneDeformationType( _topology, _memweight  ).pushTripletsDefHessian( UndeformedGeom, DeformedGeom, triplets, rowOffset, colOffset, factor );
      BendingDeformationType(  _topology, _bendWeight ).pushTripletsDefHessian( UndeformedGeom, DeformedGeom, triplets, rowOffset, colOffset, factor );
  }
   
  void pushTripletsUndefHessian ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, TripletListType& triplets, int rowOffset, int colOffset, RealType factor = 1.0 ) const {
       //TODO parallelize?!
      MembraneDeformationType( _topology, _memweight  ).pushTripletsUndefHessian( UndeformedGeom, DeformedGeom, triplets, rowOffset, colOffset, factor );
      BendingDeformationType(  _topology, _bendWeight ).pushTripletsUndefHessian( UndeformedGeom, DeformedGeom, triplets, rowOffset, colOffset, factor );
  }
  
  // mixed second derivative of deformation energy E[S_1, S_2], i.e. if "FirstDerivWRTDef" we have D_1 D_2 E[.,.], otherwise D_2 D_1 E[.,.]
  void pushTripletsMixedHessian ( const VectorType& UndeformedGeom, const VectorType& DeformedGeom, TripletListType& triplets, int rowOffset, int colOffset, const bool FirstDerivWRTDef, RealType factor = 1.0 ) const {
      //TODO parallelize?!
      MembraneDeformationType( _topology, _memweight  ).pushTripletsMixedHessian( UndeformedGeom, DeformedGeom, triplets, rowOffset, colOffset, FirstDerivWRTDef, factor );
      BendingDeformationType(  _topology, _bendWeight ).pushTripletsMixedHessian( UndeformedGeom, DeformedGeom, triplets, rowOffset, colOffset, FirstDerivWRTDef, factor );
  }  
  
  int numOfNonZeroHessianEntries () const {
    return MembraneDeformationType( _topology, _memweight  ).numOfNonZeroHessianEntries () + BendingDeformationType(  _topology, _bendWeight ).numOfNonZeroHessianEntries();
  }
  
  RealType getBendingWeight() const {
      return _bendWeight;
  }
    
};



//==========================================================================================================
// DEBUGGING CLASSES FOR SHELL DEFORMATIONS
//==========================================================================================================

//=========================================================================================================
// UNDEFORMED

//! Energy to test deformation with respect to first argument
template <typename ConfiguratorType, typename DeformationType>
class UndeformedEnergyWrapper : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:       
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;

  const MeshTopologySaver& _topology;
  RealType _weight;
  const VectorType& _defArg;
  
public:  
  UndeformedEnergyWrapper( const MeshTopologySaver& Topology, RealType Weight, const VectorType& deformedArg ) : _topology(Topology), _weight(Weight), _defArg(deformedArg) {}
  
  void apply( const VectorType& UndeformedArg, RealType& Dest ) const {
    DeformationType( _topology, _weight ).applyEnergy ( UndeformedArg, _defArg, Dest );
  }
    
};

//! Gradient to test deformation with respect to first argument
template <typename ConfiguratorType, typename DeformationType>
class UndeformedGradientWrapper : public BaseOp<typename ConfiguratorType::VectorType> {
protected:       
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;

  const MeshTopologySaver& _topology;
  RealType _weight;
  const VectorType& _defArg;
  
public:  
  UndeformedGradientWrapper( const MeshTopologySaver& Topology, RealType Weight, const VectorType& deformedArg ) : _topology(Topology), _weight(Weight), _defArg(deformedArg) {}
  
  void apply( const VectorType& UndeformedArg, VectorType& Dest ) const {
    DeformationType( _topology, _weight ).applyUndefGradient( UndeformedArg, _defArg, Dest );
  }
    
};

//! Hessian to test deformation with respect to first argument
template <typename ConfiguratorType, typename DeformationType>
class UndeformedHessianWrapper : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:     
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;  
  typedef std::vector<TripletType> TripletListType;
  
  const MeshTopologySaver& _topology;
  RealType _weight;
  const VectorType& _defArg;
  
public:  
  UndeformedHessianWrapper( const MeshTopologySaver& Topology, RealType Weight, const VectorType& deformedArg ) : _topology(Topology), _weight(Weight), _defArg(deformedArg) {}
  
  void apply( const VectorType& UndeformedArg, MatrixType& Dest ) const {
      Dest.setZero();
      DeformationType deform( _topology, _weight );
      TripletListType tripletList;
      tripletList.reserve( deform.numOfNonZeroHessianEntries() );
      deform.pushTripletsUndefHessian( UndeformedArg, _defArg, tripletList, 0, 0 );
      Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
    
};

//=========================================================================================================
// DEFORMED

//! Energy to test deformation with respect to second argument
template <typename ConfiguratorType, typename DeformationType>
class DeformedEnergyWrapper : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
protected:       
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;

  const MeshTopologySaver& _topology;
  RealType _weight;
  const VectorType& _undefArg;
  
public:  
  DeformedEnergyWrapper( const MeshTopologySaver& Topology, RealType Weight, const VectorType& undeformedArg ) : _topology(Topology), _weight(Weight), _undefArg(undeformedArg) {}
  
  void apply( const VectorType& DeformedArg, RealType& Dest ) const {
    DeformationType( _topology, _weight ).applyEnergy ( _undefArg, DeformedArg, Dest );
  }
    
};

//! Gradient to test deformation with respect to second argument
template <typename ConfiguratorType, typename DeformationType>
class DeformedGradientWrapper : public BaseOp<typename ConfiguratorType::VectorType> {
protected:       
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;

  const MeshTopologySaver& _topology;
  RealType _weight;
  const VectorType& _undefArg;
  
public:  
  DeformedGradientWrapper( const MeshTopologySaver& Topology, RealType Weight, const VectorType& undeformedArg ) : _topology(Topology), _weight(Weight), _undefArg(undeformedArg) {}
  
  void apply( const VectorType& DeformedArg, VectorType& Dest ) const {
    DeformationType( _topology, _weight ).applyDefGradient( _undefArg, DeformedArg, Dest );
  }
    
};

//! Hessian to test deformation with respect to second argument
template <typename ConfiguratorType, typename DeformationType>
class DeformedHessianWrapper : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:     
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;  
  typedef std::vector<TripletType> TripletListType;
  
  const MeshTopologySaver& _topology;
  RealType _weight;
  const VectorType& _undefArg;
  
public:  
  DeformedHessianWrapper( const MeshTopologySaver& Topology, RealType Weight, const VectorType& undeformedArg ) : _topology(Topology), _weight(Weight), _undefArg(undeformedArg) {}
  
  void apply( const VectorType& DeformedArg, MatrixType& Dest ) const {
      Dest.setZero();
      DeformationType deform( _topology, _weight );
      TripletListType tripletList;
      tripletList.reserve( deform.numOfNonZeroHessianEntries() );
      deform.pushTripletsDefHessian( _undefArg, DeformedArg, tripletList, 0, 0 );
      Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
    
};

//=========================================================================================================
// MIXED I: test derivative of x \mapsto d_2W[x,y]

//! Gradient to test deformation with respect to second argument (with undeformed as argument)
template <typename ConfiguratorType, typename DeformationType>
class DeformedGradientUndefArgWrapper : public BaseOp<typename ConfiguratorType::VectorType> {
protected:       
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;

  const MeshTopologySaver& _topology;
  RealType _weight;
  const VectorType& _defArg;
  
public:  
  DeformedGradientUndefArgWrapper( const MeshTopologySaver& Topology, RealType Weight, const VectorType& deformedArg ) : _topology(Topology), _weight(Weight), _defArg(deformedArg) {}
  
  void apply( const VectorType& UndeformedArg, VectorType& Dest ) const {
    DeformationType( _topology, _weight ).applyDefGradient( UndeformedArg, _defArg, Dest );
  }
    
};

//! Mixed Hessian of test deformation (first derivative wrt deformed argument, i.e. here the argument is undeformed)
template <typename ConfiguratorType, typename DeformationType>
class MixedHessianUndefArgWrapper : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:       
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;  
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver& _topology;
  RealType _weight;
  const VectorType& _defArg;
  
public:  
  MixedHessianUndefArgWrapper( const MeshTopologySaver& Topology, RealType Weight, const VectorType& deformedArg ) : _topology(Topology), _weight(Weight), _defArg(deformedArg) {}
  
  void apply( const VectorType& UndeformedArg, MatrixType& Dest ) const {
      int numDOFs = UndeformedArg.size();
      if( Dest.rows() != numDOFs || Dest.cols() != numDOFs )
          Dest.resize(numDOFs, numDOFs );
      Dest.setZero();
      
      DeformationType deform( _topology, _weight );
      TripletListType tripletList;
      tripletList.reserve( deform.numOfNonZeroHessianEntries() );
      // pushTripletsMixedHessian ( const VectorType& UndefGeom, const VectorType& DefGeom, TripletListType& triplets, int rowOffset, int colOffset, const bool FirstDerivWRTDef, RealType factor = 1.0 )
      bool FirstDerivWRTDef = true;
      deform.pushTripletsMixedHessian( UndeformedArg, _defArg, tripletList, 0, 0, FirstDerivWRTDef );
      Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
    
};

//=========================================================================================================
// MIXED II: test derivative of y \mapsto d_1W[x,y]

//! Gradient to test deformation with respect to first argument (with deformed as argument)
template <typename ConfiguratorType, typename DeformationType>
class UndeformedGradientDefArgWrapper : public BaseOp<typename ConfiguratorType::VectorType> {
protected:       
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;

  const MeshTopologySaver& _topology;
  RealType _weight;
  const VectorType& _undefArg;
  
public:  
  UndeformedGradientDefArgWrapper( const MeshTopologySaver& Topology, RealType Weight, const VectorType& undeformedArg ) : _topology(Topology), _weight(Weight), _undefArg(undeformedArg) {}
  
  void apply( const VectorType& DeformedArg, VectorType& Dest ) const {
    DeformationType( _topology, _weight ).applyUndefGradient( _undefArg, DeformedArg, Dest );
  }
    
};

//! Mixed Hessian of test deformation (first derivative wrt undeformed argument, i.e. here the argument is deformed)
template <typename ConfiguratorType, typename DeformationType>
class MixedHessianDefArgWrapper : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
protected:       
  typedef typename ConfiguratorType::RealType    RealType;
  typedef typename ConfiguratorType::VectorType  VectorType;
  typedef typename ConfiguratorType::SparseMatrixType  MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;  
  typedef std::vector<TripletType> TripletListType;

  const MeshTopologySaver& _topology;
  RealType _weight;
  const VectorType& _undefArg;
  
public:  
  MixedHessianDefArgWrapper( const MeshTopologySaver& Topology, RealType Weight, const VectorType& undeformedArg ) : _topology(Topology), _weight(Weight), _undefArg(undeformedArg) {}
  
  void apply( const VectorType& DeformedArg, MatrixType& Dest ) const {
      int numDOFs = DeformedArg.size();
      if( Dest.rows() != numDOFs || Dest.cols() != numDOFs )
          Dest.resize(numDOFs, numDOFs );
      Dest.setZero();
      
      DeformationType deform( _topology, _weight );
      TripletListType tripletList;
      tripletList.reserve( deform.numOfNonZeroHessianEntries() );
      // pushTripletsMixedHessian ( const VectorType& UndefGeom, const VectorType& DefGeom, TripletListType& triplets, int rowOffset, int colOffset, const bool FirstDerivWRTDef, RealType factor = 1.0 )
      bool FirstDerivWRTDef = false;
      deform.pushTripletsMixedHessian( _undefArg, DeformedArg, tripletList, 0, 0, FirstDerivWRTDef );
      Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
    
};

#endif

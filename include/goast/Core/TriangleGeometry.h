// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
//=============================================================================
//
//  Geometric functions and derivatives related to a flap of triangles
//
//=============================================================================

#ifndef TRIANGLEGEOMETRY_HH
#define TRIANGLEGEOMETRY_HH

#include "Auxiliary.h"
#include "SmallVecMat.h"

//=============================================================================
// Geometric operators in triangle flap
//=============================================================================

template<typename VecType>
void getAreaGradient( const VecType& Pi, const VecType& Pj, const VecType& Pk, VecType& grad  ) {
  VecType normal;
  normal.makeCrossProduct(Pk-Pj, Pi-Pk);
  normal /= 2. * normal.norm();
  grad.makeCrossProduct( normal, Pj - Pi);  
}
/*
// D_kN = (2A)^{-1} (Id - NN^T) (Ex)^T, with Ex \in \R^{3,3} s.t. Ex(V) is the cross prodcut between E = Pi - Pj and V, A is area of triangle
template<typename RealType>
void getNormalGradient( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, SmallMat33<RealType>& Dn  ) {
  SmallVec3<RealType> normal;
  normal.makeCrossProduct(Pk-Pj, Pi-Pk);
  RealType twoTimesArea = normal.norm();
  normal /= twoTimesArea;
  SmallVec3<RealType> edge(Pi - Pj);
  edge /= -1. * twoTimesArea;
  SmallMat33<RealType> crossOp;
  getCrossOp( edge, crossOp );
  getProjection( normal, Dn );
  Dn *= crossOp;  
}
*/

template<typename RealType>
void getNormalGradientPk( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, SmallMat33<RealType>& Dn ) {
    SmallVec3<RealType> normal;
    RealType twoTimesArea = getWeightedNormal( Pi, Pj, Pk, normal, true);
    getNormalGradientPk(Pi, Pj, normal, twoTimesArea, Dn );
}

template<typename RealType>
void getNormalGradientPk( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Nk, RealType TwoTimesArea, SmallMat33<RealType>& Dn ) {
    SmallVec3<RealType> Ek = Pi - Pj;
    SmallVec3<RealType> nCrossEk;
    nCrossEk.makeCrossProduct(Nk, Ek);
    Dn.makeTensorProduct( nCrossEk, Nk);
    Dn /= TwoTimesArea;
}

template<typename RealType, typename TensorType>
void getNormalHessian( SmallVec3<RealType> &Pi, const SmallVec3<RealType> &Pj, const SmallVec3<RealType> &Pk,
                       TensorType &D2n ) {
  D2n.resize( 3, 9, 9 );

  std::array<SmallVec3<RealType>, 3> P{ Pi, Pj, Pk };

  // Normal
  SmallVec3<RealType> normal;
  normal.makeCrossProduct(P[2] - P[1], P[0] - P[2]);

  // Area
  RealType twoTimesArea = normal.norm();
  normal /= twoTimesArea;


  // First derivatives
  std::array<SmallMat33<RealType>, 3> Dn;
  std::array<SmallVec3<RealType>, 3> DA;


  for ( const int k : { 0, 1, 2 } ) {
    SmallVec3<RealType>  edge = P[(k + 1) % 3] - P[(k + 2) % 3];
    SmallVec3<RealType> eCrossN;
    eCrossN.makeCrossProduct(edge, normal );

    Dn[k].makeTensorProduct(eCrossN, normal);
    Dn[k] *= -1. / twoTimesArea;

    Dn[k].transpose(); // we will need columns later but SmallMat only supports row access

    DA[k].makeCrossProduct(normal, edge);
    DA[k] *= -0.5;
  }


  // Second derivatives
  for ( const int k1 : { 0, 1, 2 } ) { // Inner index
    SmallVec3<RealType> edge = P[(k1 + 1) % 3] - P[(k1 + 2) % 3];

    for ( const int k2 : { 0, 1, 2 } ) { // Outer index
      SmallMat33<RealType> De;
      De.setZero();

      if ( k1 != k2 )
        De.setIdentity();
      if ( k2 == ((k1 + 2) % 3))
        De *= -1.;

      for ( const int i : { 0, 1, 2 } ) { // Coordinate of outer index
        SmallMat33<RealType> local_D2n;
        local_D2n.setZero();

        SmallVec3<RealType> auxCrossProd;
        SmallMat33<RealType> auxTensorProd;

        if ( k1 != k2 ) {
          auxCrossProd.makeCrossProduct(De[i], normal);
          auxTensorProd.makeTensorProduct(auxCrossProd, normal);
          auxTensorProd *= twoTimesArea;
          local_D2n += auxTensorProd;
        }
        auxCrossProd.makeCrossProduct(edge, Dn[k2][i]);
        auxTensorProd.makeTensorProduct(auxCrossProd, normal);
        auxTensorProd *= twoTimesArea;
        local_D2n += auxTensorProd;

        auxCrossProd.makeCrossProduct(edge, normal);
        auxTensorProd.makeTensorProduct(auxCrossProd, Dn[k2][i]);
        auxTensorProd *= twoTimesArea;
        local_D2n += auxTensorProd;

        auxCrossProd.makeCrossProduct(edge, normal);
        auxTensorProd.makeTensorProduct(auxCrossProd, normal);
        auxTensorProd *= DA[k2][i] * 2;
        local_D2n -= auxTensorProd;

        local_D2n *= -1 / (twoTimesArea * twoTimesArea);

        for ( const int r : { 0, 1, 2 } ) {
          for ( const int c : { 0, 1, 2 } ) {
            D2n( r, c * 3 + k1, i * 3 + k2 ) = local_D2n( r, c );
          }
        }
      }
    }

  }
}

//
template<typename RealType>
void getCosThetaGradK( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, const SmallVec3<RealType>& n, SmallVec3<RealType>& grad  ) {
  SmallMat33<RealType> Dn;
  getNormalGradient( Pi, Pj, Pk, Dn );
  Dn.multTransp( n, grad );
}

//
template<typename RealType>
RealType getCosOfDihedralAngle( const SmallVec3<RealType>& Pi,
                                const SmallVec3<RealType>& Pj,
                                const SmallVec3<RealType>& Pk,
                                const SmallVec3<RealType>& Pl ) {
    SmallVec3<RealType> nk, nl;
    nk.makeCrossProduct( Pk-Pj, Pi-Pk );
    nk.normalize();
    nl.makeCrossProduct( Pl-Pi, Pj-Pl );
    nl.normalize();
    return (nk*nl);
}
 
// returns dihedral angle between n1 and n2 at an edge e
template<typename RealType>
RealType getDihedralAngle( const SmallVec3<RealType>& nk,
                           const SmallVec3<RealType>& nl,
                           const SmallVec3<RealType>& e ) {
    SmallVec3<RealType> crossprod;
    crossprod.makeCrossProduct( nk, nl);
    //return std::asin( temp*e / e.norm() );
    RealType aux = std::max( std::min( nk*nl, (RealType)1. ), (RealType)-1.);
    return crossprod*e < 0. ? -std::acos( aux ) : std::acos( aux );
}

//
template<typename RealType>
RealType getDihedralAngle( const SmallVec3<RealType>& Pi,
                           const SmallVec3<RealType>& Pj,
                           const SmallVec3<RealType>& Pk,
                           const SmallVec3<RealType>& Pl ) {
    SmallVec3<RealType> nk, nl;
    nk.makeCrossProduct( Pk-Pj, Pi-Pk );
    nk.normalize();
    nl.makeCrossProduct( Pl-Pi, Pj-Pl );
    nl.normalize();
    return getDihedralAngle( nk, nl, Pj-Pi );
}

// interior angle in triangle
template<typename RealType>
RealType getInnerAngle( const SmallVec3<RealType>& p, const SmallVec3<RealType>& q ) { 
    RealType angle = std::acos( p*q / ( p.norm() * q.norm() ) );
    if( angle < 0 )
      throw BasicException ( "getInnerAngle(): interior angle in triangle is smaller than pi!!" );
    return( angle );     
}
  
// 
template<typename RealType>
void getCrossOp( const SmallVec3<RealType>& a, SmallMat33<RealType>& matrix ) {
  matrix.setZero();
  matrix.set( 0, 1, -a[2]); matrix.set( 0, 2,  a[1]);
  matrix.set( 1, 0,  a[2]); matrix.set( 1, 2, -a[0]);
  matrix.set( 2, 0, -a[1]); matrix.set( 2, 1,  a[0]);
}

template<typename RealType>
void getProjection( const SmallVec3<RealType>& x, SmallMat33<RealType>& m ) {
  m.setIdentity();
  SmallMat33<RealType> temp;
  temp.makeTensorProduct( x, x );
  m.addMultiple( temp, -1.0);
}

template<typename RealType>  
void getReflection( const SmallVec3<RealType>& x, SmallMat33<RealType>& m ) {
  m.setIdentity();
  SmallMat33<RealType> temp;
  temp.makeTensorProduct( x, x );
  m.addMultiple( temp, -2.0);
}

template<typename RealType>
RealType getAreaSqr( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk ) {
  SmallVec3<RealType> temp;   
  temp.makeCrossProduct( Pk-Pj, Pi-Pk );
  return 0.25 * temp.normSqr();
}

template<typename RealType>
RealType getArea( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk ) {
  return std::sqrt( getAreaSqr<RealType>( Pi, Pj, Pk ) );
}

template<typename RealType>
RealType getWeightedNormal( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, SmallVec3<RealType>& normal, bool normalize = false  ) {
  normal.makeCrossProduct( Pk-Pj, Pi-Pk);
  RealType val = normal.norm();
  if( normalize )
      normal /= val;
  return val;
}

template<typename RealType>
void getNormal( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, SmallVec3<RealType>& normal  ) {
  normal.makeCrossProduct( Pk-Pj, Pi-Pk);
  normal.normalize();
}

template<typename RealType>
void getAreaGradK( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, SmallVec3<RealType>& grad  ) {
  SmallVec3<RealType> a(Pi-Pk), d(Pk-Pj), e(Pj-Pi);
  RealType area = getArea( Pi, Pj, Pk );
  RealType temp1( -0.25 * dotProduct(e,a) / area ), temp2( 0.25 * dotProduct(e,d) / area );
  getWeightedVectorSum( temp1, d, temp2, a, grad );
}

template<typename RealType>
void getLengthGradk( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, SmallVec3<RealType>& grad  ) {
  grad = Pi - Pj;
  grad.normalize();
}

template<typename RealType>
void getThetaGradK( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, SmallVec3<RealType>& grad  ) {
  SmallVec3<RealType> e(Pj-Pi);
  getNormal( Pi, Pj, Pk, grad );
  grad *= -0.5 * e.norm() / getArea( Pi, Pj, Pk );
}

template<typename RealType>
void getThetaGradILeftPart( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, SmallVec3<RealType>& grad  ) {
  SmallVec3<RealType> e(Pj-Pi), d(Pk-Pj);
  getThetaGradK( Pi, Pj, Pk, grad );
  grad *= dotProduct(d,e) / dotProduct(e,e);
}

template<typename RealType>
void getThetaGradJLeftPart( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, SmallVec3<RealType>& grad  ) {
  SmallVec3<RealType> e(Pj-Pi), a(Pi-Pk);
  getThetaGradK( Pi, Pj, Pk, grad );
  grad *= dotProduct(a,e) / dotProduct(e,e);
}

template<typename RealType>
void getThetaGradI( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, const SmallVec3<RealType>& Pl, SmallVec3<RealType>& grad  ) {
  SmallVec3<RealType> temp;
  getThetaGradILeftPart( Pi, Pj, Pk, grad );
  getThetaGradILeftPart( Pi, Pj, Pl, temp );
  grad -= temp;
}

template<typename RealType>
void getThetaGradJ( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, const SmallVec3<RealType>& Pl, SmallVec3<RealType>& grad  ) {
  SmallVec3<RealType> temp;
  getThetaGradJLeftPart( Pi, Pj, Pk, grad );
  getThetaGradJLeftPart( Pi, Pj, Pl, temp );
  grad -= temp;
}


template<typename RealType>
void getHessAreaKK( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, SmallMat33<RealType>& Hess ) {
  SmallVec3<RealType> eNormalized(Pj-Pi), gradAreaK;
  SmallMat33<RealType> proj;
  eNormalized.normalize();
  getAreaGradK( Pi, Pj, Pk, gradAreaK );
  Hess.makeTensorProduct( gradAreaK, gradAreaK );
  getProjection( eNormalized, proj );
  Hess.addMultiple( proj, -0.25 * dotProduct(Pj-Pi,Pj-Pi) );
  Hess *= -1. / getArea( Pi, Pj, Pk );
}

template<typename RealType>
void getHessAreaIK( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, SmallMat33<RealType>& Hess ) {
  SmallVec3<RealType> e(Pj-Pi), d(Pk-Pj), temp1, temp2;
  getAreaGradK( Pj, Pk, Pi, temp1 );
  getAreaGradK( Pi, Pj, Pk, temp2 );
  Hess.makeTensorProduct( temp1, temp2 );
  SmallMat33<RealType> auxMat;  
  auxMat.makeTensorProduct( e, d );
  Hess.addMultiple( auxMat, 0.25 );  
  Hess.addToDiagonal( -0.25 * dotProduct(d,e) );     
  Hess *= -1. / getArea( Pi, Pj, Pk );  
  getNormal( Pi, Pj, Pk, temp1);
  getCrossOp( temp1, auxMat );
  Hess.addMultiple( auxMat, 0.5 );
}


template<typename RealType>
void getHessThetaKK( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, SmallMat33<RealType>& Hkk ) {
  
  RealType areaSqr = getArea( Pi, Pj, Pk ) * getArea( Pi, Pj, Pk );
   
  SmallVec3<RealType> e(Pj-Pi), gradArea, normal;
  getAreaGradK( Pi, Pj, Pk, gradArea );
  getNormal( Pi, Pj, Pk, normal );
        
  SmallMat33<RealType> mat1, mat2;    
  getCrossOp( e, mat1 ); 
  mat2.makeTensorProduct( gradArea, normal );
    
  getWeightedMatrixSum<RealType>( e.norm() / (4. * areaSqr), mat1,  e.norm() / areaSqr, mat2, Hkk );
}

template<typename RealType>
void getHessThetaIK( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, SmallMat33<RealType>& Hik ) {
        
  RealType area = getArea( Pi, Pj, Pk );
  RealType areaSqr = area * area;
   
  SmallVec3<RealType> e(Pj-Pi), d(Pk-Pj), gradArea, normal;    
  getAreaGradK( Pj, Pk, Pi, gradArea );
  getNormal( Pi, Pj, Pk, normal );
    
  SmallMat33<RealType> mat1, mat2, mat3;    
  mat1.makeTensorProduct( e, normal );
  getCrossOp( d, mat2 );    
  getWeightedMatrixSum<RealType>( 1. / (2.*area*e.norm()), mat1,  e.norm() / (4.*areaSqr), mat2, mat3 );
    
  mat1.makeTensorProduct( gradArea, normal );
  getWeightedMatrixSum<RealType>( 1., mat3, e.norm() / areaSqr, mat1, Hik );
}

template<typename RealType>
void getHessThetaJK( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, SmallMat33<RealType>& Hjk ) {
    
  RealType area = getArea( Pi, Pj, Pk );
  RealType areaSqr = area * area;
    
  SmallVec3<RealType> e(Pi-Pj), a(Pi-Pk), gradArea, normal;      
  getAreaGradK( Pk, Pi, Pj, gradArea );
  getNormal( Pi, Pj, Pk, normal );
    
  SmallMat33<RealType> mat1, mat2, mat3;    
  mat1.makeTensorProduct( e, normal );
  getCrossOp( a, mat2 );    
  getWeightedMatrixSum<RealType>( 1. / (2.*area*e.norm()), mat1,  e.norm() / (4.*areaSqr), mat2, mat3 );
    
  mat1.makeTensorProduct( gradArea, normal );
  getWeightedMatrixSum<RealType>( 1., mat3, e.norm() / areaSqr, mat1, Hjk );
}

template<typename RealType>
void getHessThetaILeftPartI( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, SmallMat33<RealType>& HiLeft ) {
  SmallVec3<RealType> e(Pj-Pi), d(Pk-Pj), eNormalized(Pj-Pi), gradThetaK, temp;
  eNormalized.normalize();
  SmallMat33<RealType> mat1, mat2, Refl;
  getThetaGradK( Pi, Pj, Pk, gradThetaK );
  getReflection( eNormalized, Refl );
  Refl.mult( d, temp );
  mat1.makeTensorProduct( temp, gradThetaK );
  getHessThetaIK( Pi, Pj, Pk, mat2 );
  getWeightedMatrixSum<RealType>( -1. / dotProduct(e,e), mat1, dotProduct(d,e) / dotProduct(e,e), mat2, HiLeft );
}

template<typename RealType>
void getHessThetaJLeftPartI( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, SmallMat33<RealType>& HjLeft ) {
  SmallVec3<RealType> e(Pj-Pi), d(Pk-Pj), eNormalized(Pj-Pi), gradThetaK, temp, thetak;
  eNormalized.normalize();
  SmallMat33<RealType> mat1, mat2, Refl;
  getThetaGradK( Pi, Pj, Pk, gradThetaK );
  getReflection( eNormalized, Refl );
  Refl.mult( d-e, temp );
  mat1.makeTensorProduct( temp, gradThetaK );
  getHessThetaJK( Pi, Pj, Pk, mat2 );
  getWeightedMatrixSum<RealType>( 1. / dotProduct(e,e), mat1, dotProduct(d,e) / dotProduct(e,e), mat2, HjLeft );
 }

template<typename RealType>
void getHessThetaII( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, const SmallVec3<RealType>& Pl, SmallMat33<RealType>& Hii ) {    
  SmallMat33<RealType> temp;
  getHessThetaILeftPartI(Pi, Pj, Pk, Hii);
  getHessThetaILeftPartI(Pi, Pj, Pl, temp);
  Hii -= temp;
}

template<typename RealType>
void getHessThetaJI( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, const SmallVec3<RealType>& Pl, SmallMat33<RealType>& Hji ) {
  SmallVec3<RealType> edge(Pj-Pi), d(Pk-Pj), c(Pj-Pl);
  SmallVec3<RealType> diff(d-edge), sum(c+edge);
  RealType eLengthSqr = edge.normSqr();
  
  //
  SmallVec3<RealType> thetak, thetal, grad;
  getThetaGradK( Pi, Pj, Pk, thetak );
  getThetaGradK( Pj, Pi, Pl, thetal );
  getWeightedVectorSum<RealType>( dotProduct(edge,d), thetak, -1. * dotProduct(edge,c), thetal, grad );
  
  SmallMat33<RealType> Hjk, Hjl, tensorProduct;
  getHessThetaJK( Pi, Pj, Pk, Hjk );
  getHessThetaIK( Pj, Pi, Pl, Hjl);
  
  // Hess part
  getWeightedMatrixSumTransposed<RealType>( dotProduct(edge,d) / eLengthSqr, Hjk, -1. * dotProduct(edge,c) / eLengthSqr, Hjl, Hji );
  
  tensorProduct.makeTensorProduct(grad, edge);
  Hji.addMultiple( tensorProduct, -2. / (eLengthSqr*eLengthSqr) );
  
  tensorProduct.makeTensorProduct( thetak, diff );
  Hji.addMultiple( tensorProduct, 1. / eLengthSqr  );
  tensorProduct.makeTensorProduct( thetal, sum );
  Hji.addMultiple( tensorProduct, -1. / eLengthSqr  );
}

//=============================================================================
// DIFFERENCE QUOTIENTS (ONLY FOR DEBUGGING ISSUES!)
//=============================================================================

template<typename RealType>
void getAreaIKDiffQuotient( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, SmallMat33<RealType>& Hik ) {
  RealType H = 1e-8;
  for( int i = 0; i < 3; i++ ){
    SmallVec3<RealType> PiPlusH( Pi ), grad1, grad2;    
    PiPlusH[i] += H;
    getAreaGradK( PiPlusH, Pj, Pk, grad1 );
    getAreaGradK( Pi, Pj, Pk, grad2 );
    grad1 -= grad2;
    for( int j = 0; j < 3; j++ )
      Hik.set( j, i, grad1[j] / H );
  }
}

template<typename RealType>
void getHiiDiffQuotient( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, const SmallVec3<RealType>& Pl, SmallMat33<RealType>& Hii ) {
  RealType H = 1e-8;
  for( int i = 0; i < 3; i++ ){
    SmallVec3<RealType> PiPlusH( Pi ), grad1, grad2;    
    PiPlusH[i] += H;
    getThetaGradI( PiPlusH, Pj, Pk, Pl, grad1 );
    getThetaGradI( Pi, Pj, Pk, Pl, grad2 );
    grad1 -= grad2;
    for( int j = 0; j < 3; j++ )
      Hii.set( j, i, grad1[j] / H );
  }
}

template<typename RealType>
void getHjjDiffQuotient( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, const SmallVec3<RealType>& Pl, SmallMat33<RealType>& Hjj ) {
  RealType H = 1e-8;
  for( int i = 0; i < 3; i++ ){
    SmallVec3<RealType> PjPlusH( Pj ), grad1, grad2;    
    PjPlusH[i] += H;
    getThetaGradJ( Pi, PjPlusH, Pk, Pl, grad1 );
    getThetaGradJ( Pi, Pj, Pk, Pl, grad2 );
    grad1 -= grad2;
    for( int j = 0; j < 3; j++ )
      Hjj.set( j, i, grad1[j] / H );
  }
}

template<typename RealType>
void getHijDiffQuotient( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, const SmallVec3<RealType>& Pl, SmallMat33<RealType>& Hij ) {
  RealType H = 1e-10;
  for( int i = 0; i < 3; i++ ){
    SmallVec3<RealType> PjPlusH( Pj ), PjMinusH( Pj ), grad1, grad2;    
    PjPlusH[i] += H;
    PjMinusH[i] -= H;
    getThetaGradI( Pi, PjPlusH, Pk, Pl, grad1 );
    getThetaGradI( Pi, PjMinusH, Pk, Pl, grad2 );
    grad1 -= grad2;
    for( int j = 0; j < 3; j++ )
      Hij.set( j, i, grad1[j] / (2*H) );
  }
}

template<typename RealType>
void getHjiDiffQuotient( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, const SmallVec3<RealType>& Pl, SmallMat33<RealType>& Hji ) {
  RealType H = 1e-8;
  for( int i = 0; i < 3; i++ ){
    SmallVec3<RealType> PiPlusH( Pi ), grad1, grad2;    
    PiPlusH[i] += H;
    getThetaGradJ( PiPlusH, Pj, Pk, Pl, grad1 );
    getThetaGradJ( Pi, Pj, Pk, Pl, grad2 );
    grad1 -= grad2;
    for( int j = 0; j < 3; j++ )
      Hji.set( j, i, grad1[j] / H );
  }
}

template<typename RealType>
void getHikDiffQuotient( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, const SmallVec3<RealType>& Pl, SmallMat33<RealType>& Hik ) {
  RealType H = 1e-8;
  for( int i = 0; i < 3; i++ ){
    SmallVec3<RealType> PiPlusH( Pi ), grad1, grad2;    
    PiPlusH[i] += H;
    getThetaGradK( PiPlusH, Pj, Pk, grad1 );
    getThetaGradK( Pi, Pj, Pk, grad2 );
    grad1 -= grad2;
    for( int j = 0; j < 3; j++ )
      Hik.set( j, i, grad1[j] / H );
  }
}

template<typename RealType>
void getHkiDiffQuotient( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, const SmallVec3<RealType>& Pl, SmallMat33<RealType>& Hki ) {
  RealType H = 1e-8;
  for( int i = 0; i < 3; i++ ){
    SmallVec3<RealType> PkPlusH( Pk ), grad1, grad2;    
    PkPlusH[i] += H;
    getThetaGradI( Pi, Pj, PkPlusH, Pl, grad1 );
    getThetaGradI( Pi, Pj, Pk, Pl, grad2 );
    grad1 -= grad2;
    for( int j = 0; j < 3; j++ )
      Hki.set( j, i, grad1[j] / H );
  }
}


#endif

// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
//=============================================================================
//
//  * geometric functions and derivatives related to a flap of triangles
//
//=============================================================================


#ifndef GEOMETRYEIGEN_HH
#define GEOMETRYEIGEN_HH


//== INCLUDES =================================================================
#include <iostream>
#include <EigenIncludes.h>
#include <Auxiliary.h>
//== FORWARDDECLARATIONS ======================================================

template <int dimension, typename FloatType>
class SmallEigenVec : public Eigen::Matrix<FloatType, dimension, 1> {

public:

    template<typename T>
    SmallEigenVec(const T &a) : Eigen::Matrix<FloatType, dimension, 1>(a)  {
    }

    SmallEigenVec() : Eigen::Matrix<FloatType, dimension, 1>()  {
    }

    //! add multiple of other vector
    template <typename T>
    void addMultiple ( const T& AddedVec, FloatType Factor ) {
        (*this) += Factor * AddedVec;
    }



    FloatType normSqr( ) const {
        return Eigen::Matrix<FloatType, dimension, 1>::squaredNorm();
    }

};

template <typename FloatType>
class SmallEigenVec3 : public SmallEigenVec<3, FloatType> {

public:
    template<typename T>
    SmallEigenVec3(const T &a) : SmallEigenVec<3, FloatType>(a)  {
    }

    SmallEigenVec3() : SmallEigenVec<3, FloatType>()  {
    }

    SmallEigenVec3(FloatType x, FloatType y, FloatType z) : SmallEigenVec<3, FloatType>()  {
        (*this) << x, y, z;
    }

    //! set this Vec3 to the cross product (vector product) of two other vectors
//    void makeCrossProduct ( const SmallEigenVec3<FloatType> &a, const SmallEigenVec3<FloatType> &b ) {
//        (*this) = a.cross(b);
//    }

    template <typename T1, typename T2>
    void makeCrossProduct ( const T1 &a, const T2 &b ) {
//        (*this)[0] = a[1] * b[2] - a[2] * b[1];
//        (*this)[1] = a[2] * b[0] - a[0] * b[2];
//        (*this)[2] = a[0] * b[1] - a[1] * b[0];
        auto c = a.cross(b).eval();
        (*this)[0] = c[0];
        (*this)[1] = c[1];
        (*this)[2] = c[2];

//        (*this) = SmallEigenVec3(a.cross(b));
    }

    SmallEigenVec3<FloatType> crossProduct ( const SmallEigenVec3<FloatType> &b ) const {
        SmallEigenVec3<FloatType> ret(cross(b));
        return ret;
    }

    SmallEigenVec3<FloatType>& operator= ( const Eigen::Matrix<FloatType, 3, 1> &other ) const {
        (*this) = other;
    }

};

template <int numRows, int numCols, typename FloatType>
class SmallEigenMat : public Eigen::Matrix<FloatType, numRows, numCols> {


public:

    template<typename T>
    SmallEigenMat(const T &a) : Eigen::Matrix<FloatType, numRows, numCols>(a)  {
    }

    SmallEigenMat() : Eigen::Matrix<FloatType, numRows, numCols>()  {
    }


    FloatType get ( const int I, const int J ) const {
        (*this)(I,J);
    }

    void set ( const int I, const int J, FloatType val ) {
        (*this)(I,J) = val;
    }

        // this =  A*B
    template <int dimBoth>
    void makeProduct ( const SmallEigenMat<numRows, dimBoth, FloatType> &A,
                       const SmallEigenMat<dimBoth, numCols, FloatType> &B ) {
            (*this) = A * B;
    }

    // this =  A^T*B
    template <int dimBoth>
    void makeProductAtransposedB ( const SmallEigenMat<dimBoth, numRows, FloatType> &A,
                                   const SmallEigenMat<dimBoth, numCols, FloatType> &B ) {
        (*this) = A.transpose() * B;
    }

    // this =  A*B^T
    template <int dimBoth>
    void makeProductABtransposed ( const SmallEigenMat<numRows, dimBoth, FloatType> &A,
                                   const SmallEigenMat<numCols, dimBoth, FloatType> &B ) {
        (*this) = A * B.transpose();

    }

    //! \f$ A \mapsto A + \alpha B \f$
    SmallEigenMat<numRows, numCols, FloatType> & addMultiple ( const SmallEigenMat<numRows, numCols, FloatType>& mat, const FloatType& alpha ) {
        (*this) += alpha * mat;
        return *this;
    }

    //! computes the trace of this matrix.
    FloatType tr() const {
        return Eigen::Matrix<FloatType, numRows, numCols>::trace();
    }

    //! computes the squared Frobius norm of this matrix.
    FloatType normSqr() const {
        return Eigen::Matrix<FloatType, numRows, numCols>::squaredNorm();
    }

    // this = a*b^T
    void makeTensorProduct ( const SmallEigenVec<numRows, FloatType> &a, const SmallEigenVec<numCols, FloatType> &b ) {
        (*this) = SmallEigenMat<numRows, numCols, FloatType>(a * b.transpose());
    }



    void multAdd ( const SmallEigenVec<numCols, FloatType> &Arg, SmallEigenVec<numRows, FloatType> &Dest ) const {
        Dest += (*this) * Arg;
    }

    void mult ( const SmallEigenVec<numCols, FloatType> &Arg, SmallEigenVec<numRows, FloatType> &Dest ) const {
        Dest = (*this) * Arg;
    }

    //! \todo numCols / numRows ?
    void multAddTransp ( const SmallEigenVec<numCols, FloatType> &Arg, SmallEigenVec<numRows, FloatType> &Dest ) const {
        Dest += (*this).transpose() * Arg;
    }

    void multTransp ( const SmallEigenVec<numCols, FloatType> &Arg, SmallEigenVec<numRows, FloatType> &Dest ) const {
        Dest = (*this).transpose() * Arg;
    }



};

template <typename FloatType>
class SmallEigenMat33 : public SmallEigenMat<3, 3, FloatType> {

public:
    template<typename T>
    SmallEigenMat33(const T &a) : SmallEigenMat<3, 3, FloatType>(a)  {
    }

    SmallEigenMat33() : SmallEigenMat<3, 3, FloatType>() { }

        SmallEigenMat33(FloatType xx, FloatType xy, FloatType xz, FloatType yx, FloatType yy, FloatType yz, FloatType zx,
                    FloatType zy, FloatType zz) : SmallEigenMat<3, 3, FloatType>() {
        (*this) << xx, xy, xz, yx, yy, yz, zx, zy, zz;
    }


    //! \todo what happens with inverse for non-invertible matrices?

    //! computes the determinant of this matrix.
    FloatType det() const {
        return Eigen::Matrix<FloatType, 3, 3>::determinant();
    }


//    //! \f$ A \mapsto A^T \f$
//    void transpose() {
//        Eigen::Matrix<FloatType, 3, 3>::transposeInPlace();
//    }
//
//    Eigen::Transpose< SmallEigenMat<3, 3, FloatType> > transposed() {
//        return SmallEigenMat<3, 3, FloatType>::transpose();
//    }


    //!
    void addToDiagonal( const FloatType Value ) {
        (*this).diagonal().array() += Value;
    }

    //
    void setCol( int idx, const SmallEigenVec3<FloatType>& col ) {
        (*this).col(idx) = col;
    }

    //
    void setRow( int idx, const SmallEigenVec3<FloatType>& row ) {
        (*this).row(idx) = row;
    }

//    SmallEigenMat33<FloatType>& operator= ( const Eigen::Inverse<Eigen::Matrix<double, 3, 3> > &other ) const {
//        Eigen::Matrix<double, 3, 3> temp(other);
//        (*this) = temp;
//    }

};

//
///////////////////////////////////////////////////////////
//// Input/Output operators
///////////////////////////////////////////////////////////
//template <typename T>
//inline std::ostream &operator<< ( std::ostream &os, const SmallMat33<T> &m ) {
//  return m.print ( os );
//}
//
///////////////////////////////////////////////////////////
//// Square and Cubic
///////////////////////////////////////////////////////////
//template<typename RealType>
//RealType SquareFnc( const RealType& val ){
//    return val*val;
//}
//
//template<typename RealType>
//RealType CubicFnc( const RealType& val ){
//    return val*val*val;
//}
///////////////////////////////////////////////////////////
//// Matrix-Vector operations
///////////////////////////////////////////////////////////
//
//// dot product
template<typename RealType>
RealType dotProduct(const Eigen::Matrix<RealType, 3, 1> &a, const Eigen::Matrix<RealType, 3, 1> &b) {
    return a.dot(b);
}

template<typename BinaryOp, typename LhsType, typename RhsType>
typename LhsType::Scalar dotProduct(const Eigen::CwiseBinaryOp<BinaryOp, LhsType, RhsType> &a, const Eigen::CwiseBinaryOp<BinaryOp, LhsType, RhsType> &b) {
    return a.dot(b);
}

template<typename RealType>
Eigen::Matrix<RealType, 3, 1> crossProduct(const Eigen::Matrix<RealType, 3, 1> &a,
                                           const Eigen::Matrix<RealType, 3, 1> &b) {
    return a.cross(b);
}

template<typename RealType>
Eigen::Matrix<RealType, 3, 3> tensorProduct(const Eigen::Matrix<RealType, 3, 1> &a,
                                            const Eigen::Matrix<RealType, 3, 1> &b) {
    return a * b.transpose();
}

template<typename RealType>
void addToDiagonal(Eigen::Matrix<RealType, 3, 3> &a, const RealType &b) {
    a.diagonal().array() += b;
}

//// weighted vector sum
template<typename RealType>
void getWeightedVectorSum(RealType a, const Eigen::Matrix<RealType, 3, 1> &vec1,
                          RealType b, const Eigen::Matrix<RealType, 3, 1> &vec2,
                          Eigen::Matrix<RealType, 3, 1> &res) {
    res = a * vec1 + b * vec2;
}

//// weighted matrix sum
template<typename RealType>
void getWeightedMatrixSum( RealType a, const Eigen::Matrix<RealType, 3, 3>& mat1,
                           RealType b, const Eigen::Matrix<RealType, 3, 3>& mat2, Eigen::Matrix<RealType, 3, 3>& res ) {
    res = a * mat1 + b * mat2;
}

//// weighted matrix sum (transposed)
template<typename RealType>
void getWeightedMatrixSumTransposed(RealType a, const Eigen::Matrix<RealType, 3, 3> &mat1,
                                    RealType b, const Eigen::Matrix<RealType, 3, 3> &mat2,
                                    Eigen::Matrix<RealType, 3, 3> &res) {
    res = a * mat1.transpose() + b * mat2.transpose();
}

///////////////////////////////////////////////////////////
//// Geometric operators in triangle flap
///////////////////////////////////////////////////////////

template<>
void getAreaGradient<SmallEigenVec3<double>>(const SmallEigenVec3<double> &Pi,
                                             const SmallEigenVec3<double> &Pj,
                                             const SmallEigenVec3<double> &Pk,
                                             SmallEigenVec3<double> &grad) {
//    Eigen::Matrix<double, 3, 1> normal;
//    normal = crossProduct<double>(Pk - Pj, Pi - Pk);
//    normal /= 2. * normal.norm();
    grad = 0.5 * (Pk - Pj).cross(Pi - Pk).normalized().cross(Pj - Pi);
}
//
//template<typename VecType>
//void getAreaGradient( const VecType& Pi, const VecType& Pj, const VecType& Pk, VecType& grad  ) {
//  VecType normal;
//  normal.makeCrossProduct(Pk-Pj, Pi-Pk);
//  normal /= 2. * normal.norm();
//  grad.makeCrossProduct( normal, Pj - Pi);
//}
//
//// D_kN = (2A)^{-1} (Id - NN^T) (Ex)^T, with Ex \in \R^{3,3} s.t. Ex(V) is the cross prodcut between E = Pi - Pj and V, A is area of triangle
//template<typename RealType>
//void getNormalGradient( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, SmallMat33<RealType>& Dn  ) {
//  SmallVec3<RealType> normal;
//  normal.makeCrossProduct(Pk-Pj, Pi-Pk);
//  RealType twoTimesArea = normal.norm();
//  normal /= twoTimesArea;
//  SmallVec3<RealType> edge(Pi - Pj);
//  edge /= -1. * twoTimesArea;
//  SmallMat33<RealType> crossOp;
//  getCrossOp( edge, crossOp );
//  getProjection( normal, Dn );
//  Dn *= crossOp;
//}
//
////
//template<typename RealType>
//void getCosThetaGradK( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, const SmallVec3<RealType>& n, SmallVec3<RealType>& grad  ) {
//  SmallMat33<RealType> Dn;
//  getNormalGradient( Pi, Pj, Pk, Dn );
//  Dn.multTransp( n, grad );
//}
//
////
//template<typename RealType>
//RealType getCosOfDihedralAngle( const SmallVec3<RealType>& Pi,
//                                const SmallVec3<RealType>& Pj,
//                                const SmallVec3<RealType>& Pk,
//                                const SmallVec3<RealType>& Pl ) {
//    SmallVec3<RealType> nk, nl;
//    nk.makeCrossProduct( Pk-Pj, Pi-Pk );
//    nk.normalize();
//    nl.makeCrossProduct( Pl-Pi, Pj-Pl );
//    nl.normalize();
//    return (nk*nl);
//}
//
//// returns dihedral angle between n1 and n2 at an edge e
template<typename RealType>
RealType getDihedralAngle( const Eigen::Matrix<RealType, 3, 1>& nk,
                           const Eigen::Matrix<RealType, 3, 1>& nl,
                           const Eigen::Matrix<RealType, 3, 1>& e ) {
//    Eigen::Matrix<RealType, 3, 1> crossprod = nk.cross(nl);
    //return std::asin( temp*e / e.norm() );
    RealType aux = std::max( std::min( nk.dot(nl), 1. ), -1.);
    return nk.cross(nl).dot(e) < 0. ? -std::acos( aux ) : std::acos( aux );
}

template<typename T1, typename T2, typename T3>
typename T1::Scalar getDihedralAngle(const T1 &nk,
                                     const T2 &nl,
                                     const T3 &e) {
//    Eigen::Matrix<RealType, 3, 1> crossprod = nk.cross(nl);
    //return std::asin( temp*e / e.norm() );
    typename T1::Scalar aux = std::max(std::min(nk.dot(nl), 1.), -1.);
    return nk.cross(nl).dot(e) < 0. ? -std::acos(aux) : std::acos(aux);
}

////
template<typename RealType>
RealType getDihedralAngle( const Eigen::Matrix<RealType, 3, 1>& Pi,
                           const Eigen::Matrix<RealType, 3, 1>& Pj,
                           const Eigen::Matrix<RealType, 3, 1>& Pk,
                           const Eigen::Matrix<RealType, 3, 1>& Pl ) {
    auto nk = (Pk-Pj).cross(Pi-Pk).normalized();
    auto nl = (Pl-Pi).cross(Pj-Pl).normalized();
    auto e = Pj - Pi;
    RealType aux = std::max(std::min(nk.dot(nl), 1.), -1.);
    return nk.cross(nl).dot(e) < 0. ? -std::acos(aux) : std::acos(aux);
//    return getDihedralAngle((Pk-Pj).cross(Pi-Pk).normalized(), (Pl-Pi).cross(Pj-Pl).normalized(), Pj - Pi);
}

////
template<typename RealType>
void getCrossOp(const Eigen::Matrix<RealType, 3, 1> &a, Eigen::Matrix<RealType, 3, 3> &matrix) {
    matrix.setZero();
    matrix(0, 1) = -a[2];
    matrix(0, 2) = a[1];
    matrix(1, 0) = a[2];
    matrix(1, 2) = -a[0];
    matrix(2, 0) = -a[1];
    matrix(2, 1) = a[0];
}

template<typename RealType>
void getProjection(const Eigen::Matrix<RealType, 3, 1> &x, Eigen::Matrix<RealType, 3, 3> &m) {
    m.setIdentity();
    Eigen::Matrix<RealType, 3, 3> temp;
    temp = tensorProduct(x, x);
    m -= temp;
}

template<typename RealType>
void getReflection(const Eigen::Matrix<RealType, 3, 1> &x, Eigen::Matrix<RealType, 3, 3> &m) {
    m.setIdentity();
    Eigen::Matrix<RealType, 3, 3> temp;
    temp = tensorProduct(x, x);
    m -= 2 * temp;
}

template<typename RealType>
RealType getArea(const Eigen::Matrix<RealType, 3, 1> &Pi, const Eigen::Matrix<RealType, 3, 1> &Pj,
                 const Eigen::Matrix<RealType, 3, 1> &Pk) {
    return 0.5 * (Pk - Pj).cross(Pi - Pk).norm();
}

//template<typename RealType>
//RealType getWeightedNormal( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, SmallVec3<RealType>& normal  ) {
//  normal.makeCrossProduct( Pk-Pj, Pi-Pk);
//  return normal.norm();
//}
//
template<typename RealType>
void getNormal(const Eigen::Matrix<RealType, 3, 1> &Pi, const Eigen::Matrix<RealType, 3, 1> &Pj,
               const Eigen::Matrix<RealType, 3, 1> &Pk, Eigen::Matrix<RealType, 3, 1> &normal) {
    normal = (Pk - Pj).cross(Pi - Pk);
    normal.normalize();
}

template<typename RealType>
void getAreaGradK(const Eigen::Matrix<RealType, 3, 1> &Pi, const Eigen::Matrix<RealType, 3, 1> &Pj,
                  const Eigen::Matrix<RealType, 3, 1> &Pk, Eigen::Matrix<RealType, 3, 1> &grad) {
    Eigen::Matrix<RealType, 3, 1> a(Pi - Pk), d(Pk - Pj), e(Pj - Pi);
    RealType area = getArea(Pi, Pj, Pk);
    RealType temp1(-0.25 * dotProduct(e, a) / area), temp2(0.25 * dotProduct(e, d) / area);
    grad = temp1 * d + temp2 * a;
}

template<typename RealType>
void getLengthGradk(const Eigen::Matrix<RealType, 3, 1> &Pi, const Eigen::Matrix<RealType, 3, 1> &Pj,
                    Eigen::Matrix<RealType, 3, 1> &grad) {
    grad = Pi - Pj;
    grad.normalize();
}

template<typename RealType>
void getThetaGradK(const Eigen::Matrix<RealType, 3, 1> &Pi, const Eigen::Matrix<RealType, 3, 1> &Pj,
                   const Eigen::Matrix<RealType, 3, 1> &Pk, Eigen::Matrix<RealType, 3, 1> &grad) {

    grad = (Pk - Pj).cross(Pi - Pk);
//    grad.normalize();
//    getNormal(Pi, Pj, Pk, grad);
    grad *= -1. * (Pj - Pi).norm() / std::pow(grad.norm(), 2);
}

template<typename RealType>
void getThetaGradILeftPart(const Eigen::Matrix<RealType, 3, 1> &Pi, const Eigen::Matrix<RealType, 3, 1> &Pj,
                           const Eigen::Matrix<RealType, 3, 1> &Pk, Eigen::Matrix<RealType, 3, 1> &grad) {
    Eigen::Matrix<RealType, 3, 1> e(Pj - Pi), d(Pk - Pj);
    getThetaGradK(Pi, Pj, Pk, grad);
    grad *= dotProduct(d, e) / dotProduct(e, e);
}

template<typename RealType>
void getThetaGradJLeftPart(const Eigen::Matrix<RealType, 3, 1> &Pi, const Eigen::Matrix<RealType, 3, 1> &Pj,
                           const Eigen::Matrix<RealType, 3, 1> &Pk, Eigen::Matrix<RealType, 3, 1> &grad) {
    Eigen::Matrix<RealType, 3, 1> e(Pj - Pi), a(Pi - Pk);
    getThetaGradK(Pi, Pj, Pk, grad);
    grad *= dotProduct(a, e) / dotProduct(e, e);
}

template<typename RealType>
void getThetaGradI(const Eigen::Matrix<RealType, 3, 1> &Pi, const Eigen::Matrix<RealType, 3, 1> &Pj,
                   const Eigen::Matrix<RealType, 3, 1> &Pk, const Eigen::Matrix<RealType, 3, 1> &Pl, 
                   Eigen::Matrix<RealType, 3, 1> &grad) {
    Eigen::Matrix<RealType, 3, 1> temp;
    getThetaGradILeftPart(Pi, Pj, Pk, grad);
    getThetaGradILeftPart(Pi, Pj, Pl, temp);
    grad -= temp;
}

template<typename RealType>
void getThetaGradJ(const Eigen::Matrix<RealType, 3, 1> &Pi, const Eigen::Matrix<RealType, 3, 1> &Pj,
                   const Eigen::Matrix<RealType, 3, 1> &Pk, const Eigen::Matrix<RealType, 3, 1> &Pl,
                   Eigen::Matrix<RealType, 3, 1> &grad) {
    Eigen::Matrix<RealType, 3, 1> temp;
    getThetaGradJLeftPart(Pi, Pj, Pk, grad);
    getThetaGradJLeftPart(Pi, Pj, Pl, temp);
    grad -= temp;
}


template<typename RealType>
void getHessAreaKK(const Eigen::Matrix<RealType, 3, 1> &Pi, const Eigen::Matrix<RealType, 3, 1> &Pj,
                   const Eigen::Matrix<RealType, 3, 1> &Pk, Eigen::Matrix<RealType, 3, 3> &Hess) {
    Eigen::Matrix<RealType, 3, 1> eNormalized(Pj - Pi), gradAreaK;
    Eigen::Matrix<RealType, 3, 3> proj;
    eNormalized.normalize();
    getAreaGradK(Pi, Pj, Pk, gradAreaK);
    Hess = tensorProduct(gradAreaK, gradAreaK);
    getProjection(eNormalized, proj);
    Hess += proj * -0.25 * dotProduct<RealType>(Pj - Pi, Pj - Pi);
    Hess *= -1. / getArea(Pi, Pj, Pk);
}

template<typename RealType>
void getHessAreaIK(const Eigen::Matrix<RealType, 3, 1> &Pi, const Eigen::Matrix<RealType, 3, 1> &Pj,
                   const Eigen::Matrix<RealType, 3, 1> &Pk, Eigen::Matrix<RealType, 3, 3> &Hess) {
    Eigen::Matrix<RealType, 3, 1> e(Pj - Pi), d(Pk - Pj), temp1, temp2;
    getAreaGradK(Pj, Pk, Pi, temp1);
    getAreaGradK(Pi, Pj, Pk, temp2);
    Hess = tensorProduct(temp1, temp2);
    Eigen::Matrix<RealType, 3, 3> auxMat;
    auxMat = tensorProduct(e, d);
    Hess += auxMat * 0.25;
    Hess.diagonal().array() += -0.25 * dotProduct(d, e);
    Hess *= -1. / getArea(Pi, Pj, Pk);
    getNormal(Pi, Pj, Pk, temp1);
    getCrossOp(temp1, auxMat);
    Hess += auxMat * 0.5;
}

template<typename RealType>
void getAreaIKDiffQuotient(const Eigen::Matrix<RealType, 3, 1> &Pi, const Eigen::Matrix<RealType, 3, 1> &Pj,
                           const Eigen::Matrix<RealType, 3, 1> &Pk, Eigen::Matrix<RealType, 3, 3> &Hik) {
    RealType H = 1e-8;
    for (int i = 0; i < 3; i++) {
        Eigen::Matrix<RealType, 3, 1> PiPlusH(Pi), grad1, grad2;
        PiPlusH[i] += H;
        getAreaGradK(PiPlusH, Pj, Pk, grad1);
        getAreaGradK(Pi, Pj, Pk, grad2);
        grad1 -= grad2;
        for (int j = 0; j < 3; j++)
            Hik(j, 1) = grad1[j] / H;
    }
}

template<typename RealType>
void getHessThetaKK(const Eigen::Matrix<RealType, 3, 1> &Pi, const Eigen::Matrix<RealType, 3, 1> &Pj,
                    const Eigen::Matrix<RealType, 3, 1> &Pk, Eigen::Matrix<RealType, 3, 3> &Hkk) {

    RealType areaSqr = getArea(Pi, Pj, Pk) * getArea(Pi, Pj, Pk);

    Eigen::Matrix<RealType, 3, 1> e(Pj - Pi), gradArea, normal;
    getAreaGradK(Pi, Pj, Pk, gradArea);
    getNormal(Pi, Pj, Pk, normal);

    Eigen::Matrix<RealType, 3, 3> mat1, mat2;
    getCrossOp(e, mat1);
    mat2 = tensorProduct(gradArea, normal);

    getWeightedMatrixSum(e.norm() / (4. * areaSqr), mat1, e.norm() / areaSqr, mat2, Hkk);
}

template<typename RealType>
void getHessThetaIK(const Eigen::Matrix<RealType, 3, 1> &Pi, const Eigen::Matrix<RealType, 3, 1> &Pj,
                    const Eigen::Matrix<RealType, 3, 1> &Pk, Eigen::Matrix<RealType, 3, 3> &Hik) {

    RealType area = getArea(Pi, Pj, Pk);
    RealType areaSqr = area * area;

    Eigen::Matrix<RealType, 3, 1> e(Pj - Pi), d(Pk - Pj), gradArea, normal;
    getAreaGradK(Pj, Pk, Pi, gradArea);
    getNormal(Pi, Pj, Pk, normal);

    Eigen::Matrix<RealType, 3, 3> mat1, mat2, mat3;
    mat1 = tensorProduct(e, normal);
    getCrossOp(d, mat2);
    getWeightedMatrixSum(1. / (2. * area * e.norm()), mat1, e.norm() / (4. * areaSqr), mat2, mat3);

    mat1 = tensorProduct(gradArea, normal);
    getWeightedMatrixSum(1., mat3, e.norm() / areaSqr, mat1, Hik);
}

template<typename RealType>
void getHessThetaJK(const Eigen::Matrix<RealType, 3, 1> &Pi, const Eigen::Matrix<RealType, 3, 1> &Pj,
                    const Eigen::Matrix<RealType, 3, 1> &Pk, Eigen::Matrix<RealType, 3, 3> &Hjk) {

    RealType area = getArea(Pi, Pj, Pk);
    RealType areaSqr = area * area;

    Eigen::Matrix<RealType, 3, 1> e(Pi - Pj), a(Pi - Pk), gradArea, normal;
    getAreaGradK(Pk, Pi, Pj, gradArea);
    getNormal(Pi, Pj, Pk, normal);

    Eigen::Matrix<RealType, 3, 3> mat1, mat2, mat3;
    mat1 = tensorProduct(e, normal);
    getCrossOp(a, mat2);
    getWeightedMatrixSum(1. / (2. * area * e.norm()), mat1, e.norm() / (4. * areaSqr), mat2, mat3);

    mat1 = tensorProduct(gradArea, normal);
    getWeightedMatrixSum(1., mat3, e.norm() / areaSqr, mat1, Hjk);
}

template<typename RealType>
void getHessThetaILeftPartI(const Eigen::Matrix<RealType, 3, 1> &Pi, const Eigen::Matrix<RealType, 3, 1> &Pj,
                            const Eigen::Matrix<RealType, 3, 1> &Pk, Eigen::Matrix<RealType, 3, 3> &HiLeft) {
    Eigen::Matrix<RealType, 3, 1> e(Pj - Pi), d(Pk - Pj), eNormalized(Pj - Pi), gradThetaK, temp;
    eNormalized.normalize();
    Eigen::Matrix<RealType, 3, 3> mat1, mat2, Refl;
    getThetaGradK(Pi, Pj, Pk, gradThetaK);
    getReflection(eNormalized, Refl);
    temp = Refl * d;
    mat1 = tensorProduct(temp, gradThetaK);
    getHessThetaIK(Pi, Pj, Pk, mat2);
    getWeightedMatrixSum(-1. / dotProduct(e, e), mat1, dotProduct(d, e) / dotProduct(e, e), mat2, HiLeft);
}

template<typename RealType>
void getHessThetaJLeftPartI(const Eigen::Matrix<RealType, 3, 1> &Pi, const Eigen::Matrix<RealType, 3, 1> &Pj,
                            const Eigen::Matrix<RealType, 3, 1> &Pk, Eigen::Matrix<RealType, 3, 3> &HjLeft) {
    Eigen::Matrix<RealType, 3, 1> e(Pj - Pi), d(Pk - Pj), eNormalized(Pj - Pi), gradThetaK, temp, thetak;
    eNormalized.normalize();
    Eigen::Matrix<RealType, 3, 3> mat1, mat2, Refl;
    getThetaGradK(Pi, Pj, Pk, gradThetaK);
    getReflection(eNormalized, Refl);
    temp = Refl * (d - e);
    mat1 = tensorProduct(temp, gradThetaK);
    getHessThetaJK(Pi, Pj, Pk, mat2);
    getWeightedMatrixSum(1. / dotProduct(e, e), mat1, dotProduct(d, e) / dotProduct(e, e), mat2, HjLeft);
}

template<typename RealType>
void getHessThetaII(const Eigen::Matrix<RealType, 3, 1> &Pi, const Eigen::Matrix<RealType, 3, 1> &Pj,
                    const Eigen::Matrix<RealType, 3, 1> &Pk, const Eigen::Matrix<RealType, 3, 1> &Pl,
                    Eigen::Matrix<RealType, 3, 3> &Hii) {
    Eigen::Matrix<RealType, 3, 3> temp;
    getHessThetaILeftPartI(Pi, Pj, Pk, Hii);
    getHessThetaILeftPartI(Pi, Pj, Pl, temp);
    Hii -= temp;
}

template<typename RealType>
void getHessThetaJI(const Eigen::Matrix<RealType, 3, 1> &Pi, const Eigen::Matrix<RealType, 3, 1> &Pj,
                    const Eigen::Matrix<RealType, 3, 1> &Pk, const Eigen::Matrix<RealType, 3, 1> &Pl,
                    Eigen::Matrix<RealType, 3, 3> &Hji) {
    Eigen::Matrix<RealType, 3, 1> edge(Pj - Pi), d(Pk - Pj), c(Pj - Pl);
    Eigen::Matrix<RealType, 3, 1> diff(d - edge), sum(c + edge);
    RealType eLengthSqr = edge.squaredNorm();

    //
    Eigen::Matrix<RealType, 3, 1> thetak, thetal, grad;
    getThetaGradK(Pi, Pj, Pk, thetak);
    getThetaGradK(Pj, Pi, Pl, thetal);
    getWeightedVectorSum(dotProduct(edge, d), thetak, -1. * dotProduct(edge, c), thetal, grad);

    Eigen::Matrix<RealType, 3, 3> Hjk, Hjl, tProduct;
    getHessThetaJK(Pi, Pj, Pk, Hjk);
    getHessThetaIK(Pj, Pi, Pl, Hjl);

    // Hess part
    getWeightedMatrixSumTransposed(dotProduct(edge, d) / eLengthSqr, Hjk, -1. * dotProduct(edge, c) / eLengthSqr, Hjl,
                                   Hji);

    tProduct = tensorProduct(grad, edge);
    Hji += tProduct * (-2. / (eLengthSqr * eLengthSqr));

    tProduct = tensorProduct(thetak, diff);
    Hji += tProduct * (1. / eLengthSqr);
    tProduct = tensorProduct(thetal, sum);
    Hji += tProduct * (-1. / eLengthSqr);
}

//// DIFFERENCE QUOTIENTS
//
//template<typename RealType>
//void getHiiDiffQuotient( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, const SmallVec3<RealType>& Pl, SmallMat33<RealType>& Hii ) {
//  RealType H = 1e-8;
//  for( int i = 0; i < 3; i++ ){
//    SmallVec3<RealType> PiPlusH( Pi ), grad1, grad2;
//    PiPlusH[i] += H;
//    getThetaGradI( PiPlusH, Pj, Pk, Pl, grad1 );
//    getThetaGradI( Pi, Pj, Pk, Pl, grad2 );
//    grad1 -= grad2;
//    for( int j = 0; j < 3; j++ )
//      Hii.set( j, i, grad1[j] / H );
//  }
//}
//
//template<typename RealType>
//void getHjjDiffQuotient( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, const SmallVec3<RealType>& Pl, SmallMat33<RealType>& Hjj ) {
//  RealType H = 1e-8;
//  for( int i = 0; i < 3; i++ ){
//    SmallVec3<RealType> PjPlusH( Pj ), grad1, grad2;
//    PjPlusH[i] += H;
//    getThetaGradJ( Pi, PjPlusH, Pk, Pl, grad1 );
//    getThetaGradJ( Pi, Pj, Pk, Pl, grad2 );
//    grad1 -= grad2;
//    for( int j = 0; j < 3; j++ )
//      Hjj.set( j, i, grad1[j] / H );
//  }
//}
//
//template<typename RealType>
//void getHijDiffQuotient( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, const SmallVec3<RealType>& Pl, SmallMat33<RealType>& Hij ) {
//  RealType H = 1e-10;
//  for( int i = 0; i < 3; i++ ){
//    SmallVec3<RealType> PjPlusH( Pj ), PjMinusH( Pj ), grad1, grad2;
//    PjPlusH[i] += H;
//    PjMinusH[i] -= H;
//    getThetaGradI( Pi, PjPlusH, Pk, Pl, grad1 );
//    getThetaGradI( Pi, PjMinusH, Pk, Pl, grad2 );
//    grad1 -= grad2;
//    for( int j = 0; j < 3; j++ )
//      Hij.set( j, i, grad1[j] / (2*H) );
//  }
//}
//
//template<typename RealType>
//void getHjiDiffQuotient( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, const SmallVec3<RealType>& Pl, SmallMat33<RealType>& Hji ) {
//  RealType H = 1e-8;
//  for( int i = 0; i < 3; i++ ){
//    SmallVec3<RealType> PiPlusH( Pi ), grad1, grad2;
//    PiPlusH[i] += H;
//    getThetaGradJ( PiPlusH, Pj, Pk, Pl, grad1 );
//    getThetaGradJ( Pi, Pj, Pk, Pl, grad2 );
//    grad1 -= grad2;
//    for( int j = 0; j < 3; j++ )
//      Hji.set( j, i, grad1[j] / H );
//  }
//}
//
//template<typename RealType>
//void getHikDiffQuotient( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, const SmallVec3<RealType>& Pl, SmallMat33<RealType>& Hik ) {
//  RealType H = 1e-8;
//  for( int i = 0; i < 3; i++ ){
//    SmallVec3<RealType> PiPlusH( Pi ), grad1, grad2;
//    PiPlusH[i] += H;
//    getThetaGradK( PiPlusH, Pj, Pk, grad1 );
//    getThetaGradK( Pi, Pj, Pk, grad2 );
//    grad1 -= grad2;
//    for( int j = 0; j < 3; j++ )
//      Hik.set( j, i, grad1[j] / H );
//  }
//}
//
//template<typename RealType>
//void getHkiDiffQuotient( const SmallVec3<RealType>& Pi, const SmallVec3<RealType>& Pj, const SmallVec3<RealType>& Pk, const SmallVec3<RealType>& Pl, SmallMat33<RealType>& Hki ) {
//  RealType H = 1e-8;
//  for( int i = 0; i < 3; i++ ){
//    SmallVec3<RealType> PkPlusH( Pk ), grad1, grad2;
//    PkPlusH[i] += H;
//    getThetaGradI( Pi, Pj, PkPlusH, Pl, grad1 );
//    getThetaGradI( Pi, Pj, Pk, Pl, grad2 );
//    grad1 -= grad2;
//    for( int j = 0; j < 3; j++ )
//      Hki.set( j, i, grad1[j] / H );
//  }
//}

//=============================================================================
#endif // GEOMETRYEIGEN_HH defined
//=============================================================================

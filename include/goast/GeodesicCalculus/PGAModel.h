// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Full discrete PGA model.
 * \author Heeren, Sassen
 * See Heeren et al. \cite HeZhRu18
 * 
 */
#ifndef PGAMODEL_H
#define PGAMODEL_H

#include <goast/Core.h>
#include <goast/Optimization.h>

#include "ElasticMean.h"
#include "GeodesicMean.h"
#include "GramianMatrix.h"
#include "NonlinearProjection.h"
#include "DiscreteGeodesicCalculus.h"

#include <unsupported/Eigen/SparseExtra> // For MatrixMarket IO


//======================================================================================================================
// GEODESIC (RE)SCALING
//======================================================================================================================

/**
 * \brief Geodesic scaling
 * \author Heeren
 *
 * Literally the same as doubleInterpolation in DiscreteGeodesicCalculus.h
 *
 * \note See Heeren et al. 2018 \cite HeZhRu18
 */
template<typename ConfiguratorType>
void geodesicScaling( const MeshTopologySaver& Topology,
                      const typename ConfiguratorType::VectorType& StartGeom,
                      const typename ConfiguratorType::VectorType& EndGeom,
                      const std::vector<int>& Mask,
                      const DeformationBase<ConfiguratorType>& W,
                      const OptimizationParameters<ConfiguratorType>& optPars,
                      int totalLength,
                      double rho,
                      typename ConfiguratorType::VectorType& initPath,
                      typename ConfiguratorType::VectorType& Shape,
                      bool quiet = true ){
  doubleInterpolation<ConfiguratorType>( Topology, StartGeom, EndGeom, Mask, W, optPars, totalLength, rho, initPath, Shape, quiet );
}

/**
 * \brief Geodesic rescaling
 * \author Heeren
 *
 * Inverse operation of doubleInterpolation: First inverting the convex combination than applying geodesic shooting.
 *
 *  \note See Heeren et al. 2018 \cite HeZhRu18
 */
template<typename ConfiguratorType>
void geodesicScalingInverse( const MeshTopologySaver& Topology,
                             const typename ConfiguratorType::VectorType& StartGeom,
                             const typename ConfiguratorType::VectorType& varGeom,
                             const std::vector<int>& Mask,
                             const DeformationBase<ConfiguratorType>& W,
                             const OptimizationParameters<ConfiguratorType>& optPars,
                             int totalLength,
                             double rho,
                             typename ConfiguratorType::VectorType& Shape,
                             bool quiet = true ){

  // get left and right shape
  int leftIdx = std::floor( rho * (totalLength-1) );
  double lambda = rho*(totalLength-1) - leftIdx;
  typename ConfiguratorType::VectorType intermedShape( varGeom );

  if( lambda < 0.5 ){
    if(!quiet) std::cerr << "Compute initialization for convex extrapolation..." << std::endl;
    int K = std::floor( 1. / lambda ) ;
    std::vector<typename ConfiguratorType::VectorType> tempPath;
    if(!quiet) std::cerr << "Compute " << K-1 << " integer steps of shooting..." << std::endl;
    integerExtrapolation<ConfiguratorType>( Topology,  StartGeom, varGeom, Mask, W, optPars, K-1, tempPath, true, quiet );
    intermedShape = tempPath[K];
  }

  // compute convex extrapolation
  if(!quiet) std::cerr << "Compute convex extrapolation for between lambda = " << lambda << std::endl;
  convexExtrapolation<ConfiguratorType>( Topology, StartGeom, varGeom, Mask, W, optPars, lambda, intermedShape, quiet );

  // integer shooting
  std::vector<typename ConfiguratorType::VectorType> Path;
  if(!quiet) std::cerr << "Compute " << totalLength - 2 << " integer steps of shooting..." << std::endl;
  integerExtrapolation<ConfiguratorType>( Topology,  StartGeom, intermedShape, Mask, W, optPars, totalLength - 2, Path, true, quiet );
  Shape = Path[totalLength - 1];
}

//======================================================================================================================
// GEODESIC REFLECTION
//======================================================================================================================

/**
 * \brief Compute geodesic reflection \f$ r \f$ of point \f$ p \f$ about mean \f$ s \f$, such that \f$(r,s,p)\f$ is 3-point geodesic
 * \author Heeren
 */
template<typename ConfiguratorType>
void geodesicReflection( const MeshTopologySaver& Topology,
                         const typename ConfiguratorType::VectorType& Mean,
                         const typename ConfiguratorType::VectorType& Point,
                         const std::vector<int>& Mask,
                         const DeformationBase<ConfiguratorType>& W,
                         const OptimizationParameters<ConfiguratorType>& optPars,
                         typename ConfiguratorType::VectorType& Reflection,
                         bool quiet = true ){
  GeodesicExtrapolation<ConfiguratorType> extrapolationOp( Topology, Mean, Point, W, optPars, quiet );
  if( Mask.size() > 0 )
    extrapolationOp.setBoundaryMask( Mask );
  std::vector<typename ConfiguratorType::VectorType> Path;
  extrapolationOp.execute( 1, Path, false );
  Reflection = Path[2];
}



//======================================================================================================================
// DISCRETE GEODESIC PGA MODEL
//======================================================================================================================

/**
 * \brief Class that represents a full discrete PGA model (see Heeren et al., \cite HeZhRu18)
 * \author Heeren, Sassen
 *
 * A discrete PGA model based on \f$ n \f$ input shapes \f$ s^1, ..., s^n \f$ is built in several steps:
 *
 * First, a discrete geodesic mean \f$ \bar s\f$ is computed via \f[ \bar s = \arg\min_s \sum_{i=1}^n \sum_{k=1}^K W[s_{k-1}^i, s_k^i] \f]
 * subject to \f$ s_0^i = s \f$ and \f$ s_K^i=s^i \f$ for \f$ i = 1, ..., n \f$ (see also documenation in GeodesicMean.h).
 *
 * The mean comes with a sequence of geodesic segments \f$ (s_0^i = \bar s, s_1^i, ..., s_K^i = s^i) \f$ connecting  \f$ \bar s \f$ with the ith input shape \f$ s^i \f$.
 *
 * In particular, we denote \f$ (s_1^1, s_1^2, ..., s_1^n) \f$ as the set of inner shapes that represent shape evaluations of discrete logarithms (of \f$ s^i \f$ wrt. \f$ \bar s \f$).
 *
 * Second, an approximative Gram matrix is assembled based on the discrete mean and the inner shapes (see e.g. documenation in GramianMatrix.h).
 *
 * Third, based on a SVD of the approximate Gram's matrix we compute weighted elastic averages of the inner shapes to determine discrete principal variations (pvs)
 *
 * Finally, we compute geodesic reflections of all pvs \f$ p_j \f$ about the mean \f$ p_0 := \bar s \f$ for \f$ j = 1, ..., J \f$,
 * where \f$ J \f$ is the number of variations to be considered.
 *
 * Note the ordering convention for the full set of pvs: \f$ (p_0, p_1, p_2, ..., p_J, p_{-1}, p_{-2},..., p_{-J}) \f$,
 * that means \f$ (p[j], p[0], p[j+J]) \f$ is a short discrete geodesic for \f$ j = 1,...,J \f$ and \f$ p[0] \f$ is the geodesic mean.
 */
template<typename ConfiguratorType, typename DeformationType>
class DiscretePGAModel {
    
protected:
    typedef typename ConfiguratorType::RealType RealType;

    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType MatrixType;
    typedef typename ConfiguratorType::TripletType TripletType;
    typedef std::vector<TripletType> TripletListType;

    const ParameterParser& _pparser;
    
    TriMesh _refMesh;
    MeshTopologySaver* _topolP;
    DeformationType* _approxSqrDistP;

    std::vector<int> const* _maskP;

    VectorType _inputShapes;
    VectorType _mean;
    VectorType _geodMeanAndSegments;
    
    VectorType _eigenvalues;
    Eigen::MatrixXd _eigenvectors;
    VectorType _principalVariations;
    
    const unsigned int _K, _numInputData;
    unsigned int _numPrincVar, _numLocalDofs;    
    
    //std::vector<Eigen::Ref<const VectorType> > _princVarsRefs;
    //std::vector<Eigen::Ref<const VectorType> > _inputShapesRefs;
    
    bool _initMeanSegs, _initPrincVars, _quiet;

public:
    // load discrete PGA model based on parameter parser
    DiscretePGAModel( const ParameterParser& pparser, bool initMeanAndSegments = true, bool initPrincVariations = true ) 
    : _pparser(pparser), _topolP(NULL), _approxSqrDistP(NULL), _maskP(NULL), _K( _pparser.getInt("K") ), _numInputData( _pparser.getInt("numOfInputData") ), _quiet(false) { 
        // initialize
        clear();
        
        // load precomputed data?
        if( _pparser.checkAndGetBool("loadMeanAndSegments") ){
            loadGeodesicMeanAndSegments( _pparser.getString("loadnameStem") );
        }
        else{        
          // load input data
          if( !_pparser.hasVariable("inputFilenameStem"))
            throw BasicException("DiscretePGAModel::DiscretePGAModel() no inputFilenameStem given in parser!");
          loadInputData(  pparser.getString("inputFilenameStem") );
        }
        
        // set W
        _approxSqrDistP = new DeformationType( *_topolP, _pparser.getDouble("bendWeight") );          
          
        // initialize geodesic mean and geodesic segments (by loading or computing)
        if( !initMeanAndSegments )
            return;
        
        // initialize mean and segments (if they have not been loaded before)
        if( !_pparser.checkAndGetBool("loadMeanAndSegments") ){
            // compute mean and segments            
            computeGeodesicMeanAndSegments( _pparser.getInt("numOfAlternatingStepsGeodMean"), _pparser.getInt("numOfFullStepsGeodMean") );
        }     
        else{
            if(!_quiet){ 
              std::cerr << "Check segments..." << std::endl;
              VectorType alpha = VectorType::Constant(_numInputData, 1. / _numInputData); 
              OptimizationParameters<ConfiguratorType> dummyParams;
              GeodesicMean<ConfiguratorType>(*_topolP, *_approxSqrDistP, _K, _inputShapes, alpha, _numInputData, dummyParams, true).checkIfAllSegmentsAreRelaxed( _geodMeanAndSegments );
            }
        }         
        
        // initialize principal variations (by loading or computing)
        if( !initPrincVariations )
            return;
        
        // load principal variations 
        _numPrincVar = _pparser.getInt("numPrincipalVariations");
        if( _pparser.checkAndGetBool("loadPrincipalVariations") )
            loadPrincipalVariations( _pparser.getString("loadnameStem") );
        else
            computePrincipalVariations( _numPrincVar );
        
        // double-check quality of principal variations
        if(!_quiet) {
          for (int j = 0; j < _numPrincVar; j++){
            RealType tempEnergy;
            WeightedGeod3Energy<ConfiguratorType>( *_approxSqrDistP, _principalVariations.segment( (1 + _numPrincVar + j) * _numLocalDofs, _numLocalDofs ), _principalVariations.segment( (j + 1) * _numLocalDofs, _numLocalDofs),  0.5 ).apply( _mean, tempEnergy );
            std::cerr << "Path energy  in " << j << "th reflection is " << tempEnergy << std::endl;
      
            VectorType tempGrad;
            WeightedGeod3Gradient<ConfiguratorType>( *_approxSqrDistP, _principalVariations.segment( (1 + _numPrincVar + j) * _numLocalDofs, _numLocalDofs ), _principalVariations.segment( (j + 1) * _numLocalDofs, _numLocalDofs),  0.5 ).apply( _mean, tempGrad );
            std::cerr << "Norm of path energy gradient in " << j << "th reflection is " << tempGrad.norm() << std::endl;            
          }  
        }

    }
    
    ~DiscretePGAModel(){
        if( _topolP )
            delete _topolP;
        if( _approxSqrDistP )
            delete _approxSqrDistP;
    }
    
    int getK() const {
        return _K;
    }
    
    int getNumInputShapes() const {
        return _numInputData;
    }
    
    int getNumPrincipalVariations() const {
        return _numPrincVar;
    }
    
    const VectorType& getMean() const {
        return _mean;
    }
    
    const VectorType& getGeodesicMeanAndSegments() const {
        return _geodMeanAndSegments;
    }
    
    const VectorType& getPrincipalVariations() const {
        return _principalVariations;
    }
    
    const TriMesh& getReferenceMesh() const {
        return _refMesh;
    }

    const VectorType &getInputShapes() const {
      return _inputShapes;
    }

    // compute geodesic mean \bar s and geodesic segments by
    // first applying some steps of alternting scheme, i.e. (a) fix mean and update segments in parallel (b) update mean as elastic average of inner shapes,
    // and finally relaxing the full energy with n * (K-1) + 1 free shapes
    void computeGeodesicMeanAndSegments( int numOfAlternatingSteps, int numOfFullSteps ) {
        
        auto t_start_total = std::chrono::high_resolution_clock::now();

        if( !_topolP )
            throw BasicException("DiscretePGAModel::computeGeodesicMeanAndSegments(): Topology has not been defined!");
        if( _numInputData == 0 )
            throw BasicException("DiscretePGAModel::computeGeodesicMeanAndSegments(): No input shapes!");
        if( _inputShapes.size() != _numInputData * _numLocalDofs )
            throw BasicException("DiscretePGAModel::computeGeodesicMeanAndSegments(): dimension of input shapes is wrong!");
        if( !_approxSqrDistP )
            throw BasicException("DiscretePGAModel::computeGeodesicMeanAndSegments(): approximative functional W has not been defined!");

       // define optimization parameters
       OptimizationParameters<ConfiguratorType> optParsElasticMean;
       optParsElasticMean.setGradientIterations( _pparser.getIntOrDefault("elastMeanGradDescIters", 0) );
       optParsElasticMean.setNewtonIterations( _pparser.getIntOrDefault("elastMeanNewtonIters", 0) );
       optParsElasticMean.setQuietMode( SHOW_ALL );
       bool elastOpQuiet = false;
       
       OptimizationParameters<ConfiguratorType> optParsGeodesicMean;
       optParsGeodesicMean.setGradientIterations( _pparser.getIntOrDefault("geodMeanGradDescIters", 0) );
       optParsGeodesicMean.setBFGSIterations( _pparser.getIntOrDefault("geodMeanBFGSIters" , 0) );
       optParsGeodesicMean.setNewtonIterations( _pparser.getIntOrDefault("geodMeanNewtonIters", 0) );
       optParsGeodesicMean.setInitializationScheme( _pparser.getIntOrDefault("geodMeanInitializationScheme", 4) );
       optParsGeodesicMean.setQuiet();       
       //optParsGeodesicMean.setQuietMode( SHOW_ALL ); 
       bool geodOpQuiet = false;
       
       // setup functionals
       VectorType alpha = VectorType::Constant(_numInputData, 1. / _numInputData);       
       ElasticMean<ConfiguratorType> elasticMeanOp(  *_topolP, *_approxSqrDistP, _inputShapes, alpha, _numInputData, optParsElasticMean, elastOpQuiet );
       GeodesicMean<ConfiguratorType> geodesicMeanOp(*_topolP, *_approxSqrDistP, _K, _inputShapes, alpha, _numInputData, optParsGeodesicMean, geodOpQuiet);
       
       // set boundary mask
       if( _maskP ){
         if(!_quiet) std::cerr << "Set bounday mask for geodesic mean computation" << std::endl;
         geodesicMeanOp.setBoundaryMask( *_maskP );
         elasticMeanOp.setBoundaryMask( *_maskP );
       }       
       
       // precompute elastic mean as initialization
       if(!_quiet) std::cerr << "\n\n====================================================================================================" << std::endl; 
       if(!_quiet) std::cerr << "COMPUTE ELASTIC MEAN AS INITIALIZATION..." << std::endl;
       auto t_start_elast = std::chrono::high_resolution_clock::now();
       _mean = _inputShapes.segment( 0, _numLocalDofs );
       elasticMeanOp.execute(_mean);
       auto t_end_elast = std::chrono::high_resolution_clock::now();
       if(!_quiet) std::cout << std::fixed << "Elastic mean computed in " << std::chrono::duration<double, std::ratio<1> >(t_end_elast - t_start_elast).count() << "seconds." << std::endl << std::endl;
       if(!_quiet) std::cerr << "===================================="<< std::endl << std::endl;

       // initialize geodesic segments
       int numTotalFreeShapes = _numInputData * (_K-1) + 1;
       _geodMeanAndSegments.resize( numTotalFreeShapes * _numLocalDofs );
       _geodMeanAndSegments.segment( 0, _numLocalDofs ) = _mean;
       for (int k = 0; k < _numInputData; k++)
            for (int j = 0; j < _K - 1; j++)
                _geodMeanAndSegments.block((1 + k * (_K - 1) + j)* _numLocalDofs, 0, _numLocalDofs, 1) = _inputShapes.block(k * _numLocalDofs, 0, _numLocalDofs, 1);

       // alternating optimization
       if( numOfAlternatingSteps > 0 ){
         auto t_start = std::chrono::high_resolution_clock::now();
         if(!_quiet) std::cerr << "\n\n====================================================================================================" << std::endl; 
         if(!_quiet) std::cerr << "PERFORM ALTERNATING SCHEME TO COMPUTE INITIALIZATION FOR GEODESIC MEAN..." << std::endl;
         geodesicMeanOp.executeAlternating( _geodMeanAndSegments, numOfAlternatingSteps );
         auto t_end = std::chrono::high_resolution_clock::now();   
         if(!_quiet) std::cout << std::fixed << "Alternating optim. done in " << std::chrono::duration<double, std::ratio<1> >(t_end - t_start).count() << "seconds." << std::endl << std::endl;
       }
   
       // full optimization
       if( numOfFullSteps > 0 ){
         auto t_start = std::chrono::high_resolution_clock::now();
         if(!_quiet) std::cerr << "\n\n====================================================================================================" << std::endl; 
         if(!_quiet) std::cerr << "COMPUTE GEODESIC MEAN..." << std::endl;
         optParsGeodesicMean.setGradientIterations( numOfFullSteps );
         optParsGeodesicMean.setBFGSIterations( 0 );
         optParsGeodesicMean.setNewtonIterations( 0 );
         //optParsGeodesicMean.setQuietMode( SHOW_ALL );
         geodesicMeanOp.execute( _geodMeanAndSegments );
         auto t_end = std::chrono::high_resolution_clock::now();   
         if(!_quiet) std::cout << std::fixed << "Full optim. done in " << std::chrono::duration<double, std::ratio<1> >(t_end - t_start).count() << "seconds." << std::endl << std::endl;
       }
       
       if(!_quiet) {
         std::cerr << "Check segments..." << std::endl;
         geodesicMeanOp.checkIfAllSegmentsAreRelaxed( _geodMeanAndSegments );
       }
       
       // get mean separately
       _mean = _geodMeanAndSegments.segment( 0, _numLocalDofs );
       
       //
       _initMeanSegs = true;

       auto t_end_total = std::chrono::high_resolution_clock::now();
       if(!_quiet) std::cout << std::fixed << "Total geodesic mean computation done in " << std::chrono::duration<double, std::ratio<1> >(t_end_total - t_start_total).count() << "seconds." << std::endl << std::endl;
    }

    
    //
    void computePrincipalVariations( int numPrincipalVariations ) {
        
        int numTotalFreeShapes = _numInputData * (_K-1) + 1;
        if( _geodMeanAndSegments.size() != numTotalFreeShapes * _numLocalDofs )
            throw BasicException("DiscretePGAModel::computePrincipalVariations: wrong size in geodesic segments!");
        if( _eigenvalues.size() == 0 )
            computeSVDGramMatrix();
        
        // #### Principal Variations ####
        if(!_quiet) std::cerr << "\n\n====================================================================================================" << std::endl;   
        if(!_quiet) std::cerr << "START COMPUTATION OF PRONCIPAL VARIATIONS." << std::endl; 
        auto t_start_pv = std::chrono::high_resolution_clock::now();        
        
        _numPrincVar = numPrincipalVariations;
        if(!_quiet) std::cerr << "Consider J = " << _numPrincVar << " variations." << std::endl; 
        bool outerQuiet = true;
        
        // set optimization params
        OptimizationParameters<ConfiguratorType> optParsPV;
        optParsPV.setGradientIterations( _pparser.getIntOrDefault("elastMeanPVGradDescIters", 0) );
        optParsPV.setBFGSIterations(_pparser.getIntOrDefault("elastMeanPVBFGSIters", 0));
        optParsPV.setNewtonIterations(_pparser.getIntOrDefault("elastMeanPVNewtonIters", 0));
        optParsPV.setQuiet();
        
        // set optimization params
        OptimizationParameters<ConfiguratorType> optParsReflec;
        optParsReflec.setGradientIterations(_pparser.getIntOrDefault("geodReflecGradDescIters", 0));
        optParsReflec.setBFGSIterations(_pparser.getIntOrDefault("geodReflecBFGSIters", 0));
        optParsReflec.setNewtonIterations(_pparser.getIntOrDefault("geodReflecNewtonIters", 0));
        optParsReflec.setQuiet();

        // principal variations ordering convention: {p_0, p_1, p_2, ..., p_J, p_{-1}, p_{-2}, ..., p_{-J}},
        // i.e. (p_j, p_0, p_{-j}) is geodesic for all j = 1, ..., J
        _principalVariations.resize( (2 * _numPrincVar + 1) * _numLocalDofs );
        _principalVariations.segment( 0, _numLocalDofs ) = _mean;     

        // get inner shapes of the geodesics, i.e. their points closes to the mean
        const int numGeodesicsDofs = _numLocalDofs * (_K - 1);
        std::vector<VectorType> innerShapes;
        innerShapes.reserve(_numInputData);
        for (int k = 0; k < _numInputData; k++)
            innerShapes.push_back( _geodMeanAndSegments.segment( (1 + k * (_K - 1)) * _numLocalDofs, _numLocalDofs) );
        
        // compute reflections
        if(!_quiet) std::cerr << "\nCompute reflections..." << std::endl;
        std::vector<VectorType> innerShapesReflected(_numInputData);
        
        //TODO change this!
        std::vector<int> dummyMask;
        if( _maskP )
            throw BasicException("DiscretePGAModel::computePrincipalVariations: bdry mask is not taken into account yet!");
        
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int k = 0; k < _numInputData; k++)
            geodesicReflection<ConfiguratorType>(*_topolP, _mean, innerShapes[k], dummyMask, *_approxSqrDistP, optParsReflec, innerShapesReflected[k], outerQuiet );


        
        if(!_quiet) std::cerr << "\nFinally start to compute " << _numPrincVar << " principal variations..." << std::endl;
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int j = 0; j < _numPrincVar; j++) {
            //if(!_quiet) std::cerr << "\n-----------------------------------------------------------------" << std::endl;
            if(!_quiet) std::cerr << "Computing principal variation " << j << std::endl;
            // Local cache for the principal Variation
            VectorType principalVariation( _mean );

            // Normalize s.t. the sum of coefficients is one
            VectorType alphaWeights = _eigenvectors.col( _numInputData - 1 - j );
            alphaWeights /= alphaWeights.sum();
            //if(!_quiet)  std::cerr << "weights = " << alphaWeights << std::endl;

            // Shapes to use for the elastic mean, as we may have to reflect some of the endpoints of the geodesics
            VectorType localShapes( _numInputData * _numLocalDofs );
            for (int k = 0; k < _numInputData; k++) {

                // If the sign of the coefficient is negative, we have to reflect the endpoint at the mean and change
                // the sign of the coefficient in the eigenvector
                if (alphaWeights[k] < 0) {
                    //std::cerr << "PV " << j << " reflecting pt " << k << std::endl;
                    localShapes.segment(k * _numLocalDofs, _numLocalDofs) = innerShapesReflected[k];
                    alphaWeights[k] *= -1;
                }
                else
                    localShapes.segment(k * _numLocalDofs, _numLocalDofs) = innerShapes[k];
            }

            //if(!_quiet)  std::cerr << "PV " << j << ": Computing elastic mean" << std::endl;
            // Compute the weighted elastic mean of the (possibly reflected) closest shapes to get the principal variation
            ElasticMean<ConfiguratorType> elasticMeanOp( *_topolP, *_approxSqrDistP, localShapes, alphaWeights, _numInputData, optParsPV, outerQuiet );
            if( _maskP )
              elasticMeanOp.setBoundaryMask( *_maskP );
            elasticMeanOp.execute(principalVariation);

            // store it in the global vector (cf ordering convention above!)
            _principalVariations.segment( (j + 1) * _numLocalDofs, _numLocalDofs) = principalVariation;
   
            //if(!_quiet) std::cerr << "Reflecting principal variation " << j << std::endl;
            VectorType reflectedPrincipalVariation( _mean );
            geodesicReflection<ConfiguratorType>(*_topolP, _mean, principalVariation, dummyMask, *_approxSqrDistP, optParsReflec, reflectedPrincipalVariation, outerQuiet );
            
            // store it in the global vector (cf ordering convention above!)
            _principalVariations.segment( (1 + _numPrincVar + j) * _numLocalDofs, _numLocalDofs ) = reflectedPrincipalVariation;
        }

        auto t_end_pv = std::chrono::high_resolution_clock::now();
        if(!_quiet) std::cout << std::fixed << "Principal variations done in " << std::chrono::duration<double, std::ratio<1> >(t_end_pv - t_start_pv).count() << "seconds." << std::endl;
        
        _initPrincVars = true;
    }
    
    // Project an unseen shape s onto the submanifold spanned by the principal variations
    // The projection is done in three steps:
    // (1) scaling, i.e. compute a geodesic scaling of s wrt. the mean \bar s - denoted by s_loc
    // (2) local projection of s_loc onto the convex polyhedron spannend by the pvs  - denoted by P_loc[s_loc]
    // (3) re-scaling of P_loc[s_loc] to get projection P[s] of s
    // The user can specify the number of pvs to be taken into account.
    // The result is given as the global projection P[s] as well as the (alpha) coefficients of P_loc[s_loc]
    void projectShape( const VectorType &unseenGeom, int numPrincVars, VectorType &Alpha, VectorType& globalProjection ) const {
        
        if( !_initMeanSegs || !_initPrincVars )
            throw BasicException("DiscretePGAModel::projectShape(): model has not been initialized yet!");
        
        if( unseenGeom.size() != _numLocalDofs )
            throw BasicException("DiscretePGAModel::projectShape(): argument has wrong size!");
        
        if( numPrincVars > _numPrincVar )
            throw BasicException("DiscretePGAModel::projectShape(): the model does not have so many dimensions!");
        
        // #### Principal Variations ####
        if(!_quiet) std::cerr << "\n\n====================================================================================================" << std::endl;   
        if(!_quiet) std::cerr << "START NONLINEAR PROJECTION USING J = " << numPrincVars << " PRINCIPAL VARIATIONS." << std::endl; 
        auto t_start_proj = std::chrono::high_resolution_clock::now();

        RealType tempE;
        _approxSqrDistP->applyEnergy( _mean, unseenGeom, tempE);
        if(!_quiet) std::cerr << "Initial global projection error = " << tempE << std::endl;
   
        // local projection optim. params
        OptimizationParameters<ConfiguratorType> optParsProjInner;
        optParsProjInner.setGradientIterations( _pparser.getIntOrDefault("localProjInnerGradDescIters", 10));
        optParsProjInner.setNewtonIterations( _pparser.getIntOrDefault("localProjInnerNewtonIters", 0) );
        optParsProjInner.setQuiet();
        //optParsProjInner.setQuietMode( SHOW_ALL );
   
        OptimizationParameters<ConfiguratorType> optParsProjOuter;
        optParsProjOuter.setGradientIterations( _pparser.getIntOrDefault("localProjOuterGradDescIters", 10));
        optParsProjOuter.setBFGSIterations( _pparser.getIntOrDefault("localProjOuterBFGSIters", 0) );
        optParsProjOuter.setQuietMode( SHOW_ALL );
 
        // compute scaling factor
        RealType rho = computeScalingFactor( unseenGeom, numPrincVars );

        if (std::isnan(rho) || rho > 1) { // If rho is nan, e.g. if we try to project the mean, or rho is really big
          // Set alphas to a constant value
          Alpha.resize(2 * numPrincVars);
          Alpha.setConstant(1. / (2*numPrincVars));

          // Return mean
          Alpha = _mean;

          return;
        }
        
        // set boundary mask
        std::vector<int> mask;
        if( _maskP )
            mask = *_maskP;
   
        // optim. params for scaling
        OptimizationParameters<ConfiguratorType> optParsScale;
        optParsScale.setGradientIterations( _pparser.getIntOrDefault("geodScalingGradDescIters", 10) );
        optParsScale.setBFGSIterations( _pparser.getIntOrDefault("geodScalingBFGSIters", 0) );
        optParsScale.setNewtonIterations( _pparser.getIntOrDefault("geodScalingNewtonIters", 0) );
        optParsScale.setInitializationScheme( _pparser.getIntOrDefault("geodScalingInitializationScheme", 4) );
        //optParsScale.setQuietMode( SHOW_ALL );
        optParsScale.setQuiet();   
   
        // scale geometry  
        VectorType geodesicPathToUnseenShape;
        if(!_quiet) std::cerr << "1) SCALING." << std::endl;
        if(!_quiet) std::cerr << "Scaling factor rho = " << rho << std::endl; 
        VectorType scaledGeom( _mean );
        geodesicScaling<ConfiguratorType>( *_topolP, _mean, unseenGeom, mask, *_approxSqrDistP, optParsScale, _K + 1, rho, geodesicPathToUnseenShape, scaledGeom, false );


        auto t_start_loc_proj = std::chrono::high_resolution_clock::now(); 
        if(!_quiet) std::cerr << "\n2) LOCAL PROJECTION." << std::endl; 
        
        VectorType localProjection( _mean );
        VectorType activePVs( _principalVariations.segment( 0, (1 + 2*numPrincVars) * _numLocalDofs) );
        NonlinearLocalProjection<ConfiguratorType> localProjOp( *_topolP, *_approxSqrDistP, activePVs, scaledGeom, optParsProjOuter, optParsProjInner, _quiet );
        if (Alpha.size() != 2 * numPrincVars) {
          Alpha.resize( 2 * numPrincVars );
          for ( int i = 0; i < 2 * numPrincVars; i++ )
            Alpha[i] = 1. / (2. * numPrincVars + 1.);
        }
        
        if(!_quiet) std::cerr << "Initial weights = " << Alpha << std::endl;
        localProjOp.execute( Alpha, localProjection );
        if(!_quiet) std::cerr << "Optimal weights = " << Alpha << std::endl;
        if(!_quiet) std::cerr << "sum = " << Alpha.sum() << std::endl << std::endl;

        _approxSqrDistP->applyEnergy( scaledGeom, localProjection, tempE);
        if(!_quiet) std::cerr << "Local projection error = " << tempE << std::endl;
        
        auto t_end_loc_proj = std::chrono::high_resolution_clock::now();
        if(!_quiet) std::cout << std::fixed << "Local projection done in " << std::chrono::duration<double, std::ratio<1> >(t_end_loc_proj - t_start_loc_proj).count() << "seconds." << std::endl;
   
        // rescale 
        if(!_quiet) std::cerr << "3) RESCALING." << std::endl;
        if(!_quiet) std::cerr << "Inverse scaling factor 1/rho = " << 1./rho << std::endl; 
        globalProjection = localProjection;
        geodesicScalingInverse<ConfiguratorType>( *_topolP, _mean, localProjection, mask, *_approxSqrDistP, optParsScale, _K + 1, rho, globalProjection );
        _approxSqrDistP->applyEnergy( globalProjection, unseenGeom, tempE);
        if(!_quiet) std::cerr << "Global projection error = " << tempE << std::endl;
   
        auto t_end_proj = std::chrono::high_resolution_clock::now();
        if(!_quiet) std::cout << std::fixed << "Nonlinear projection done in " << std::chrono::duration<double, std::ratio<1> >(t_end_proj - t_start_proj).count() << "seconds." << std::endl;
    }          

    
    //
    void clear() {
        if (! _quiet ) std::cerr << "Clear discrete PGA model." << std::endl;
        if( _topolP )
            delete _topolP;
        _refMesh.clear();
        
        _inputShapes.resize(0);
        _mean.resize(0);
        _geodMeanAndSegments.resize(0);
        
        _eigenvalues.resize(0);
        _eigenvectors.resize(0,0);
        _principalVariations.resize(0);

        _numPrincVar = 0;
        _numLocalDofs = 0;
        
        _initMeanSegs = false;
        _initPrincVars = false;
    }
        
    //
    void loadInputData( const std::string& inputFilenameStem ) {
        
      if( !_quiet ) std::cerr << "Load input data from disk." << std::endl;
      TriMesh mesh;
      VectorType Geometry;

      if( _numLocalDofs > 0 )
        _inputShapes.resize( _numInputData * _numLocalDofs );
      
      for (int k = 0; k < _numInputData; k++) {            
            std::ostringstream filename;
            filename << inputFilenameStem << k << ".ply";
            if ( !OpenMesh::IO::read_mesh( mesh, filename.str() ) )
                std::cerr << "Read error for " << filename.str() << std::endl;             
            
            // initialize topology if this hasn't been done before
            if( initTopology( mesh ) )
                _inputShapes.resize( _numInputData * _numLocalDofs );
            
            // 
            getGeometry(mesh, Geometry);
            _inputShapes.segment(k * _numLocalDofs, _numLocalDofs) = Geometry;
      }
    }
    
    // save geodesic mean and segments
    void saveGeodesicMeanAndSegments( const std::string& savenameStem ) const {
        
        if( !_initMeanSegs )
            throw BasicException("DiscretePGAModel::saveGeodesicMeanAndSegments(): geodesic mean and segments have not been computed yet!");
        if( _numInputData == 0 )
            throw BasicException("DiscretePGAModel::saveGeodesicMeanAndSegments(): input data empty!");
        if( _numLocalDofs == 0 )
            throw BasicException("DiscretePGAModel::saveGeodesicMeanAndSegments(): no local dofs!");
        
        const int numOfTotalFreeShapes = 1 + _numInputData * (_K-1);
        const int numOfTotalDofs  = _numLocalDofs * numOfTotalFreeShapes;
        
        if( _geodMeanAndSegments.size() != numOfTotalDofs )
            throw BasicException("DiscretePGAModel::saveGeodesicMeanAndSegments(): wrong size in geodesic segments!");
        
        TriMesh outputMesh ( _refMesh );

        // save geodesic mean
        std::ostringstream meanName;
        meanName << savenameStem << "_geodMean.ply";
        setGeometry( outputMesh, _mean );
        OpenMesh::IO::write_mesh(outputMesh, meanName.str() );
        
        // save segments
        for( int i = 0; i < _numInputData; i++ ){
            // mean shape
            std::ostringstream startName;
            startName << savenameStem << "_segment" << i << "_shape0.ply";
            setGeometry( outputMesh, _mean );
            OpenMesh::IO::write_mesh(outputMesh, startName.str() );
            
            // segments            
            for( int j = 0; j < _K - 1; j++ ){
              std::ostringstream saveName;
              saveName << savenameStem << "_segment" << i << "_shape" << j+1 << ".ply";
              setGeometry( outputMesh, _geodMeanAndSegments.segment( (1 + i * (_K-1) + j) * _numLocalDofs , _numLocalDofs) );
              OpenMesh::IO::write_mesh(outputMesh, saveName.str() );
            }
            
            // input shape
            std::ostringstream inputName;
            inputName << savenameStem << "_segment" << i << "_shape" << _K << ".ply";
            setGeometry( outputMesh, _inputShapes.segment(i * _numLocalDofs, _numLocalDofs) );
            OpenMesh::IO::write_mesh(outputMesh, inputName.str() );
        }
    }
    
    //
    void loadGeodesicMeanAndSegments( const std::string& loadnameStem ) {
        
       // load mean
       TriMesh loadMesh;       
       if(!_quiet) std::cerr << "Load mean from disk..." << std::endl;
       std::ostringstream loadMean;
       loadMean << loadnameStem << "_geodMean.ply";
       OpenMesh::IO::read_mesh( loadMesh, loadMean.str() );
       getGeometry( loadMesh, _mean );       
       
       // check whether mean is consistent with topology
       initTopology( loadMesh );
       
       // resize 
       _inputShapes.resize( _numInputData * _numLocalDofs );
       _geodMeanAndSegments.resize( (1 + (_K-1)*_numInputData) * _numLocalDofs );
       _geodMeanAndSegments.segment(0, _numLocalDofs ) = _mean;
       
       // load segments
       VectorType loadGeom;
       if(!_quiet) std::cerr << "Load geodesic segments from disk!" << std::endl;
       for( int k = 0; k < _numInputData; k++ ){
         for( int i = 0; i < _K-1; i++ ){
           std::ostringstream loadSegment;
           loadSegment << loadnameStem << "_segment" << k << "_shape" << i+1 << ".ply";
           OpenMesh::IO::read_mesh( loadMesh, loadSegment.str() );
           getGeometry( loadMesh, loadGeom );
           _geodMeanAndSegments.segment( (1 + k*(_K-1) + i) * _numLocalDofs , _numLocalDofs ) = loadGeom;
         }
       }
       
       // load input data
       if(!_quiet) std::cerr << "Load input data!" << std::endl;       
       for( int k = 0; k < _numInputData; k++ ){
           std::ostringstream loadInputShape;
           loadInputShape << loadnameStem << "_segment" << k << "_shape" << _K << ".ply";
           OpenMesh::IO::read_mesh( loadMesh, loadInputShape.str() );
           getGeometry( loadMesh, loadGeom );
           _inputShapes.segment( k * _numLocalDofs, _numLocalDofs ) = loadGeom;
       }
       
       //
       _initMeanSegs = true;
    }
    
    
    //
    void savePrincipalVariations( const std::string& savenameStem ) const {
        
        if( !_initPrincVars )
           throw BasicException("DiscretePGAModel::savePrincipalVariations(): principal variations have not been initialized!");
        
        if( _principalVariations.size() != (1 + 2 * _numPrincVar) * _numLocalDofs )
            throw BasicException("DiscretePGAModel::savePrincipalVariations(): wrong size of principal variations!");
        
        // principal variations ordering convention: {p_0, p_1, p_2, ..., p_J, p_{-1}, p_{-2}, ..., p_{-J}},
        TriMesh outputMesh( _refMesh );
        for (int j = 0; j < _numPrincVar; j++) {     
            std::ostringstream pvName;
            pvName << savenameStem << "_posPV" <<j << ".ply";
            setGeometry( outputMesh, _principalVariations.segment( (1 + j) * _numLocalDofs, _numLocalDofs) );
            OpenMesh::IO::write_mesh( outputMesh, pvName.str() );

            std::ostringstream pvNameReflec;
            pvNameReflec << savenameStem << "_negPV" <<j << ".ply";            
            setGeometry( outputMesh, _principalVariations.segment( (1 + _numPrincVar + j) * _numLocalDofs, _numLocalDofs) );
            OpenMesh::IO::write_mesh( outputMesh, pvNameReflec.str() );
        }
    }
    
    //
    void loadPrincipalVariations( const std::string& loadnameStem ) {
        
       if( !_initMeanSegs )
           throw BasicException("DiscretePGAModel::loadPrincipalVariations(): pga model has not been initialized!");
        
       TriMesh loadMesh;
       VectorType loadGeom;       
       _principalVariations.resize( (1 + 2 * _numPrincVar) * _numLocalDofs );
       _principalVariations.segment( 0, _numLocalDofs) = _mean;
       
       // load positive variations
       if(!_quiet) std::cerr << "Load positive pcs from disk!" << std::endl;
       for( int j = 0; j < _numPrincVar; j++ ){
           std::ostringstream pvName;
            pvName << loadnameStem << "_posPV" << j << ".ply";
           OpenMesh::IO::read_mesh( loadMesh, pvName.str() );
           getGeometry( loadMesh, loadGeom );
           _principalVariations.segment( (j + 1) * _numLocalDofs, _numLocalDofs) = loadGeom;
       }
       
       // load negative variations
       if(!_quiet) std::cerr << "Load negative pcs from disk!" << std::endl;
       for( int j = 0; j < _numPrincVar; j++ ){
           std::ostringstream pvName;
           pvName << loadnameStem << "_negPV" << j << ".ply";
           OpenMesh::IO::read_mesh( loadMesh, pvName.str() );
           getGeometry( loadMesh, loadGeom );
           _principalVariations.segment( (1 + _numPrincVar + j) * _numLocalDofs, _numLocalDofs ) = loadGeom;
       }
       
       _initPrincVars = true;
    }
            
protected:
    //
    bool initTopology( const TriMesh& refMesh ) {
        // check whether topology is set
        if( _topolP ){
            if( _topolP->getNumVertices() == refMesh.n_vertices() )
                return false;
            clear();
        }
        
        //
        if (! _quiet ) std::cerr << "Reset topology in discrete PGA model." << std::endl;
        _refMesh = refMesh;
        _topolP = new MeshTopologySaver( _refMesh );
        _numLocalDofs = 3 * _topolP->getNumVertices();
        
        return true;
    }
    
       
    //
    void computeSVDGramMatrix() {
        
        int numTotalFreeShapes = _numInputData * (_K-1) + 1;
        if( _geodMeanAndSegments.size() != numTotalFreeShapes * _numLocalDofs )
            throw BasicException("DiscretePGAModel::computeSVDGramMatrix(): wrong size in geodesic segments!");
        
        if(!_quiet) std::cerr << "\n\n====================================================================================================" << std::endl;   
        if(!_quiet) std::cerr << "START GRAM MATRIX COMPUTATION." << std::endl;   
        auto t_start_gram = std::chrono::high_resolution_clock::now();
        OptimizationParameters<ConfiguratorType> dummyParams;
        GramianMatrixAssembler<ConfiguratorType> gma( *_topolP, *_approxSqrDistP, _inputShapes, _mean, _numInputData, _K, dummyParams, false );
        if( _maskP )
          gma.setBoundaryMask( *_maskP );

        MatrixType gramianMatrix;
        if(!_quiet) std::cerr << "Assemble Gram matrix" << std::endl;
        
        gma.assemble( _geodMeanAndSegments.segment(_numLocalDofs, (numTotalFreeShapes-1) * _numLocalDofs), gramianMatrix );
        Eigen::MatrixXd denseGramianMatrix( gramianMatrix );
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es( denseGramianMatrix );
        auto t_end_gram = std::chrono::high_resolution_clock::now();
        if(!_quiet) std::cout << "Assembled gramian matrix in " << std::chrono::duration<double, std::ratio<1> >(t_end_gram - t_start_gram).count() << "seconds."  << std::endl;       
        
        //
        _eigenvalues = es.eigenvalues();
        _eigenvectors = es.eigenvectors();
        if(!_quiet) std::cerr << "Eigenvalues = " << _eigenvalues << std::endl;
    }

    //
    RealType computeScalingFactor( const VectorType& unseenGeom, int numPrincVars ) const {
        if(!_quiet) std::cerr << "Compute scaling factor" << std::endl;
        RealType kappa = _pparser.getDoubleOrDefault("kappa", 0.5);
        RealType rho = FLT_MAX;
        
        if( numPrincVars > _numPrincVar )
            throw BasicException("DiscretePGAModel::computeScalingFactor(): the model does not have so many dimensions!");
   
        // compare against all active pcs
        RealType tempE;
        for( int i = 0; i < numPrincVars; i++ ){       
          _approxSqrDistP->applyEnergy ( _mean, _principalVariations.segment( (i + 1) * _numLocalDofs, _numLocalDofs), tempE);
          if(!_quiet) std::cerr << "dist[ mean, pc_" << i << "] approx " << std::sqrt( tempE ) << std::endl;
          if( rho > std::sqrt( tempE ) )
            rho = std::sqrt( tempE );
        }
        _approxSqrDistP->applyEnergy ( _mean, unseenGeom, tempE);
        if(!_quiet) std::cerr << "dist[ mean, s]  approx " << std::sqrt( tempE ) << std::endl;
        return rho * kappa / std::sqrt( tempE );
    }
    
    
};

#endif

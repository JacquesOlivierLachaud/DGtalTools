/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
/**
 * @file 3dFastCorrectedNormalCurrent.cpp
 * @ingroup surfaceTools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2017/06/01
 *
 *
 * This file is part of the DGtalTools.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#define GL_SILENCE_DEPRECATION
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

//#include "EstimatorHelpers.h"
#include "DGtal/helpers/Shortcuts.h"
#include "DGtal/helpers/ShortcutsGeometry.h"
#include "EstimatorHelpers.h"
#include "FastCorrectedNormalCurrent.h"


using namespace DGtal;
using namespace functors;
using namespace std;

/**
 @page Doc3DCorrectedNormalCurrent 3d-fcnc

 @brief  Computes and visualizes the 3d corrected normal current of digital surfaces.

 @b Usage:  3d-fcnc -i file.vol

 @b Allowed @b options @b are:

 @code
  -h [ --help ]                         display this message
  -i [ --input ] arg                    .vol file
 @endcode

 @b Example:
 @see
 @ref 3d-fcnc.cpp
 @ref Doc3DCorrectedNormalCurrent
 */

///////////////////////////////////////////////////////////////////////////////

template <typename Scalar>
GradientColorMap<Scalar> 
getErrorColorMap( Scalar max )
{
  GradientColorMap<Scalar> gradcmap( 0.0, max );
  gradcmap.addColor( Color( 255, 255, 255 ) );
  gradcmap.addColor( Color( 255,   0,   0 ) );
  gradcmap.addColor( Color( 0,   0,   0 ) );
  return gradcmap;
}

template <typename Colors>
Colors getSimplifiedColorMap( const Colors& colors, unsigned char lost_mask = 0x03 )
{
  auto    new_colors = colors;
  unsigned char mask = 0xff ^ lost_mask;
  for ( Color& c : new_colors )
    {
      c.red  ( c.red()   & mask );
      c.green( c.green() & mask );
      c.blue ( c.blue()  & mask );
      c.alpha( c.alpha() & mask );
    }
  return new_colors;
}

/**
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam( std::string param )
{
  trace.error() << " Parameter: " << param << " is required.";
  trace.info() << std::endl;
}

//namespace po = DGtal::po;

int main( int argc, char** argv )
{  
  typedef Z3i::Space                            Space;
  typedef Z3i::KSpace                           KSpace;
  typedef Shortcuts<KSpace>                     SH;
  typedef ShortcutsGeometry<KSpace>             SHG;
  typedef EstimatorHelpers<KSpace>              EH;
  typedef SH::IdxDigitalSurface                 Surface;
  typedef Surface::DigitalSurfaceContainer      Container;
  typedef SH::RealVector                        RealVector;
  typedef SH::BinaryImage                       BinaryImage;
  typedef SH::ImplicitShape3D                   ImplicitShape3D;
  typedef FastCorrectedNormalCurrent<Container> Current;
  typedef Current::Vertex                       Vertex;
  
  // parse command line ----------------------------------------------
  po::options_description general_opt( "Allowed options are" );
  general_opt.add_options()
    ( "help,h", "display this message" )
    ( "m-coef", po::value<double>()->default_value( 3.0 ), "the coefficient k that defines the radius of the ball used in measures, that is r := k h^b" )
    ( "m-pow", po::value<double>()->default_value( 0.5 ), "the coefficient b that defines the radius of the ball used in measures, that is r := k h^b" );
  
  EH::optionsImplicitShape   ( general_opt );
  EH::optionsDigitizedShape  ( general_opt );
  EH::optionsVolFile         ( general_opt );
  EH::optionsNoisyImage      ( general_opt );
  EH::optionsNormalEstimators( general_opt );
  general_opt.add_options()
    ( "quantity,Q", po::value<std::string>()->default_value( "Mu1" ), "the quantity that is evaluated in Mu0|Mu1|Mu2|MuOmega|H|G|Omega|MuXY|HII|GII, with H := Mu1/(2Mu0), G := Mu2/Mu0, Omega := MuOmega/sqrt(Mu0), MuXY is the anisotropic curvature tensor and HII and GII are the mean and gaussian curvatures estimated by II." )
    ( "anisotropy", po::value<std::string>()->default_value( "NAdd" ), "tells how is symmetrized the anisotropic measure mu_XY, in Mult|Add|NMult|NAdd: Mult forces symmetry by M*M^t, Add forces symmetry by 0.5*(M+M^t)+NxN, NMult and NAdd normalized by the area.")
    ( "crisp,C", "when specified, when computing measures in a ball, do not approximate the relative intersection of cells with the ball but only consider if the cell centroid is in the ball (faster by 30%, but less accurate)." )
    ( "interpolate,I", "when specified, it interpolate the given corrected normal vector field and uses the corresponding measures." );
  EH::optionsDisplayValues   ( general_opt );
  general_opt.add_options()
    ( "zero-tic", po::value<double>()->default_value( 0.0 ), "adds a black band around zero of given thickness in colormaps." );

  //#endif
  general_opt.add_options()
    ( "error", po::value<std::string>()->default_value( "error.txt" ), "the name of the output file that sum up l2 and loo errors in estimation." );
  general_opt.add_options()
    ( "max-error", po::value<double>()->default_value( 0.2 ), "the error value corresponding to black." );
  general_opt.add_options()
    ( "output,o", po::value<std::string>()->default_value( "none" ), "the basename for output obj files or none if no output obj is wanted." );
  
  po::variables_map vm;
  bool parseOK = EH::args2vm( general_opt, argc, argv, vm );
  bool neededArgsGiven=true;

  if ( vm.count( "polynomial-list" ) )
    {
      trace.info() << "List of predefined polynomials:" << std::endl;
      auto L = SH::getPolynomialList();
      for ( auto p : L ) {
	trace.info() << "  " << p.first << " -> " << p.second << std::endl;
      }
    }

  if (parseOK && ( ! vm.count("polynomial") ) && ( ! vm.count( "input" ) ) ) {
    missingParam("--polynomial or --input");
    neededArgsGiven = false;
  }
  
  if ( !neededArgsGiven || !parseOK || vm.count("help") || argc <= 1 )
    {
      trace.info()<< "Builds the 3d corrected normal currents" <<std::endl
                  << general_opt << "\n"
                  << "Basic usage: "<<std::endl
		  << "\t 3d-fcnc -e II -t 3 -r 3 -Q H --tics 1 --minValue -0.3 --maxValue 0.3 --polynomial-list -p goursat -g 0.5 -V Measure --m-coef 3 -N 0.2" << std::endl
                  << "\t 3d-fcnc -p \"3x^2+5y^2+7z^2-1\" "<<std::endl
                  << "\t 3d-fcnc -i \"file.vol\" "<<std::endl
		  << "\t 3d-fcnc -e II -t 3 -r 3 -Q H --tics 1 --minValue -0.3 --maxValue 0.3 -V Measure --m-coef 3  -i ~/Images/3d/vol/fandisk-128.vol -m 0 -M 255 -N 0.4 -g 0.25" << std::endl
                  << "\t 3d-fcnc --polynomial-list  // to get the list of predefined polynomials. "<<std::endl
                  << std::endl;
      return 0;
    }
  auto quantity   = vm[ "quantity"   ].as<std::string>();
  auto anisotropy = vm[ "anisotropy" ].as<std::string>();
  std::vector< std::string > quantities = { "Mu0", "Mu1", "Mu2", "MuOmega", "H", "G", "Omega", "MuXY", "HII", "GII" };
  if ( std::count( quantities.begin(), quantities.end(), quantity ) == 0 ) {
    trace.error() << "Quantity should be in Mu0|Mu1|Mu2|MuOmega|H|G|Omega|MuXY|HII|GII.";
    trace.info()  << " I read quantity=" << quantity << std::endl;
    return 0;
  }

  // Digital space.
  KSpace                      K;
  unsigned int                nb = 0;
  CountedPtr<BinaryImage>     bimage( nullptr );
  CountedPtr<ImplicitShape3D> shape ( nullptr );
  auto params = SH::defaultParameters() | SHG::defaultParameters();
  trace.beginBlock( "Make Shape" );
  // Generic parameters.
  params( "gridstep",          vm[ "gridstep" ].as<double>() );
  params( "noise",             vm[ "noise"    ].as<double>() );
  params( "surfelAdjacency",          0 ); // 0:interior
  params( "nbTriesToFindABel",   100000 ); // number of tries in method Surfaces::findABel
  params( "surfaceComponents", "AnyBig" ); // "AnyBig"|"All"
  params( "projectionMaxIter",       20 ); // the maximum number of iter for the projection.
  params( "projectionAccuracy",  0.0001 ); // zero-proximity stop crit. during projection.
  params( "projectionGamma",        0.5 ); // the displacement coef. of the projection.
  params( "verbose",                  1 );
  params( "t-ring",            vm[ "trivial-ring" ].as<double>() );
  params( "kernel",            vm[ "kernel"       ].as<string>() );
  params( "R-radius",          vm[ "R-radius"     ].as<double>() );
  params( "r-radius",          vm[ "r-radius"     ].as<double>() );
  params( "alpha",             vm[ "alpha"        ].as<double>() );
  params( "surfelEmbedding",   vm[ "embedding"    ].as<int>()    );
  trace.info() << params << std::endl;
  if ( vm.count( "polynomial" ) )
    {
      // Fill useful parameters
      params( "polynomial", vm[ "polynomial" ].as<string>() );
      params( "minAABB",    vm[ "minAABB"    ].as<double>() );
      params( "maxAABB",    vm[ "maxAABB"    ].as<double>() );
      params( "offset",     5.0 );
      shape        = SH::makeImplicitShape3D( params );
      K            = SH::getKSpace( params );
      auto dshape  = SH::makeDigitizedImplicitShape3D( shape, params );
      bimage       = SH::makeBinaryImage( dshape, params );
    }
  else if ( vm.count( "input" ) )
    {
      // Fill useful parameters
      params( "thresholdMin", vm[ "thresholdMin" ].as<int>() );
      params( "thresholdMax", vm[ "thresholdMax" ].as<int>() );
      params( "closed",       vm[ "closed" ].as<int>() );
      auto volfile = vm[ "input" ].as<string>();
      bimage       = SH::makeBinaryImage( volfile, params );
      K            = SH::getKSpace( bimage, params );
    }
  auto size    = K.upperBound() - K.lowerBound();
  trace.info() << "- Domain size is " << ( size[ 0 ] + 1 )
	       << " x " << ( size[ 1 ] + 1 )
	       << " x " << ( size[ 2 ] + 1 ) << std::endl;

  std::for_each( bimage->cbegin(), bimage->cend(),
		 [&nb] ( bool v ) { nb += v ? 1 : 0; } );
  trace.info() << "- digital shape has " << nb << " voxels." << std::endl;
  
  auto sembedder   = SH::getSCellEmbedder( K );
  auto embedder    = SH::getCellEmbedder( K );
  auto surface     = SH::makeDigitalSurface( bimage, K, params );
  auto idx_surface = SH::makeIdxDigitalSurface( surface );
  if ( surface == 0 ) {
    trace.error() << "- surface is empty (either empty or full volume). ";
    trace.info()  << std::endl;
    trace.endBlock();
    return 1;
  }
  trace.info() << "- surface has " << surface->size()<< " surfels." << std::endl;
  trace.endBlock();

  const auto minValue         = vm[ "minValue" ].as<double>();
  const auto maxValue         = vm[ "maxValue" ].as<double>();
  const auto colormap_name    = vm[ "colormap" ].as<string>();
  const auto zt               = vm[ "zero-tic" ].as<double>();
  const double vthickness      = 0.1;
  const double vlength         = 1.5;
  auto outputfile    = vm[ "output"   ].as<string>();
  auto estimator     = vm[ "estimator" ].as<string>();
  const double h     = vm[ "gridstep"  ].as<double>();
  const double mcoef = vm[ "m-coef"    ].as<double>();
  const double mpow  = vm[ "m-pow"     ].as<double>();
  const double mr    = mcoef * pow( h, mpow );
  SH::Scalars       displayed_values;
  SH::Scalars       expected_values;
  SH::Scalars       measured_values;
  SH::RealVectors   expected_normals;
  SH::RealVectors   measured_normals;
  Statistic<double> stat_expected_curv;
  Statistic<double> stat_measured_curv;
  double time_curv_ground_truth  = 0.0;
  double time_curv_estimations   = 0.0;
  double time_mu_estimations     = 0.0;
  double time_normal_estimations = 0.0;
  
  trace.beginBlock( "Compute surfels" );
  params( "surfaceTraversal", "Default" );
  const auto     surfels = SH::getSurfelRange( surface, params );
  params( "surfaceTraversal", "DepthFirst" );
  //params( "surfaceTraversal", "Default" ); //< Force default (bug in II ?)
  const auto dft_surfels = SH::getSurfelRange( surface, params );
  const auto       match = SH::getRangeMatch ( surfels, dft_surfels );
  trace.endBlock();

  // Compute true and/or estimated normal vectors
  if ( vm.count( "polynomial" ) )
    {
      trace.beginBlock( "Compute true normals" );
      expected_normals = SHG::getNormalVectors( shape, K, dft_surfels, params );
      trace.endBlock();
    }
  trace.beginBlock( "Estimate normals" );
  if ( estimator == "True" && vm.count( "polynomial" ) )
    measured_normals = expected_normals;
  else if ( estimator == "CTrivial" )
    measured_normals = SHG::getCTrivialNormalVectors( surface, dft_surfels, params );
  else if ( estimator == "VCM" )
    measured_normals = SHG::getVCMNormalVectors( surface, dft_surfels, params );
  else if ( estimator == "II" ) {
    measured_normals = SHG::getIINormalVectors( bimage, dft_surfels, params );
    auto oriented_normals = SHG::getCTrivialNormalVectors( surface, dft_surfels, params );
    SHG::orientVectors( measured_normals, oriented_normals );
  }
  else // default is "Trivial"
    measured_normals = SHG::getTrivialNormalVectors( K, dft_surfels );
  time_normal_estimations = trace.endBlock();

  // Compute true or estimated curvatures
  if ( vm.count( "polynomial" ) )
    {
      trace.beginBlock( "Compute true curvatures" );
      expected_values = ( ( quantity == "H" ) || ( quantity == "Mu1" )
			  || ( quantity == "HII" ) )
	? SHG::getMeanCurvatures    ( shape, K, dft_surfels, params )
	: SHG::getGaussianCurvatures( shape, K, dft_surfels, params );
      stat_expected_curv.addValues( expected_values.cbegin(), expected_values.cend() );
      stat_expected_curv.terminate();
      trace.info() << "- truth curv: avg = " << stat_expected_curv.mean() << std::endl;
      trace.info() << "- truth curv: min = " << stat_expected_curv.min() << std::endl;
      trace.info() << "- truth curv: max = " << stat_expected_curv.max() << std::endl;
      time_curv_ground_truth = trace.endBlock();
    }
  if ( ( quantity == "HII" ) || ( quantity == "GII" ) )
    {
      trace.beginBlock( "Compute II curvature estimations" );
      measured_values = ( quantity == "HII" )
	? SHG::getIIMeanCurvatures    ( bimage, dft_surfels, params )
	: SHG::getIIGaussianCurvatures( bimage, dft_surfels, params );
      stat_measured_curv.addValues( measured_values.cbegin(), measured_values.cend() );
      stat_measured_curv.terminate();
      trace.info() << "- II curv: avg = " << stat_measured_curv.mean() << std::endl;
      trace.info() << "- II curv: min = " << stat_measured_curv.min() << std::endl;
      trace.info() << "- II curv: max = " << stat_measured_curv.max() << std::endl;
      time_curv_estimations = trace.endBlock();
    }
  else 
    {
      trace.beginBlock( "Compute corrected normal current" );
      Current C( *idx_surface, h, vm.count( "crisp" ) );
      C.setCorrectedNormals( dft_surfels.begin(), dft_surfels.end(), measured_normals.begin() );
      C.setInterpolationMode( vm.count( "interpolate" ) );
      trace.info() << C << " m-ball-r = " << mr << "(continuous)"
		   << " " << (mr/h) << " (discrete)" << std::endl;
      trace.info() << "Interpolation mode is "
		   << ( vm.count( "interpolate" ) ? "ON" : "OFF" ) << std::endl;
      double              area = 0.0;
      double              intG = 0.0;
      std::vector<double> mu0( dft_surfels.size() );
      std::vector<double> mu1( dft_surfels.size() );
      std::vector<double> mu2( dft_surfels.size() );
      std::vector<double> muOmega( dft_surfels.size() );
      SH::RealVectors     muDir0( dft_surfels.size() );
      SH::RealVectors     muDir1( dft_surfels.size() );
      SH::RealVectors     muDir2( dft_surfels.size() );
      SH::Scalars         kappa0( dft_surfels.size() );
      SH::Scalars         kappa1( dft_surfels.size() );
      SH::Scalars         kappa2( dft_surfels.size() );
      bool       mu0_needed = true;
      bool       mu1_needed = false;
      bool       mu2_needed = false;
      bool   muOmega_needed = false;
      bool      muXY_needed = false;
      if ( quantity == "Mu0" )     mu0_needed = true;
      if ( quantity == "Mu1" )     mu1_needed = true;
      if ( quantity == "Mu2" )     mu2_needed = true;
      if ( quantity == "MuOmega" ) muOmega_needed = true;
      if ( quantity == "H" )       mu0_needed = mu1_needed = true;
      if ( quantity == "G" )       mu0_needed = mu2_needed = true;
      if ( quantity == "Omega" )   mu0_needed = muOmega_needed = true;
      if ( quantity == "MuXY" )    muXY_needed = true;
      if ( anisotropy == "NMult" || anisotropy == "NAdd" ) mu0_needed = true;
      trace.beginBlock( "Compute mu_k everywhere" );
      trace.info() << "computeAllMu0" << std::endl;
      if ( mu0_needed ) C.computeAllMu0();
      trace.info() << "computeAllMu1" << std::endl;
      if ( mu1_needed ) C.computeAllMu1();
      trace.info() << "computeAllMu2" << std::endl;
      if ( mu2_needed ) C.computeAllMu2();
      trace.info() << "computeAllMuOmega" << std::endl;
      if ( muOmega_needed ) C.computeAllMuOmega();
      trace.info() << "computeAllMuXYAnisotropic" << std::endl;
      if ( muXY_needed ) C.computeAllAnisotropicMu();
      time_mu_estimations = trace.endBlock();
      //#pragma omp parallel for schedule(dynamic)
      trace.info() << "compute measures" << std::endl;
      Vertex              i = 0;
      Vertex              j = dft_surfels.size();
      const bool normalized = ( anisotropy == "NMult" ) || ( anisotropy == "NAdd" );
      const bool M_times_Mt = ( anisotropy == "Mult" )  || ( anisotropy == "NMult" );
      for ( auto aSurfel : dft_surfels )
	{
	  // std::cout << i << " / " << j << std::endl;
	  trace.progressBar( i, j );
	  Vertex v = C.getVertex( aSurfel );
	  area    += C.mu0( v );
	  if ( mu0_needed )     mu0[ i ]     = C.mu0Ball( v, mr );
	  if ( mu1_needed )     mu1[ i ]     = C.mu1Ball( v, mr );
	  if ( mu2_needed )     mu2[ i ]     = C.mu2Ball( v, mr );
	  if ( muOmega_needed ) muOmega[ i ] = C.muOmegaBall( v, mr );
	  if ( muXY_needed && ! M_times_Mt ) {
	    auto M = C.anisotropicMuBall( v, mr );
	    auto N = measured_normals[ i ];
	    M += M.transpose();
	    M *= 0.5;
	    const double   coef_N = normalized ? 100.0 * mu0[ i ] : 100.0;
	    for ( int j = 0; j < 3; j++ )
	       for ( int k = 0; k < 3; k++ )
		 M( j, k ) += coef_N * N[ j ] * N[ k ];
	    auto V = M;
	    RealVector L;
	    EigenDecomposition< 3, double>::getEigenDecomposition( M, V, L );
	    muDir0[ i ] = V.column( 0 );
	    muDir1[ i ] = V.column( 1 );
	    muDir2[ i ] = V.column( 2 );
            kappa0[ i ] = normalized ? -L[ 0 ] / mu0[ i ] : -L[ 0 ];
            kappa1[ i ] = normalized ? -L[ 1 ] / mu0[ i ] : -L[ 1 ];
            kappa2[ i ] = normalized ? -L[ 2 ] / mu0[ i ] : -L[ 2 ];
	  } else if ( muXY_needed && M_times_Mt ) {
	    auto M = C.anisotropicMuBall( v, mr );
	    auto MMt = M * M.transpose();
	    auto V = M;
	    RealVector L;
	    EigenDecomposition< 3, double>::getEigenDecomposition( MMt, V, L );
	    muDir0[ i ] = V.column( 0 );
	    muDir1[ i ] = V.column( 1 );
	    muDir2[ i ] = V.column( 2 );
            kappa0[ i ] = sqrt( fabs( L[ 0 ] ) ) / ( normalized ? mu0[ i ] : 1.0 );
            kappa1[ i ] = sqrt( fabs( L[ 1 ] ) ) / ( normalized ? mu0[ i ] : 1.0 );
            kappa2[ i ] = sqrt( fabs( L[ 2 ] ) ) / ( normalized ? mu0[ i ] : 1.0 );
	  }
	  ++i;
	}
      // Computing total Gauss curvature.
      if ( mu2_needed )
	{
	  trace.info() << "compute total Gauss curvature" << std::endl;
	  for ( auto f : idx_surface->allFaces() ) intG += C.mu2( f );
	}
      if ( quantity == "Mu0" )          measured_values = mu0;
      else if ( quantity == "Mu1" )     measured_values = mu1;
      else if ( quantity == "Mu2" )     measured_values = mu2;
      else if ( quantity == "MuOmega" ) measured_values = muOmega;
      else if ( quantity == "H" )
	{
	  measured_values.resize( dft_surfels.size() );
	  std::transform( mu0.cbegin(), mu0.cend(),
			  mu1.cbegin(), measured_values.begin(),
			  [] ( double m0, double m1 ) { return m1 / (2.0*m0); } );
	}
      else if ( quantity == "G" )
	{
	  measured_values.resize( dft_surfels.size() );
	  std::transform( mu0.cbegin(), mu0.cend(),
			  mu2.cbegin(), measured_values.begin(),
			  [] ( double m0, double m2 ) { return m2 / m0; } );
	}
      else if ( quantity == "Omega" )
	{
	  measured_values.resize( dft_surfels.size() );
	  std::transform( mu0.cbegin(), mu0.cend(),
			  muOmega.cbegin(), measured_values.begin(),
			  [] ( double m0, double m2 ) { return m2 / sqrt( m0 ); } );
	}
      if ( quantity != "MuXY" )
	{
	  for ( i = 0; i < j; ++i ) stat_measured_curv.addValue( measured_values[ i ] );
	  stat_measured_curv.terminate();
	  trace.info() << "- CNC area      = " << area << std::endl;
	  trace.info() << "- CNC total G   = " << intG << std::endl;
	  trace.info() << "- CNC mu: avg = " << stat_measured_curv.mean() << std::endl;
	  trace.info() << "- CNC mu: min = " << stat_measured_curv.min() << std::endl;
	  trace.info() << "- CNC mu: max = " << stat_measured_curv.max() << std::endl;
	}
      time_curv_estimations = trace.endBlock();
      if ( quantity == "MuXY" ) {
        trace.beginBlock( "Output principal curvatures and directions files." );
	SH::RealPoints positions( dft_surfels.size() );
	std::transform( dft_surfels.cbegin(), dft_surfels.cend(), positions.begin(),
			[&] (const SH::SCell& c) { return sembedder( c ); } );
	auto pos0 = positions,   pos1 = positions,   pos2 = positions;
        auto max0 = kappa0[ 0 ], max1 = kappa1[ 0 ], max2 = kappa2[ 0 ];
        auto min0 = kappa0[ 0 ], min1 = kappa1[ 0 ], min2 = kappa2[ 0 ];
	for ( SH::Idx i = 0; i < pos0.size(); ++i ) {
	  muDir0[ i ] *= vlength;
	  muDir1[ i ] *= vlength;
	  muDir2[ i ] *= vlength;
	  pos0[ i ] -= 0.5 * muDir0[ i ];
	  pos1[ i ] -= 0.5 * muDir1[ i ];
	  pos2[ i ] -= 0.5 * muDir2[ i ];
          max0 = std::max( max0, kappa0[ i ] );
          max1 = std::max( max1, kappa1[ i ] );
          max2 = std::max( max2, kappa2[ i ] );
          min0 = std::min( min0, kappa0[ i ] );
          min1 = std::min( min1, kappa1[ i ] );
          min2 = std::min( min2, kappa2[ i ] );
	}
        trace.info() << min0 << " <= kappa0 <= " << max0 << std::endl;
        trace.info() << min1 << " <= kappa1 <= " << max1 << std::endl;
        trace.info() << min2 << " <= kappa2 <= " << max2 << std::endl;
	SH::saveOBJ( surface, SH::RealVectors(), SH::Colors(),
		     outputfile+"-primal.obj" );
        auto colors = SH::Colors( dft_surfels.size() );
	// TODO 
	const auto cmap_kappa
	  = SH::getZeroTickedColorMap( minValue, maxValue,
				       params( "colormap", colormap_name )
				       ( "zero-tic", zt ) );
        for ( SH::Idx i = 0; i < colors.size(); i++ )
	  colors[ i ] = cmap_kappa( kappa0[ i ] );
	SH::saveVectorFieldOBJ( pos0, muDir0, vthickness, getSimplifiedColorMap( colors ), outputfile+"-dir0.obj" );
        for ( SH::Idx i = 0; i < colors.size(); i++ )
	  colors[ i ] = cmap_kappa( kappa1[ i ] );
	SH::saveVectorFieldOBJ( pos1, muDir1, vthickness, getSimplifiedColorMap( colors ), outputfile+"-dir1.obj" );
	// normal direction, color has no meaning
	SH::saveVectorFieldOBJ( pos2, muDir2, vthickness, SH::Colors(),
		     outputfile+"-dir2.obj", SH::Color::Black, SH::Color::Blue );
	const auto colormap0 = SH::getColorMap( min0, max0,
						params( "colormap", colormap_name ) );
        for ( SH::Idx i = 0; i < colors.size(); i++ )
	  colors[ i ] = colormap0( kappa0[ match[ i ] ] );
        SH::saveOBJ( surface, SH::RealVectors(), getSimplifiedColorMap( colors ), outputfile+"-kappa0.obj" );
        const auto colormap1 = SH::getColorMap( min1, max1,
						params( "colormap", colormap_name ) );
        for ( SH::Idx i = 0; i < colors.size(); i++ )
	  colors[ i ] = colormap1( kappa1[ match[ i ] ] );
        SH::saveOBJ( surface, SH::RealVectors(), getSimplifiedColorMap( colors ), outputfile+"-kappa1.obj" );
        const auto colormap2 = SH::getColorMap( min2, max2,
						params( "colormap", colormap_name ) );
        for ( SH::Idx i = 0; i < colors.size(); i++ )
	  colors[ i ] = colormap2( kappa2[ match[ i ] ] );
        SH::saveOBJ( surface, SH::RealVectors(), getSimplifiedColorMap( colors ), outputfile+"-kappa2.obj" );

	SH::saveOBJ( surface, SH::RealVectors(), SH::Colors(),
		     outputfile+"-primal.obj" );
	if ( vm.count( "polynomial" ) ) {
	  SH::Cell2Index c2i;
	  auto pointels = SH::getPointelRange( c2i, surface );
	  SH::RealPoints pos( pointels.size() );
	  std::transform( pointels.cbegin(), pointels.cend(), pos.begin(),
			  [&] (const SH::Cell& c) { return h * embedder( c ); } ); 
	  auto ppos     = SHG::getPositions( shape, pos, params );
	  std::transform( ppos.cbegin(), ppos.cend(), pos.begin(),
			  [&] (const SH::RealPoint& p) { return p / h; } );
	  SH::saveOBJ( surface, [&] (const SH::Cell& c){ return pos[ c2i[ c ] ];},
		       SH::RealVectors(), SH::Colors(),
		       outputfile+"-qproj.obj" );
	  SH::RealPoints pos_surf( dft_surfels.size() );
	  std::transform( dft_surfels.cbegin(), dft_surfels.cend(), pos_surf.begin(),
			  [&] (const SH::SCell& c) { return h * sembedder( c ); } ); 
	  auto ppos_surf = SHG::getPositions( shape, pos_surf, params );
	  std::transform( ppos_surf.cbegin(), ppos_surf.cend(), pos_surf.begin(),
			  [&] (const SH::RealPoint& p) { return p / h; } );
	  // auto real_pos = SH::getMatchedRange( pos_surf, match );
	  pos0 = pos1 = pos2 = pos_surf; //real_pos;
	  for ( SH::Idx i = 0; i < pos0.size(); ++i ) {
	    pos0[ i ] -= 0.5 * muDir0[ i ];
	    pos1[ i ] -= 0.5 * muDir1[ i ];
	    pos2[ i ] -= 0.5 * muDir2[ i ];
	  }
        for ( SH::Idx i = 0; i < colors.size(); i++ )
	  colors[ i ] = cmap_kappa( kappa0[ i ] );
	SH::saveVectorFieldOBJ( pos0, muDir0, vthickness, colors, outputfile+"-qproj-dir0.obj" );
        for ( SH::Idx i = 0; i < colors.size(); i++ )
	  colors[ i ] = cmap_kappa( kappa1[ i ] );
	SH::saveVectorFieldOBJ( pos1, muDir1, vthickness, colors, outputfile+"-qproj-dir1.obj" );
	SH::saveVectorFieldOBJ( pos2, muDir2, vthickness, SH::Colors(),
		     outputfile+"-qproj-dir2.obj", SH::Color::Black,
		     SH::Color::Blue );
	}
	trace.endBlock();
      }
    }
  
  if ( outputfile != "none" )
    {
      trace.beginBlock( "Save results as OBJ" );
      const bool has_ground_truth = expected_values.size() != 0;
      const bool has_estimations  = measured_values.size() != 0;
      const auto colormap         = SH::getZeroTickedColorMap( minValue, maxValue,
                                                               params( "colormap", colormap_name )
                                                               ( "zero-tic", zt ) );
      
      trace.info() << "#mvalues=" << measured_values.size() << std::endl;
      trace.info() << "#evalues=" << expected_values.size() << std::endl;

      
      auto colors = SH::Colors( surfels.size() );
      if ( has_ground_truth )
	{
	  for ( SH::Idx i = 0; i < colors.size(); i++ )
	    colors[ i ] = colormap( expected_values[ match[ i ] ] ); 
	  SH::saveOBJ( surface, SH::getMatchedRange( expected_normals, match ), colors,
		       outputfile+"-truth.obj" );
	}
      if ( has_estimations )
	{
	  for ( SH::Idx i = 0; i < colors.size(); i++ )
	    colors[ i ] = colormap( measured_values[ match[ i ] ] ); 
	  SH::saveOBJ( surface, SH::getMatchedRange( measured_normals, match ), colors,
		       outputfile+"-estimation.obj" );
	}
      if ( has_ground_truth && has_estimations )
	{
	  const auto error_values = SHG::getScalarsAbsoluteDifference( measured_values, expected_values );
	  const auto max_error    = vm[ "max-error" ].as<double>();
	  const auto stat_error   = SHG::getStatistic( error_values );
	  const auto error_cmap   = getErrorColorMap( max_error );
	  for ( SH::Idx i = 0; i < colors.size(); i++ )
	    colors[ i ] = error_cmap( error_values[ match[ i ] ] ); 
	  SH::saveOBJ( surface, SH::getMatchedRange( measured_normals, match ), colors,
		       outputfile+"-error.obj" );
	}


      trace.endBlock();
    }
  
  
  if ( ! vm.count( "polynomial" ) ) return 0;
  
  trace.beginBlock( "Output statistics" );
  auto error_fname = vm[ "error" ].as<std::string>();
  std::ofstream ferr;
  ferr.open ( error_fname.c_str(),
	      std::ofstream::out | std::ofstream::app );
  if ( ! ferr.good() )
    trace.warning() << "Unable to open file " << error_fname << std::endl;
  ferr << "######################################################################"
       << std::endl;
  ferr << "# ";
  for ( int i = 0; i < argc; ++i ) ferr << " " << argv[ i ];
  ferr << std::endl;
  ferr << "#---------------------------------------------------------------------"
       << std::endl;
  ferr << "# Q=" << quantity
       << " P="  << vm[ "polynomial" ].as<std::string>()
       << " mr= " << mr << "(continuous)"
       << " mrd=" << (mr/h) << " (discrete)" << std::endl;
  ferr << "# time_curv_ground_truth  = " << time_curv_ground_truth << " ms" << std::endl
       << "# time_normal_estimations = " << time_normal_estimations << " ms" << std::endl
       << "# time_curv_estimations   = " << time_curv_estimations << " ms" << std::endl
       << "# time_mu_estimations     = " << time_mu_estimations << " ms" << std::endl;
  ferr << "# h(1) size(2) l1(3) l2(4) loo(5)"
       << " m_mean(6) m_dev(7) m_min(8) m_max(9)"
       << " exp_mean(10) exp_dev(11) exp_min(12) exp_max(13) " << std::endl;
  ferr << "# mr(14) mrd(15) t_curv_gt(16) t_normal_est(17) t_curv_est(18) t_mu_est(19)"
       << " r-radius(20) m-coef(21)"
       << std::endl;
  ferr << h << " " << measured_values.size()
       << " " << SHG::getScalarsNormL1 ( measured_values, expected_values )
       << " " << SHG::getScalarsNormL2 ( measured_values, expected_values )
       << " " << SHG::getScalarsNormLoo( measured_values, expected_values )
       << " " << stat_measured_curv.mean()
       << " " << sqrt( stat_measured_curv.variance() )
       << " " << stat_measured_curv.min()
       << " " << stat_measured_curv.max()
       << " " << stat_expected_curv.mean()
       << " " << sqrt( stat_expected_curv.variance() )
       << " " << stat_expected_curv.min()
       << " " << stat_expected_curv.max();
  ferr << " " << mr << " " << (mr/h) << " " << time_curv_ground_truth
       << " " << time_normal_estimations << " " << time_curv_estimations
       << " " << time_mu_estimations
       << " " << vm[ "r-radius" ].as<double>()
       << " " << vm[ "m-coef" ].as<double>()
       << std::endl;
  ferr << "#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
       << std::endl;
  ferr.close();
  trace.endBlock();
  
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
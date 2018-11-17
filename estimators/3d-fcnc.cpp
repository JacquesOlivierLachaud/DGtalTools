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
    ( "quantity,Q", po::value<std::string>()->default_value( "Mu1" ), "the quantity that is evaluated in Mu0|Mu1|Mu2|MuOmega|H|G|Omega|HII|GII, with H := Mu1/(2Mu0), G := Mu2/Mu0, Omega := MuOmega/sqrt(Mu0), and HII and GII are the mean and gaussian curvatures estimated by II." )
    ( "crisp,C", "when specified, when computing measures in a ball, do not approximate the relative intersection of cells with the ball but only consider if the cell centroid is in the ball (faster by 30%, but less accurate)." );
  EH::optionsDisplayValues   ( general_opt );
  //#endif
  general_opt.add_options()
    ( "error", po::value<std::string>()->default_value( "error.txt" ), "the name of the output file that sum up l2 and loo errors in estimation." );
  general_opt.add_options()
    ( "max-error", po::value<double>()->default_value( 0.2 ), "the error value corresponding to black." );
  general_opt.add_options()
    ( "output,o", po::value<std::string>()->default_value( "fcnc" ), "the basename for output obj files." );
  
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
  auto quantity = vm[ "quantity" ].as<std::string>();
  std::vector< std::string > quantities = { "Mu0", "Mu1", "Mu2", "MuOmega", "H", "G", "Omega", "HII", "GII" };
  if ( std::count( quantities.begin(), quantities.end(), quantity ) == 0 ) {
    trace.error() << "Quantity should be in Mu0|Mu1|Mu2|MuOmega|H|G|Omega|HII|GII.";
    trace.info() << std::endl;
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
      params( "closed",       vm[ "closed"       ].as<int>()    );
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
  params( "surfaceTraversal", "DepthFirst" );
  auto dft_surfels  = SH::getSurfelRange( surface, params );
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
  else if ( estimator == "II" )
    measured_normals = SHG::getIINormalVectors( bimage, dft_surfels, params );
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
      trace.info() << C << " m-ball-r = " << mr << "(continuous)"
		   << " " << (mr/h) << " (discrete)" << std::endl;
      double              area = 0.0;
      double              intG = 0.0;
      std::vector<double> mu0( dft_surfels.size() );
      std::vector<double> mu1( dft_surfels.size() );
      std::vector<double> mu2( dft_surfels.size() );
      std::vector<double> muOmega( dft_surfels.size() );
      bool       mu0_needed = true;
      bool       mu1_needed = false;
      bool       mu2_needed = false;
      bool   muOmega_needed = false;
      if ( quantity == "Mu0" )     mu0_needed = true;
      if ( quantity == "Mu1" )     mu1_needed = true;
      if ( quantity == "Mu2" )     mu2_needed = true;
      if ( quantity == "MuOmega" ) muOmega_needed = true;
      if ( quantity == "H" )       mu0_needed = mu1_needed = true;
      if ( quantity == "G" )       mu0_needed = mu2_needed = true;
      if ( quantity == "Omega" )   mu0_needed = muOmega_needed = true;
      trace.beginBlock( "Compute mu_k everywhere" );
      trace.info() << "computeAllMu0" << std::endl;
      if ( mu0_needed ) C.computeAllMu0();
      trace.info() << "computeAllMu1" << std::endl;
      if ( mu1_needed ) C.computeAllMu1();
      trace.info() << "computeAllMu2" << std::endl;
      if ( mu2_needed ) C.computeAllMu2();
      trace.info() << "computeAllMuOmega" << std::endl;
      if ( muOmega_needed ) C.computeAllMuOmega();
      time_mu_estimations = trace.endBlock();
      //#pragma omp parallel for schedule(dynamic)
      trace.info() << "compute measures" << std::endl;
      Vertex              i = 0;
      Vertex              j = dft_surfels.size();
      for ( auto aSurfel : dft_surfels )
	{
	  // std::cout << i << " / " << j << std::endl;
	  trace.progressBar( i, j );
	  Vertex v = C.getVertex( aSurfel );
	  area    += C.mu0( v );
	  if ( mu0_needed ) mu0[ i ] = C.mu0Ball( v, mr );
	  if ( mu1_needed ) mu1[ i ] = C.mu1Ball( v, mr );
	  if ( mu2_needed ) mu2[ i ] = C.mu2Ball( v, mr );
	  if ( muOmega_needed ) muOmega[ i ] = C.muOmegaBall( v, mr );
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
      for ( i = 0; i < j; ++i ) stat_measured_curv.addValue( measured_values[ i ] );
      stat_measured_curv.terminate();
      trace.info() << "- CNC area      = " << area << std::endl;
      trace.info() << "- CNC total G   = " << intG << std::endl;
      trace.info() << "- CNC curv: avg = " << stat_measured_curv.mean() << std::endl;
      trace.info() << "- CNC curv: min = " << stat_measured_curv.min() << std::endl;
      trace.info() << "- CNC curv: max = " << stat_measured_curv.max() << std::endl;
      time_curv_estimations = trace.endBlock();
    }

  trace.beginBlock( "Save results as OBJ" );
  const bool has_ground_truth = expected_values.size() != 0;
  const bool has_estimations  = measured_values.size() != 0;
  const auto minValue         = vm[ "minValue" ].as<double>();
  const auto maxValue         = vm[ "maxValue" ].as<double>();
  const auto colormap_name    = vm[ "colormap" ].as<string>();
  const auto outputfile       = vm[ "output"   ].as<string>();
  const auto colormap         = SH::getColorMap( minValue, maxValue, colormap_name );

  trace.info() << "#mvalues=" << measured_values.size() << std::endl;
  trace.info() << "#evalues=" << expected_values.size() << std::endl;

  trace.beginBlock( "Compute surfels" );
  params( "surfaceTraversal", "Default" );
  const auto surfels = SH::getSurfelRange( surface, params );
  const auto match   = SH::getRangeMatch ( surfels, dft_surfels );
  trace.endBlock();

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
      const auto error_cmap   = SH::getColorMap( 0.0, max_error,
						 params( "colormap", "Custom" ) );
      for ( SH::Idx i = 0; i < colors.size(); i++ )
	colors[ i ] = error_cmap( error_values[ match[ i ] ] ); 
      SH::saveOBJ( surface, SH::getMatchedRange( measured_normals, match ), colors,
		   outputfile+"-error.obj" );
    }
  trace.endBlock();
  
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
       << std::endl;
  ferr << "#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
       << std::endl;
  ferr.close();
  trace.endBlock();
  
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

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
  //#ifdef WITH_VISU3D_QGLVIEWER
  EH::optionsDisplayValues   ( general_opt );
  general_opt.add_options()
    ( "view,V", po::value<std::string>()->default_value( "Measure" ), "the display mode in Measure|Truth|Error|None" );
  //#endif
  general_opt.add_options()
    ( "error", po::value<std::string>()->default_value( "error.txt" ), "the name of the output file that sum up l2 and loo errors in estimation." );
  
  po::variables_map vm;
  bool parseOK = EH::args2vm( general_opt, argc, argv, vm );
  bool neededArgsGiven=true;

  if ( vm.count( "polynomial-list" ) )
    {
      trace.info() << "List of predefined polynomials:" << std::endl;
      auto L = SG::getPolynomialList();
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
  KSpace                        K;
  unsigned int                 nb = 0;
  CountedPtr<BinaryImage>  bimage( nullptr );
  CountedPtr<ImplicitShape> shape( nullptr );
  auto params = SH::defaultParameters() | SHG::defaultParameters();
  trace.beginBlock( "Make Shape" );
  // Generic parameters.
  params( "gristep",           vm[ "gridstep" ].as<double>() );
  params( "noise",             vm[ "noise"    ].as<double>() );
  params( "surfelAdjacency",          0 ); // 0:interior
  params( "nbTriesToFindABel",   100000 ); // number of tries in method Surfaces::findABel
  params( "surfaceComponents", "AnyBig" ); // "AnyBig"|"All"
  params( "projectionMaxIter",       20 ); // the maximum number of iter for the projection.
  params( "projectionAccuracy",  0.0001 ); // zero-proximity stop crit. during projection.
  params( "projectionGamma",        0.5 ); // the displacement coef. of the projection.
  params( "verbose",                  1 );
  params( "t-ring",            vm[ "trivial-ring" ].as<double>() );
  params( "kernel",            vm[ "kernel"       ].as<double>() );
  params( "R-radius",          vm[ "R-radius"     ].as<double>() );
  params( "r-radius",          vm[ "r-radius"     ].as<double>() );
  params( "alpha",             vm[ "alpha"        ].as<double>() );
  params( "surfelEmbedding",   vm[ "embedding"    ].as<double>() );
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
      bimage       = EH::makeBinaryImage( dshape, params );
    }
  else if ( vm.count( "input" ) )
    {
      // Fill useful parameters
      params( "thresholdMin", vm[ "thresholdMin" ].as<string>() );
      params( "thresholdMax", vm[ "thresholdMax" ].as<string>() );
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
  
  auto surface = EH::makeIdxDigitalSurface( bimage, params );
  if ( surface == 0 ) {
    trace.error() << "- surface is empty (either empty or full volume). ";
    trace.info()  << std::endl;
    trace.endBlock();
    return 1;
  }
  trace.info() << "- surface has " << surface->size()<< " surfels." << std::endl;
  trace.endBlock();

  auto estimator     = vm[ "estimator" ].as<string>();
  auto view          = vm[ "view"      ].as<string>();
  const double h     = vm[ "gridstep"  ].as<double>();
  const double mcoef = vm[ "m-coef"    ].as<double>();
  const double mpow  = vm[ "m-pow"     ].as<double>();
  const double mr    = mcoef * pow( h, mpow );
  std::vector<double>     displayed_values;
  std::vector<double>     measured_values;
  std::vector<double>     expected_values;
  std::vector<RealVector> normals;
  Statistic<double>       stat_exp_curv;
  Statistic<double>       stat_meas_curv;
  
  trace.beginBlock( "Compute surfels" );
  params( "surfaceTraversal", "DepthFirst" );
  auto dft_surfels = SH::getSurfelRange( surface, params );
  trace.endBlock();

  if ( vm.count( "polynomial" ) ) {
    trace.beginBlock( "Compute true normals" );
    normals         = SHG::getNormalVectors( shape, K, surfels, params );
    trace.endBlock();
    trace.beginBlock( "Compute true curvatures" );
    expected_values = ( ( quantity == "H" ) || ( quantity == "Mu1" )
			|| ( quantity == "HII" ) )
      ? SHG::getMeanCurvatures    ( shape, K, surfels, params )
      : SHG::getGaussianCurvatures( shape, K, surfels, params );
    stat_exp_curv.addValues( expected_values.cbegin(), expected_values.cend() );
    stat_exp_curv.terminate();
    trace.info() << "- truth curv: avg = " << stat_exp_curv.mean() << std::endl;
    trace.info() << "- truth curv: min = " << stat_exp_curv.min() << std::endl;
    trace.info() << "- truth curv: max = " << stat_exp_curv.max() << std::endl;
    trace.endBlock();
  } else {
    trace.beginBlock( "Estimate normals" );
      if ( estimator == "Trivial" )
	normals = SHG::getTrivialNormalVectors( K, surfels );
      else if ( estimator == "CTrivial" )
	normals = SHG::getCTrivialNormalVectors( surface, surfels, params );
      else if ( estimator == "VCM" )
	normals = SHG::getVCMNormalVectors( surface, surfels, params );
      else if ( estimator == "II" )
	normals = SHG::getIINormalVectors( bimage, surfels, params );
    trace.endBlock();
  }
  if ( ( quantity == "HII" ) || ( quantity == "GII" ) )
    {
      trace.beginBlock( "Compute II curvature estimations" );
      measured_values = ( quantity == "HII" )
	? SHG::getIIMeanCurvatures    ( bimage, surfels, params )
	: SHG::getIIGaussianCurvatures( bimage, surfels, params );
      stat_meas_curv.addValues( measured_values.cbegin(), measured_values.cend() );
      stat_meas_curv.terminate();
      trace.info() << "- II curv: avg = " << stat_meas_curv.mean() << std::endl;
      trace.info() << "- II curv: min = " << stat_meas_curv.min() << std::endl;
      trace.info() << "- II curv: max = " << stat_meas_curv.max() << std::endl;
      trace.endBlock();
    }
  else // if ( view == "Measure" || view == "Error" )
    {
      trace.beginBlock( "Computing corrected normal current" );
      Current C( *surface, h, vm.count( "crisp" ) );
      C.setCorrectedNormals( surfels.begin(), surfels.end(), normals.begin() );
      trace.info() << C << " m-ball-r = " << mr << "(continuous)"
		   << " " << (mr/h) << " (discrete)" << std::endl;
      double              area = 0.0;
      double              intG = 0.0;
      std::vector<double> mu0( surfels.size() );
      std::vector<double> mu1( surfels.size() );
      std::vector<double> mu2( surfels.size() );
      std::vector<double> muOmega( surfels.size() );
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
      trace.info() << "computeAllMu0" << std::endl;
      if ( mu0_needed ) C.computeAllMu0();
      trace.info() << "computeAllMu1" << std::endl;
      if ( mu1_needed ) C.computeAllMu1();
      trace.info() << "computeAllMu2" << std::endl;
      if ( mu2_needed ) C.computeAllMu2();
      trace.info() << "computeAllMuOmega" << std::endl;
      if ( muOmega_needed ) C.computeAllMuOmega();
      //#pragma omp parallel for schedule(dynamic)
      trace.info() << "compute measures" << std::endl;
      Vertex              i = 0;
      Vertex              j = surfels.size();
      for ( auto aSurfel : surfels )
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
	  for ( auto f : surface->allFaces() ) intG += C.mu2( f );
	}
      if ( quantity == "Mu0" )          measured_values = mu0;
      else if ( quantity == "Mu1" )     measured_values = mu1;
      else if ( quantity == "Mu2" )     measured_values = mu2;
      else if ( quantity == "MuOmega" ) measured_values = muOmega;
      else if ( quantity == "H" )
	{
	  measured_values.resize( surfels.size() );
	  std::transform( mu0.cbegin(), mu0.cend(),
			  mu1.cbegin(), measured_values.begin(),
			  [] ( double m0, double m1 ) { return m1 / (2.0*m0); } );
	}
      else if ( quantity == "G" )
	{
	  measured_values.resize( surfels.size() );
	  std::transform( mu0.cbegin(), mu0.cend(),
			  mu2.cbegin(), measured_values.begin(),
			  [] ( double m0, double m2 ) { return m2 / m0; } );
	}
      else if ( quantity == "Omega" )
	{
	  measured_values.resize( surfels.size() );
	  std::transform( mu0.cbegin(), mu0.cend(),
			  muOmega.cbegin(), measured_values.begin(),
			  [] ( double m0, double m2 ) { return m2 / sqrt( m0 ); } );
	}
      for ( i = 0; i < j; ++i ) stat_meas_curv.addValue( measured_values[ i ] );
      stat_meas_curv.terminate();
      trace.info() << "- CNC area      = " << area << std::endl;
      trace.info() << "- CNC total G   = " << intG << std::endl;
      trace.info() << "- CNC curv: avg = " << stat_meas_curv.mean() << std::endl;
      trace.info() << "- CNC curv: min = " << stat_meas_curv.min() << std::endl;
      trace.info() << "- CNC curv: max = " << stat_meas_curv.max() << std::endl;
      trace.endBlock();
    }

#ifdef WITH_VISU3D_QGLVIEWER
  if ( view != "None" ) {
    typedef Viewer3D<Space,KSpace> MyViewever3D;
    typedef Display3DFactory<Space,KSpace> MyDisplay3DFactory;
    
    trace.beginBlock( "View measure" );
    trace.info() << "view mode is " << view << std::endl;
    trace.info() << "#mvalues=" << measured_values.size() << std::endl;
    trace.info() << "#evalues=" << expected_values.size() << std::endl;
    MyViewever3D viewer( K );
    viewer.show();
    if ( view == "Measure" ) 
      displayed_values = measured_values;
    else if ( view == "Truth" )
      displayed_values = expected_values;
    else if ( view == "Error" )
      displayed_values = EH::absoluteDifference( measured_values, expected_values );
    trace.info() << "#surfels=" << surfels.size() << std::endl;
    trace.info() << "#dvalues=" << displayed_values.size() << std::endl;
    EH::viewSurfelValues( viewer, vm, surfels, displayed_values, normals );
    EH::viewSurfaceIsolines( viewer, vm, surface, surfels, displayed_values );
    viewer << MyViewever3D::updateDisplay;
    application.exec();
    trace.endBlock();
  }
#endif
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
  ferr << "# h size l1 l2 loo"
       << " m_mean m_dev m_min m_max"
       << " exp_mean exp_dev exp_min exp_max " << std::endl;
  ferr << h << " " << measured_values.size()
       << " " << EH::normL1 ( measured_values, expected_values )
       << " " << EH::normL2 ( measured_values, expected_values )
       << " " << EH::normLoo( measured_values, expected_values )
       << " " << stat_meas_curv.mean()
       << " " << sqrt( stat_meas_curv.variance() )
       << " " << stat_meas_curv.min()
       << " " << stat_meas_curv.max()
       << " " << stat_exp_curv.mean()
       << " " << sqrt( stat_exp_curv.variance() )
       << " " << stat_exp_curv.min()
       << " " << stat_exp_curv.max()
       << std::endl;
  ferr << "#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
       << std::endl;
  ferr.close();
  trace.endBlock();
  
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

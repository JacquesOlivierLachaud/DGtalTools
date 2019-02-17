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
 * @file measure-normals.cpp
 * @ingroup surfaceTools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2019/02/17
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


using namespace DGtal;
using namespace functors;
using namespace std;

/**
 @page DocMeasureNormals measure-normals

 @brief Computes statistics about digital normal estimators (accuracy and Hölder/Lipschitz property).

 @b Usage:  measure-normals -p polynomial -e II

 @b Allowed @b options @b are:

 @code
  -h [ --help ]                         display this message
  -p [ --polynomial ] arg               polynomial
 @endcode

 @b Example:
 @see
 @ref measure-normals.cpp
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
  typedef SH::Surfel                            Surfel;
  
  // parse command line ----------------------------------------------
  po::options_description general_opt( "Allowed options are" );
  general_opt.add_options()  ( "help,h", "display this message" );
  EH::optionsImplicitShape   ( general_opt );
  EH::optionsDigitizedShape  ( general_opt );
  EH::optionsNormalEstimators( general_opt );
  general_opt.add_options()
    ( "error", po::value<std::string>()->default_value( "error.txt" ), "the name of the output file that sum up l2 and loo errors in estimation." );
  po::variables_map vm;
  bool parseOK         = EH::args2vm( general_opt, argc, argv, vm );
  bool neededArgsGiven = true;
  if ( vm.count( "polynomial-list" ) )
    {
      trace.info() << "List of predefined polynomials:" << std::endl;
      auto L = SH::getPolynomialList();
      for ( auto p : L ) {
	trace.info() << "  " << p.first << " -> " << p.second << std::endl;
      }
    }

  if (parseOK && ( ! vm.count("polynomial") ) ) {
    missingParam("--polynomial");
    neededArgsGiven = false;
  }
  
  if ( !neededArgsGiven || !parseOK || vm.count("help") || argc <= 1 )
    {
      trace.info()<< "Computes statistics about digital normal estimators (accuracy and Hölder/Lipschitz property)." <<std::endl
                  << general_opt << "\n"
                  << "Basic usage: "<<std::endl
		  << "\t measure-normals -e II -t 3 -r 3 -p goursat -g 0.5" << std::endl
                  << std::endl;
      return 0;
    }

  // Building space, implicit shape, digitized shape.
  KSpace                      K;
  unsigned int                nb = 0;
  CountedPtr<BinaryImage>     bimage( nullptr );
  CountedPtr<ImplicitShape3D> shape ( nullptr );
  auto params = SH::defaultParameters() | SHG::defaultParameters();
  trace.beginBlock( "Make Shape" );
  // Generic parameters.
  params( "gridstep",          vm[ "gridstep" ].as<double>() );
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
  params( "polynomial",        vm[ "polynomial" ].as<string>() );
  params( "minAABB",           vm[ "minAABB"    ].as<double>() );
  params( "maxAABB",           vm[ "maxAABB"    ].as<double>() );
  params( "offset",            5.0 );
  shape        = SH::makeImplicitShape3D( params );
  K            = SH::getKSpace( params );
  auto dshape  = SH::makeDigitizedImplicitShape3D( shape, params );
  bimage       = SH::makeBinaryImage( dshape, params );
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
  if ( surface == 0 ) {
    trace.error() << "- surface is empty (either empty or full volume). ";
    trace.info()  << std::endl;
    trace.endBlock();
    return 1;
  }
  trace.info() << "- surface has " << surface->size()<< " surfels." << std::endl;
  trace.endBlock();

  // Compute surfels.
  trace.beginBlock( "Compute surfels" );
  params( "surfaceTraversal", "DepthFirst" );
  const auto dft_surfels = SH::getSurfelRange( surface, params );
  trace.endBlock();

  
  // Compute true and estimated normal vectors
  SH::RealVectors expected_normals;
  SH::RealVectors measured_normals;
  auto            estimator = vm[ "estimator" ].as<string>();
  trace.beginBlock( "Compute true and estimated normals" );
  expected_normals = SHG::getNormalVectors( shape, K, dft_surfels, params );
  if ( estimator == "True" )
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
  auto time_normal_estimations = trace.endBlock();

  trace.beginBlock( "Compute statistics" );
  auto normal_angle_dev = SHG::getVectorsAngleDeviation( expected_normals,
							 measured_normals );
  auto zeroes           = SHG::Scalars( normal_angle_dev.size(), 0.0 );
  auto stat_angle_error = SHG::getStatistic( normal_angle_dev );
  trace.info() << "n-u "
	       << " l1=" << SHG::getScalarsNormL1 ( normal_angle_dev, zeroes )
	       << " l2=" << SHG::getScalarsNormL2 ( normal_angle_dev, zeroes )
	       << " loo=" << SHG::getScalarsNormLoo ( normal_angle_dev, zeroes );
  Statistic<double> delta_u;
  std::map<Surfel,RealVector> est_normal;
  for ( unsigned int i = 0; i < dft_surfels.size(); ++i )
    est_normal[ dft_surfels[ i ] ] = measured_normals[ i ];
  for ( auto v : *surface ) {
    std::vector<SH::Vertex> N;
    auto est_normal_v = est_normal[ v ];
    auto out_it       = std::back_inserter( N );
    surface->writeNeighbors( out_it, v );
    for ( auto sn : N )
      delta_u.addValue( ( est_normal_v - est_normal[ sn ] ).norm() );
  }
  delta_u.terminate();
  trace.info() << "delta_u "
	       << " mean=" << delta_u.mean()
	       << " loo=" << delta_u.max();
  trace.endBlock();
  return 0;
}

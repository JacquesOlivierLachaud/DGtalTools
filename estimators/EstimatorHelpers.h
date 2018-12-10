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

#pragma once

/**
 * @file EstimatorHelpers.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 * @author Marion Foare (\c marion.foare@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2017/06/17
 *
 * Header file for module EstimatorHelpers.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(EstimatorHelpers_RECURSES)
#error Recursive header files inclusion detected in EstimatorHelpers.h
#else // defined(EstimatorHelpers_RECURSES)
/** Prevents recursive inclusion of headers. */
#define EstimatorHelpers_RECURSES

#if !defined EstimatorHelpers_h
/** Prevents repeated inclusion of headers. */
#define EstimatorHelpers_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <sstream>
#include <string>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
// always include EigenSupport.h before any other Eigen headers
// #include "DGtal/math/linalg/EigenSupport.h"
#include "DGtal/base/Common.h"
#include "DGtal/base/Clone.h"
#include "DGtal/math/MPolynomial.h"
#include "DGtal/math/Statistic.h"
#include "DGtal/io/readers/MPolynomialReader.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/shapes/ShapeGeometricFunctors.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/SurfelAdjacency.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/geometry/volumes/KanungoNoise.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/geometry/surfaces/estimation/TrueDigitalSurfaceLocalEstimator.h"
#include "DGtal/geometry/surfaces/estimation/VoronoiCovarianceMeasureOnDigitalSurface.h"
#include "DGtal/geometry/surfaces/estimation/VCMDigitalSurfaceLocalEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IIGeometricFunctors.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantVolumeEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantCovarianceEstimator.h"
#include "DGtal/io/readers/GenericReader.h"
#ifdef WITH_VISU3D_QGLVIEWER
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/Display3DFactory.h"
#endif


//////////////////////////////////////////////////////////////////////////////
/* namespace DGtal { */
/*   namespace functors { */
/*     namespace ShapeGeometricFunctors { */
/*     /\** */
/*      * Description of template class 'ShapeGaussianCurvatureFunctor' <p> */
/*      * \brief Aim: A functor RealPoint -> Quantity that returns the */
/*      * gaussian curvature at given point. */
/*      * */
/*      * @tparam TShape the type of the shape where geometric estimation */
/*      * are made. It must have method \a gaussianCurvature. */
/*      *\/ */
/*     template <typename TShape> */
/*     struct ShapeGaussianCurvatureFunctor { */
/*       typedef TShape Shape; */
/*       typedef typename Shape::RealPoint RealPoint; */
/*       typedef typename Shape::RealVector RealVector; */
/*       typedef typename RealVector::Component Scalar; */
/*       typedef RealPoint Argument; */
/*       typedef Scalar Quantity; */
/*       typedef Quantity Value; */
      
/*       /\** */
/*        * Constructor. A shape may also be attached at construction. */
/*        * */
/*        * @param aShape the shape of interest. The alias can be secured */
/*        * if a some counted pointer is handed. */
/*        *\/ */
/*       ShapeGaussianCurvatureFunctor( ConstAlias<Shape> aShape = 0 ) : myShape( aShape ) {} */
      
/*       /\** */
/*        * Attach a shape. */
/*        * */
/*        * @param aShape the shape of interest. The alias can be secured */
/*        * if a some counted pointer is handed. */
/*        *\/ */
/*       void attach( ConstAlias<Shape> aShape ) */
/*       { */
/*         myShape = aShape; */
/*       } */

/*       /\** */
/*          Map operator RealPoint -> RealVector giving the normal vector. */
/*          @param p any point on the shape. */
/*          @return the normal at point p (as the normalized gradient). */
/*       *\/ */
/*       Quantity operator()( const RealPoint & p ) const */
/*       { */
/*         return myShape->gaussianCurvature( p ); */
/*       } */

/*     private: */
/*       /// The shape of interest. */
/*       CountedConstPtrOrConstPtr<Shape> myShape; */
/*     }; */
/*     } */
/*   } */
/* } */
namespace DGtal
{
  namespace po  = boost::program_options;
  namespace sgf = functors::ShapeGeometricFunctors;
  /////////////////////////////////////////////////////////////////////////////
  // template class EstimatorHelpers
  template < typename TKSpace >
  struct EstimatorHelpers
  {
    typedef TKSpace                                  KSpace;
    typedef typename KSpace::Space                   Space;
    typedef typename Space::Integer                  Integer;
    typedef typename Space::Point                    Point;
    typedef typename Space::Vector                   Vector;
    typedef typename Space::RealVector               RealVector;
    typedef typename Space::RealPoint                RealPoint;
    typedef typename RealVector::Component           Scalar;
    typedef MPolynomial< Space::dimension, Scalar >  PolynomialN;
    typedef ImplicitPolynomial3Shape<Space>          ImplicitShape;
    typedef GaussDigitizer< Space, ImplicitShape >   ImplicitDigitalShape;
    typedef typename ImplicitDigitalShape::Domain    Domain;
    typedef ImageContainerBySTLVector<Domain, bool>  BinaryImage;
    typedef ImageContainerBySTLVector<Domain, unsigned char> GrayScaleImage;
    typedef LightImplicitDigitalSurface< KSpace, BinaryImage > SurfaceContainer;
    typedef DigitalSurface< SurfaceContainer >       Surface;
    typedef typename Surface::Surfel                 Surfel;
    typedef typename Surface::Cell                   Cell;
    typedef typename Surface::SCell                  SCell;
    typedef typename Surface::Vertex                 Vertex;
    typedef typename Surface::Arc                    Arc;
    typedef typename Surface::ArcRange               ArcRange;
    typedef DGtal::Statistic<Scalar>                 AngleDevStatistic;
    typedef sgf::ShapeNormalVectorFunctor<ImplicitShape>      NormalFunctor;
    typedef sgf::ShapeMeanCurvatureFunctor<ImplicitShape>     MeanCurvatureFunctor;
    typedef sgf::ShapeGaussianCurvatureFunctor<ImplicitShape> GaussianCurvatureFunctor;
    typedef TrueDigitalSurfaceLocalEstimator
    < KSpace, ImplicitShape, NormalFunctor >             TrueNormalEstimator;
    typedef TrueDigitalSurfaceLocalEstimator
    < KSpace, ImplicitShape, MeanCurvatureFunctor >      TrueMeanCurvatureEstimator;
    typedef TrueDigitalSurfaceLocalEstimator
    < KSpace, ImplicitShape, GaussianCurvatureFunctor >  TrueGaussianCurvatureEstimator;
#ifdef WITH_VISU3D_QGLVIEWER
    typedef GradientColorMap<Scalar>                 ColorMap;
    typedef Viewer3D<Space,KSpace>                   Viewer;
    typedef Display3DFactory<Space,KSpace>           DisplayFactory;
#endif    
    // ------------------- parsing related functions -----------------------------
    
    /// Parses command line given as \a argc, \a argv according to \a
    /// opts and fills in storage map \a vm.
    ///
    /// @param[in] opts the allowed options.
    /// @param[in] argc the size of array \a argv.
    ///
    /// @param[in] argv an array of C-style strings containing the
    /// desired command-line options.
    ///
    /// @param[inout] vm the map options to variables that is updated.
    static bool args2vm( const po::options_description& opts,
			 int argc, char** argv,
			 po::variables_map& vm )
    {
      bool parseOK = true;
      try
	{
	  po::store( po::parse_command_line( argc, argv, opts ), vm );
	}
      catch( const std::exception & ex )
	{
	  parseOK = false;
	  trace.error() << "[EstimatorHelpers::args2vm]"
			<< " Error checking program options: " << ex.what() << std::endl;
	}
      return parseOK;
    }
    
    /// Parses command line given as a vector of strings \a args
    /// according to \a opts and fills in storage map \a vm.
    ///
    /// @param[in] opts the allowed options.
    /// @param[in] args a vector of strings containing the
    /// desired command-line options cut by words.
    ///
    /// @param[inout] vm the map options to variables that is updated.
    static bool vectorString2vm( const po::options_description& opts,
				 std::vector< std::string > args,
				 po::variables_map& vm )
    {
      bool parseOK = true;
      try
	{
	  po::store( po::command_line_parser( args ).options( opts ).run(), vm );
	}
      catch( const std::exception & ex )
	{
	  parseOK = false;
	  trace.error() << "[EstimatorHelpers::vectorString2vm]"
			<< " Error checking program options: " << ex.what() << std::endl;
	}
      return parseOK;
    }
    
    /// Parses a string, cuts it into words, analyzes it
    /// according to \a opts and fills in storage map \a vm.
    ///
    /// @param[in] opts the allowed options.
    /// @param[in] s a string containing the
    /// desired command-line options separated by spaces.
    ///
    /// @param[inout] vm the map options to variables that is updated.
    static bool string2vm( const po::options_description& opts,
			   const std::string& s,
			   po::variables_map& vm )
    {
      // construct a stream from the string
      std::stringstream sstr( s );
      // use stream iterators to copy the stream to the vector as
      // whitespace separated strings
      std::istream_iterator<std::string> it(sstr);
      std::istream_iterator<std::string> end;
      std::vector<std::string> results(it, end);
      return vectorString2vm( opts, results, vm );
    }

    
    // ------------------- setting options related functions ---------------------------    
    /// Add options for implicit shape given by a polynomial.
    static void optionsImplicitShape( po::options_description& desc )
    {
      desc.add_options()
	( "polynomial,p", po::value<std::string>(), "the implicit polynomial whose zero-level defines the shape of interest." )
	( "polynomial-list", "displays the list of predefined polynomials (names that you can use with '-p' options).");
    }

    /// Add options for implicit shape digitization.
    static void optionsDigitizedShape( po::options_description& desc )
    {
      desc.add_options()
	("minAABB,a",  po::value<Scalar>()->default_value( -10.0 ), "the min value of the AABB bounding box (domain)" )
	("maxAABB,A",  po::value<Scalar>()->default_value( 10.0 ), "the max value of the AABB bounding box (domain)" )
	("gridstep,g", po::value< Scalar >()->default_value( 1.0 ), "the gridstep that defines the digitization (often called h). " );
    }

    /// Add options for vol file.
    static void optionsVolFile( po::options_description& desc )
    {
      desc.add_options()
        ("input,i", po::value<std::string>(), "vol file (.vol) , pgm3d (.p3d or .pgm3d, pgm (with 3 dims)) file or sdp (sequence of discrete points)" )
	("thresholdMin,m",  po::value<int>()->default_value(0), "threshold min (excluded) to define binary shape" )
	("thresholdMax,M",  po::value<int>()->default_value(255), "threshold max (included) to define binary shape" )
	("closed",  po::value<int>()->default_value(1), "tells if the cellular space is closed (1:default) or open." );
    }

    /// Add options for implicit shape digitization.
    static void optionsNoisyImage( po::options_description& desc )
    {
      desc.add_options()
	("noise,N", po::value<Scalar>()->default_value( 0.0 ), "the Kanungo noise level l=arg, with l^d the probability that a point at distance d is flipped inside/outside." );
    }

    /// Add options for implicit shape digitization.
    static void optionsNormalEstimators( po::options_description& desc )
    {
      desc.add_options()
	("estimator,e", po::value<std::string>()->default_value( "True" ), "the chosen normal estimator: True | VCM | II | Trivial | CTrivial" )
	("R-radius,R", po::value<Scalar>()->default_value( 5 ), "the constant for parameter R in R(h)=R h^alpha (VCM)." )
	("r-radius,r", po::value<Scalar>()->default_value( 3 ), "the constant for parameter r in r(h)=r h^alpha (VCM,II,Trivial)." )
	("kernel,k", po::value<std::string>()->default_value( "hat" ), "the function chi_r, either hat or ball." )
	("alpha", po::value<Scalar>()->default_value( 0.5 ), "the parameter alpha in r(h)=r h^alpha (VCM)." )
	("trivial-ring,t", po::value<Scalar>()->default_value( 3 ), "the parameter t defining the radius for the convolved Trivial estimator. It corresponds to a discrete distance on the surface graph and is not related to the digitization gridstep.Used for reorienting normal directions." )
	("embedding,E", po::value<int>()->default_value( 0 ), "the surfel -> point embedding for VCM estimator: 0: Pointels, 1: InnerSpel, 2: OuterSpel." );
    }

    /// Add options for implicit shape digitization.
    static void optionsDisplayValues( po::options_description& desc )
    {
      desc.add_options()
	("colormap", po::value<std::string>()->default_value( "Tics" ), "the chosen colormap for displaying values." )
	("zero", po::value<double>()->default_value( 0.0 ), "the value of reference, displayed in black." )
	("tics", po::value<double>()->default_value( 1.0 ), "the spacing between values with a basis at the reference value, displayed in grey." )
	("minValue", po::value<double>()->default_value( -0.3 ), "a specified min value associated with lowest color in colormap ." )
	("maxValue", po::value<double>()->default_value( 0.3 ), "a specified max value associated with highest color in colormap ." );
    }
    
    // ------------------- Shapes related functions -----------------------------
    
    /// Returns a map associating a name and a polynomial,
    /// e.g. "sphere1", "x^2+y^2+z^2-1".
    ///
    /// @return the map associating a polynomial to a name.
    static std::map< std::string, std::string >
    getPolynomialList()
    {
      std::vector< std::pair< std::string, std::string > >
	Ps = { { "sphere1", "x^2+y^2+z^2-1" },
	       { "sphere9", "x^2+y^2+z^2-81" },
	       { "ellipsoid", "3*x^2+2*y^2+z^2-90" },
	       { "cylinder", "x^2+2*z^2-90" },
	       { "torus",   "(x^2+y^2+z^2+6*6-2*2)^2-4*6*6*(x^2+y^2)" },
	       { "rcube",   "x^4+y^4+z^4-6561" },
	       { "goursat", "-1*(8-0.03*x^4-0.03*y^4-0.03*z^4+2*x^2+2*y^2+2*z^2)" },
	       { "distel",  "10000-(x^2+y^2+z^2+1000*(x^2+y^2)*(x^2+z^2)*(y^2+z^2))"},
	       { "leopold", "(x^2*y^2*z^2+4*x^2+4*y^2+3*z^2)-100" },
	       { "diabolo", "x^2-(y^2+z^2)^2" },
	       { "heart",   "-1*(x^2+2.25*y^2+z^2-1)^3+x^2*z^3+0.1125*y^2*z^3" },
	       { "crixxi",  "-0.9*(y^2+z^2-1)^2-(x^2+y^2-1)^3" } };
      std::map< std::string, std::string > L;
      for ( auto p : Ps )
	L[ p.first ] = p.second;
      return L;
    }

    /// Builds a 3D implicit shape from argument "-polynomial".
    ///
    /// @param[in] vm the options sets in the variable map (arguments
    /// given to the program). Recognized parameters are given in \ref
    /// optionsImplicitShape.
    ///
    /// @return a smart pointer on the created implicit shape.
    static CountedPtr<ImplicitShape>
    makeImplicitShape( const po::variables_map& vm )
    {
      typedef MPolynomialReader< Space::dimension, Scalar> Polynomial3Reader;
      std::string poly_str = vm[ "polynomial" ].as<std::string>();
      // Recognizes some strings:
      auto PL = getPolynomialList();
      if ( PL[ poly_str ] != "" ) poly_str = PL[ poly_str ];
      PolynomialN poly;
      Polynomial3Reader reader;
      std::string::const_iterator iter
	= reader.read( poly, poly_str.begin(), poly_str.end() );
      if ( iter != poly_str.end() )
	{
	  trace.error() << "[EstimatorHelpers::makeImplicitShape]"
			<< " ERROR reading polynomial: I read only <"
			<< poly_str.substr( 0, iter - poly_str.begin() )
			<< ">, and I built P=" << poly << std::endl;
	}
      return CountedPtr<ImplicitShape>( new ImplicitShape( poly ) );
    }

    /// Makes the Gauss digitization of the given implicit shape
    /// according to parameters given as arguments to the program.
    ///
    /// @param[in] vm the options sets in the variable map (arguments
    /// given to the program). Recognized parameters are given in \ref
    /// optionsDigitizedShape.
    ///
    /// @param[in] shape the implicit shape.
    /// @param[out] K the Khalimsky space whose domain encompasses the digital shape.
    /// @return a smart pointer on the created implicit digital shape.
    static CountedPtr<ImplicitDigitalShape>
    makeImplicitDigitalShapeFromImplicitShape
    ( const po::variables_map& vm,
      CountedPtr<ImplicitShape> shape,
      KSpace& K )
    {
      Scalar min_x = vm[ "minAABB" ].as<Scalar>();
      Scalar max_x = vm[ "maxAABB" ].as<Scalar>();
      Scalar h     = vm[ "gridstep" ].as<Scalar>();
      RealPoint p1( min_x - 5.0 * h, min_x - 5.0 * h, min_x - 5.0 * h );
      RealPoint p2( max_x + 5.0 * h, max_x + 5.0 * h, max_x + 5.0 * h );
      CountedPtr<ImplicitDigitalShape> dshape( new ImplicitDigitalShape() );
      dshape->attach( *shape );
      dshape->init( p1, p2, h );
      Domain domain = dshape->getDomain();
      if ( ! K.init( domain.lowerBound(), domain.upperBound(), true ) )
	trace.error() << "[EstimatorHelpers::makeImplicitDigitalShapeFromImplicitShape]"
		      << " Error building Khalimsky space K=" << K << std::endl;
      return dshape;
    }

    /// Vectorizes an implicitly defined digital shape into a binary image.
    ///
    /// @param[in] dshape a smart pointer on an implicit digital shape.
    /// @return a smart pointer on a binary image that samples the digital shape.
    static CountedPtr<BinaryImage>
    makeBinaryImageFromImplicitDigitalShape
    ( CountedPtr<ImplicitDigitalShape> dshape )
    {
      const Domain shapeDomain    = dshape->getDomain();
      CountedPtr<BinaryImage> img ( new BinaryImage( shapeDomain ) );
      std::transform( shapeDomain.begin(), shapeDomain.end(),
		      img->begin(),
		      [dshape] ( const Point& p ) { return (*dshape)(p); } );
      return img;
    }

    /// Vectorizes an implicitly defined digital shape into a binary
    /// image with added Kanungo noise.
    ///
    /// @param[in] dshape a smart pointer on an implicit digital shape.
    ///
    /// @param[in] noiseLevel the parameter that tunes the amount of
    /// Kanungo noise (0.0 is none, 0.5 is quite a lot, 1.0 means
    /// random image).
    ///
    /// @return a smart pointer on a binary image that samples the
    /// digital shape with added noise.
    static CountedPtr<BinaryImage>
    makeNoisyBinaryImageFromImplicitDigitalShape
    ( CountedPtr<ImplicitDigitalShape> dshape,
      Scalar noiseLevel )
    {
      typedef KanungoNoise< ImplicitDigitalShape, Domain > KanungoPredicate;
      const Domain shapeDomain    = dshape->getDomain();
      CountedPtr<BinaryImage> img ( new BinaryImage( shapeDomain ) );
      KanungoPredicate noisy_dshape( *dshape, shapeDomain, noiseLevel );
      std::transform( shapeDomain.begin(), shapeDomain.end(),
		      img->begin(),
		      [&noisy_dshape] ( const Point& p ) { return noisy_dshape(p); } );
      return img;
    }

    /// Adds Kanungo noise to a binary image.
    ///
    /// @param[in] bimage a smart pointer on a binary image.
    ///
    /// @param[in] noiseLevel the parameter that tunes the amount of
    /// Kanungo noise (0.0 is none, 0.5 is quite a lot, 1.0 means
    /// random image).
    ///
    /// @return a smart pointer on the noisified binary image.
    static CountedPtr<BinaryImage>
    makeNoisyBinaryImageFromBinaryImage
    ( CountedPtr<BinaryImage> bimage,
      Scalar noiseLevel )
    {
      typedef KanungoNoise< BinaryImage, Domain > KanungoPredicate;
      const Domain shapeDomain    = bimage->domain();
      CountedPtr<BinaryImage> img ( new BinaryImage( shapeDomain ) );
      KanungoPredicate noisy_dshape( *bimage, shapeDomain, noiseLevel );
      std::transform( shapeDomain.begin(), shapeDomain.end(),
		      img->begin(),
		      [&noisy_dshape] ( const Point& p ) { return noisy_dshape(p); } );
      return img;
    }

    /// Vectorizes an implicitly defined digital shape into a binary
    /// image, and possibly add Kanungo noise to the result depending
    /// on parameters given in \a vm.
    ///
    /// @param[in] vm the options sets in the variable map (arguments
    /// given to the program). Recognized parameters are given in \ref
    /// optionsNoisyImage.
    ///
    /// @param[in] dshape a smart pointer on an implicit digital shape.
    ///
    /// @return a smart pointer on a binary image that samples the
    /// digital shape with added noise.
    static CountedPtr<BinaryImage>
    makeNoisyOrNotBinaryImageFromImplicitDigitalShape
    ( const po::variables_map& vm,
      CountedPtr<ImplicitDigitalShape> dshape )
    {
      Scalar noiseLevel = vm[ "noise" ].as<Scalar>();
      if ( noiseLevel == 0.0 )
	return makeBinaryImageFromImplicitDigitalShape( dshape );
      else
	return makeNoisyBinaryImageFromImplicitDigitalShape( dshape, noiseLevel );
    }

    /// Loads a vol file (specified with option --input) and returns
    /// the corresponding binary image.
    ///
    /// @param[in] vm the options sets in the variable map (arguments
    /// given to the program). Recognized parameters are given in \ref
    /// optionsNoisyImage.
    ///
    /// @param[out] outputs the Khalimsky space whose domain contains
    /// the loaded vol file.
    ///
    /// @return a smart pointer on a binary image that represents the
    /// (thresholded) vol file.
    static CountedPtr<BinaryImage>
    makeBinaryImageFromVolFile
    ( const po::variables_map& vm,
      KSpace&                  K )
    {
      string inputFilename = vm["input"].as<std::string>();
      int thresholdMin     = vm["thresholdMin"].as<int>();
      int thresholdMax     = vm["thresholdMax"].as<int>();
      GrayScaleImage image = GenericReader<GrayScaleImage>::import( inputFilename );
      Domain        domain = image.domain();
      bool        space_ok = K.init( domain.lowerBound(), domain.upperBound(), true );
      if ( ! space_ok )
	{
	  trace.error() << "[EstimatorHelpers::makeBinaryImageFromVolFile]"
			<< " error in space initialization." << std::endl;
	  return CountedPtr<BinaryImage>( nullptr );
	}
      typedef functors::IntervalForegroundPredicate<GrayScaleImage> ThresholdedImage;
      ThresholdedImage tImage( image, thresholdMin, thresholdMax );
      CountedPtr<BinaryImage> img ( new BinaryImage( domain ) );
      std::transform( domain.begin(), domain.end(),
		      img->begin(),
		      [tImage] ( const Point& p ) { return tImage(p); } );
      return img;
    }

    /// Loads a vol file (specified with option --input) and returns
    /// the corresponding binary image.
    ///
    /// @param[in] vm the options sets in the variable map (arguments
    /// given to the program). Recognized parameters are given in \ref
    /// optionsNoisyImage.
    ///
    /// @param[out] outputs the Khalimsky space whose domain contains
    /// the loaded vol file.
    ///
    /// @return a smart pointer on a binary image that represents the
    /// (thresholded) vol file.
    static CountedPtr<BinaryImage>
    makeNoisyOrNotBinaryImageFromVolFile
    ( const po::variables_map& vm,
      KSpace&                  K )
    {
      Scalar noiseLevel = vm[ "noise" ].as<Scalar>();
      if ( noiseLevel == 0.0 )
	return makeBinaryImageFromVolFile( vm, K );
      else {
	KSpace tmp_K;
	auto   tmp_bimage = makeBinaryImageFromVolFile( vm, tmp_K );
	auto   tmp_domain = tmp_bimage->domain();
	Domain ext_domain( tmp_domain.lowerBound() - Point::diagonal( 5 ),
			   tmp_domain.upperBound() + Point::diagonal( 5 ) );
	CountedPtr<BinaryImage> bimage ( new BinaryImage( ext_domain ) );
	unsigned int nb0 = 0;
	unsigned int nb1 = 0;
	for ( auto  p : tmp_domain )
	  if ( (*tmp_bimage)( p ) ) { bimage->setValue( p, true  ); nb1++; }
	  else                      { bimage->setValue( p, false ); nb0++; }
	trace.info() << "input image: nb0=" << nb0 << " nb1=" << nb1 << std::endl;
	bool space_ok = K.init( ext_domain.lowerBound(),
				ext_domain.upperBound(), true );
	if ( ! space_ok )
	  {
	    trace.error() << "[EstimatorHelpers::makeNoisyOrNotBinaryImageFromVolFile]"
			  << " error in space initialization." << std::endl;
	  }
	return makeNoisyBinaryImageFromBinaryImage( bimage, noiseLevel );
      }
    }
    
    /// Builds a digital surface from a space \a K and a binary image \a bimage.
    ///
    /// @param[in] K the Khalimsky space whose domain encompasses the digital shape.
    /// @param[in] bimage a binary image representing the characteristic function of a digital shape.
    ///
    /// @return a smart pointer on a (light) digital surface that
    /// represents the boundary of the digital shape (at least a big
    /// component).
    static CountedPtr<Surface>
    makeDigitalSurfaceFromBinaryImage
    ( const KSpace& K,
      CountedPtr<BinaryImage> bimage )
    {
      SurfelAdjacency< KSpace::dimension > surfAdj( true );

      // We have to search for a surfel that belong to a big connected component.
      CountedPtr<Surface> ptrSurface;
      Surfel              bel;
      Scalar              minsize    = bimage->extent().norm();
      unsigned int        nb_surfels = 0;
      unsigned int        tries      = 0;
      do {
        try { // Search initial bel
          bel = Surfaces<KSpace>::findABel( K, *bimage, 100000 );
        } catch (DGtal::InputException e) {
          trace.error() << "[EstimatorHelpers::makeDigitalSurfaceFromBinaryImage]"
			<< " ERROR Unable to find bel." << std::endl;
          return ptrSurface;
        }
	// this pointer will be acquired by the surface.
        SurfaceContainer* surfContainer = new SurfaceContainer( K, *bimage, surfAdj, bel );
        ptrSurface = CountedPtr<Surface>( new Surface( surfContainer ) ); // acquired
        nb_surfels = ptrSurface->size();
      } while ( ( nb_surfels < 2 * minsize ) && ( tries++ < 150 ) );
      if( tries >= 150 ) {
	trace.warning() << "[EstimatorHelpers::makeDigitalSurfaceFromBinaryImage]"
			<< "ERROR cannot find a proper bel in a big enough component."
			<< std::endl;
      }
      return ptrSurface;
    }

    /// Given a digital surface, returns a vector of surfels in a depth-first
    /// traversal order (best for Integral Invariant).
    ///
    /// @param[in] surface a smort pointer on a digital surface.
    /// @return a range of surfels as a vector.
    static std::vector< Surfel >
    computeDepthFirstSurfelRange( CountedPtr<Surface> surface )
    {
      typedef DepthFirstVisitor< Surface > Visitor;
      typedef GraphVisitorRange< Visitor > VisitorRange;
      std::vector< Surfel > result;
      VisitorRange range( new Visitor( *surface, *( surface->begin() ) ) );
      std::for_each( range.begin(), range.end(),
		     [&result] ( Surfel s ) { result.push_back( s ); } );
      return result;
    }
       

    /// Given a space \a K, an implicit \a shape, a sequence of \a
    /// surfels, and a gridstep \a h, returns the true normals at the
    /// specified surfels, in the same order.
    ///
    /// @param[in] K the Khalimsky space whose domain encompasses the digital shape.
    /// @param[in] shape the implicit shape.
    /// @param[in] h the grid step to embed surfels.
    /// @param[in] surfels the sequence of surfels at which we compute the normals
    ///
    /// @return the vector containing the true normals, in the same
    /// order as \a surfels.
    static std::vector< RealVector >
    computeTrueNormals( const KSpace&                K,
			CountedPtr<ImplicitShape>    shape,
			const Scalar                 h,
			const std::vector< Surfel >& surfels )
    {
      std::vector< RealVector > n_true_estimations;
      TrueNormalEstimator       true_estimator;
      true_estimator.attach( *shape );
      true_estimator.setParams( K, NormalFunctor(), 20, 0.0001, 0.5 );
      true_estimator.init( h, surfels.begin(), surfels.end() );
      true_estimator.eval( surfels.begin(), surfels.end(),
			   std::back_inserter( n_true_estimations ) );
      return n_true_estimations;
    }

    /// Given a space \a K, an implicit \a shape, a sequence of \a
    /// surfels, and a gridstep \a h, returns the true mean curvatures at the
    /// specified surfels, in the same order.
    ///
    /// @param[in] K the Khalimsky space whose domain encompasses the digital shape.
    /// @param[in] shape the implicit shape.
    /// @param[in] h the grid step to embed surfels.
    /// @param[in] surfels the sequence of surfels at which we compute the mean curvatures.
    ///
    /// @return the vector containing the true mean curvatures, in the same
    /// order as \a surfels.
    static std::vector< Scalar >
    computeMeanCurvatures( const KSpace&                K,
			   CountedPtr<ImplicitShape>    shape,
			   const Scalar                 h,
			   const std::vector< Surfel >& surfels )
    {
      std::vector< Scalar >      n_true_estimations;
      TrueMeanCurvatureEstimator true_estimator;
      true_estimator.attach( *shape );
      true_estimator.setParams( K, MeanCurvatureFunctor(), 20, 0.0001, 0.5 );
      true_estimator.init( h, surfels.begin(), surfels.end() );
      true_estimator.eval( surfels.begin(), surfels.end(),
			   std::back_inserter( n_true_estimations ) );
      return n_true_estimations;
    }


    /// Given a space \a K, an implicit \a shape, a sequence of \a
    /// surfels, and a gridstep \a h, returns the true gaussian curvatures at the
    /// specified surfels, in the same order.
    ///
    /// @param[in] K the Khalimsky space whose domain encompasses the digital shape.
    /// @param[in] shape the implicit shape.
    /// @param[in] h the grid step to embed surfels.
    /// @param[in] surfels the sequence of surfels at which we compute the gaussian curvatures.
    ///
    /// @return the vector containing the true gaussian curvatures, in the same
    /// order as \a surfels.
    static std::vector< Scalar >
    computeGaussianCurvatures( const KSpace&                K,
			       CountedPtr<ImplicitShape>    shape,
			       const Scalar                 h,
			       const std::vector< Surfel >& surfels )
    {
      std::vector< Scalar >          n_true_estimations;
      TrueGaussianCurvatureEstimator true_estimator;
      true_estimator.attach( *shape );
      true_estimator.setParams( K, GaussianCurvatureFunctor(), 20, 0.0001, 0.5 );
      true_estimator.init( h, surfels.begin(), surfels.end() );
      true_estimator.eval( surfels.begin(), surfels.end(),
			   std::back_inserter( n_true_estimations ) );
      return n_true_estimations;
    }
    
    /// Given a digital space \a K and a sequence of \a surfels,
    /// returns the trivial normals at the specified surfels, in the
    /// same order.
    ///
    /// @param[in] K the Khalimsky space whose domain encompasses the digital shape.
    /// @param[in] surfels the sequence of surfels at which we compute the normals
    ///
    /// @return the vector containing the estimated normals, in the
    /// same order as \a surfels.
    static std::vector< RealVector >
    computeTrivialNormals( const KSpace&                K,
			   const std::vector< Surfel >& surfels )
    {
      std::vector< RealVector > result;
      for ( auto s : surfels )
	{
	  Dimension  k = K.sOrthDir( s );
	  bool  direct = K.sDirect( s, k );
	  RealVector t = RealVector::zero;
	  t[ k ]       = direct ? -1.0 : 1.0;
	  result.push_back( t );
	}
      return result;
    }

    /// Given a digital surface \a surface, a sequence of \a surfels,
    /// and some parameters \a vm, returns the normal Voronoi
    /// Covariance Measure (VCM) estimation at the specified surfels,
    /// in the same order.
    ///
    /// @param[in] vm the options sets in the variable map (arguments
    /// given to the program). Recognized parameters are given in \ref
    /// optionsNormalEstimators (trivial-ring)
    /// @param[in] surface the digital surface
    /// @param[in] surfels the sequence of surfels at which we compute the normals
    ///
    /// @return the vector containing the estimated normals, in the
    /// same order as \a surfels.
    static std::vector< RealVector >
    computeConvolvedTrivialNormals( const po::variables_map&     vm,
				    CountedPtr<Surface>          surface,
				    const std::vector< Surfel >& surfels )
    {
      Scalar t = vm[ "trivial-ring" ].as<double>();
      typedef ExactPredicateLpSeparableMetric<Space,2>              Metric;
      typedef functors::HatFunction<Scalar>                         Functor;
      typedef functors::ElementaryConvolutionNormalVectorEstimator
	< Surfel, CanonicSCellEmbedder<KSpace> >                    SurfelFunctor;
      typedef LocalEstimatorFromSurfelFunctorAdapter
	< SurfaceContainer, Metric, SurfelFunctor, Functor>         NormalEstimator;
      trace.info() << " CTrivial normal t=" << t << " (discrete)" << std::endl;
      const Functor fct( 1.0, t );
      const KSpace &  K = surface->container().space();
      Metric    aMetric;
      CanonicSCellEmbedder<KSpace> canonic_embedder( K );
      std::vector< RealVector >    n_estimations;
      SurfelFunctor                surfelFct( canonic_embedder, 1.0 );
      NormalEstimator              estimator;
      estimator.attach( *surface);
      estimator.setParams( aMetric, surfelFct, fct, t );
      estimator.init( 1.0, surfels.begin(), surfels.end());
      estimator.eval( surfels.begin(), surfels.end(),
		      std::back_inserter( n_estimations ) );
      std::transform( n_estimations.cbegin(), n_estimations.cend(), n_estimations.begin(),
		      [] ( RealVector v ) { return -v; } );
      return n_estimations;
    }
    
    /// Given a digital surface \a surface, a sequence of \a surfels,
    /// and some parameters \a vm, returns the normal Voronoi
    /// Covariance Measure (VCM) estimation at the specified surfels,
    /// in the same order.
    ///
    /// @param[in] vm the options sets in the variable map (arguments
    /// given to the program). Recognized parameters are given in \ref
    /// optionsNormalEstimators.
    /// @param[in] surface the digital surface
    /// @param[in] surfels the sequence of surfels at which we compute the normals
    ///
    /// @return the vector containing the estimated normals, in the
    /// same order as \a surfels.
    static std::vector< RealVector >
    computeVCMNormals( const po::variables_map&     vm,
		       CountedPtr<Surface>          surface,
		       const std::vector< Surfel >& surfels )
    {
      typedef ExactPredicateLpSeparableMetric<Space,2> Metric;
      std::vector< RealVector > n_estimations;
      std::string kernel = vm[ "kernel" ].as<std::string>();
      Scalar      h      = vm[ "gridstep" ].as<Scalar>();
      Scalar      R      = vm[ "R-radius" ].as<Scalar>();
      Scalar      r      = vm[ "r-radius" ].as<Scalar>();
      Scalar      t      = vm[ "trivial-ring" ].as<Scalar>();
      Scalar      alpha  = vm[ "alpha" ].as<Scalar>();
      int      embedding = vm[ "embedding" ].as<int>();
      // Adjust parameters according to gridstep if specified.
      if ( alpha != 1.0 ) R *= pow( h, alpha-1.0 );
      if ( alpha != 1.0 ) r *= pow( h, alpha-1.0 );
      Surfel2PointEmbedding embType = embedding == 0 ? Pointels :
                                      embedding == 1 ? InnerSpel : OuterSpel;
      trace.info() << "- VCM normal kernel=" << kernel << " emb=" << embedding
		   << " alpha=" << alpha << std::endl;
      trace.info() << "- VCM normal r=" << (r*h)  << " (continuous) "
		   << r << " (discrete)" << std::endl;
      trace.info() << "- VCM normal R=" << (R*h)  << " (continuous) "
		   << R << " (discrete)" << std::endl;
      trace.info() << "- VCM normal t=" << t << " (discrete)" << std::endl;
      if ( kernel == "hat" ) {
	typedef functors::HatPointFunction<Point,Scalar>             KernelFunction;
	typedef VoronoiCovarianceMeasureOnDigitalSurface
	  < SurfaceContainer, Metric, KernelFunction >               VCMOnSurface;
	typedef functors::VCMNormalVectorFunctor<VCMOnSurface>       NormalFunctor;
	typedef VCMDigitalSurfaceLocalEstimator
	  < SurfaceContainer, Metric, KernelFunction, NormalFunctor> VCMNormalEstimator;
	KernelFunction chi_r( 1.0, r );
	VCMNormalEstimator estimator;
	estimator.attach( *surface );
	estimator.setParams( embType, R, r, chi_r, t, Metric(), true );
	estimator.init( h, surfels.begin(), surfels.end() );
	estimator.eval( surfels.begin(), surfels.end(),
			std::back_inserter( n_estimations ) );
      } else if ( kernel == "ball" ) {
	typedef functors::BallConstantPointFunction<Point,Scalar>    KernelFunction;
	typedef VoronoiCovarianceMeasureOnDigitalSurface
	  < SurfaceContainer, Metric, KernelFunction >               VCMOnSurface;
	typedef functors::VCMNormalVectorFunctor<VCMOnSurface>       NormalFunctor;
	typedef VCMDigitalSurfaceLocalEstimator
	  < SurfaceContainer, Metric, KernelFunction, NormalFunctor> VCMNormalEstimator;
	KernelFunction chi_r( 1.0, r );
	VCMNormalEstimator estimator;
	estimator.attach( *surface );
	estimator.setParams( embType, R, r, chi_r, t, Metric(), true );
	estimator.init( h, surfels.begin(), surfels.end() );
	estimator.eval( surfels.begin(), surfels.end(),
			std::back_inserter( n_estimations ) );
      }
      return n_estimations;
    }

    /// Given a digital shape \a bimage, a sequence of \a surfels,
    /// and some parameters \a vm, returns the normal Integral
    /// Invariant (VCM) estimation at the specified surfels, in the
    /// same order.
    ///
    /// @param[in] vm the options sets in the variable map (arguments
    /// given to the program). Recognized parameters are given in \ref
    /// optionsNormalEstimators.
    /// @param[in] K the digital space where the shape lives.
    /// @param[in] bimage the characteristic function of the shape as a binary image (inside is true, outside is false).
    /// @param[in] surfels the sequence of surfels at which we compute the normals
    ///
    /// @return the vector containing the estimated normals, in the
    /// same order as \a surfels.
    ///
    /// @note It is better to have surfels in a specific order, as
    /// given for instance by computeDepthFirstSurfelRange.
    static std::vector< RealVector >
    computeIINormals( const po::variables_map&     vm,
		      const KSpace&                K,
		      CountedPtr<BinaryImage>      bimage,
		      const std::vector< Surfel >& surfels )
    {
      typedef functors::IINormalDirectionFunctor<Space> IINormalFunctor;
      typedef IntegralInvariantCovarianceEstimator
	<KSpace, BinaryImage, IINormalFunctor>          IINormalEstimator;
      std::vector< RealVector > n_estimations;
      Scalar h     = vm[ "gridstep" ].as<Scalar>();
      Scalar r     = vm[ "r-radius" ].as<Scalar>();
      Scalar alpha = vm[ "alpha"    ].as<Scalar>();
      if ( alpha != 1.0 ) r *= pow( h, alpha-1.0 );
      trace.info() << " II normal alpha=" << alpha << std::endl;
      trace.info() << " II normal r=" << (r*h)  << " (continuous) "
		   << r << " (discrete)" << std::endl;
      IINormalFunctor     functor;
      functor.init( h, r*h );
      IINormalEstimator   ii_estimator( functor );
      ii_estimator.attach( K, *bimage );
      ii_estimator.setParams( r );
      ii_estimator.init( h, surfels.begin(), surfels.end() );
      ii_estimator.eval( surfels.begin(), surfels.end(),
			 std::back_inserter( n_estimations ) );
      const std::vector< RealVector > n_trivial = computeTrivialNormals( K, surfels );
      orientVectors( n_trivial, n_estimations );
      return n_estimations;
    }

    /// Given a digital surface \a surface, its corresponding
    /// characteristic function as a binary image \a bimage, a
    /// sequence of \a surfels, and some parameters \a vm, returns
    /// some normal estimation at the specified surfels, in the same
    /// order. The estimator is chosen according to parameter
    /// "-estimator" and can be any of "True", "Trivial", "II" or "VCM".
    ///
    /// @param[in] vm the options sets in the variable map (arguments
    /// given to the program). Recognized parameters are given in \ref
    /// optionsNormalEstimators.
    /// @param[in] shape the implicit shape (if 0, "True" is not possible).
    /// @param[in] surface the digital surface (if 0, "VCM" is not possible).
    /// @param[in] bimage the characteristic function of the shape as a binary image (inside is true, outside is false) (if 0, "II" is not possible).
    /// @param[in] surfels the sequence of surfels at which we compute the normals
    ///
    /// @return the vector containing the estimated normals, in the
    /// same order as \a surfels.
    static std::vector< RealVector >
    computeNormals( const po::variables_map&     vm,
		    const KSpace&                K,
		    CountedPtr<ImplicitShape>    shape,
		    CountedPtr<Surface>          surface,
		    CountedPtr<BinaryImage>      bimage,
		    const std::vector< Surfel >& surfels )
    {
      std::string estimator = vm[ "estimator" ].as<std::string>();
      const Scalar        h = vm[ "gridstep" ].as<Scalar>();
      if ( estimator == "True" )
	{
	  if ( shape != 0 ) return computeTrueNormals( K, shape, h, surfels );
	  trace.error() << "[EstimatorHelpers::computeNormals]"
			<< " Error. Cannot estimate true normals if implicit shape is not given."
			<< std::endl;
	}
      else if ( estimator == "Trivial" )
	return computeTrivialNormals( K, surfels );
      else if ( estimator == "CTrivial" )
	{
	  if ( surface != 0 ) return computeConvolvedTrivialNormals( vm, surface, surfels );
	  trace.error() << "[EstimatorHelpers::computeNormals]"
			<< " Error. Cannot estimate CTrivial normals if digital surface is not given."
			<< std::endl;
	}
      else if ( estimator == "VCM" )
	{
	  if ( surface != 0 ) return computeVCMNormals( vm, surface, surfels );
	  trace.error() << "[EstimatorHelpers::computeNormals]"
			<< " Error. Cannot estimate VCM normals if digital surface is not given."
			<< std::endl;
	}
      else if ( estimator == "II" )
	{
	  if ( bimage != 0 )
	    {
	      const Scalar               t = vm[ "trivial-ring" ].as<Scalar>();
	      std::vector< RealVector > ii = computeIINormals( vm, K, bimage, surfels );
	      if ( ( t >= 1.0 ) && ( surface != 0 ) )
		{ // Forces realignment.
		  std::vector< RealVector > cn =
		    computeConvolvedTrivialNormals( vm, surface, surfels );
		  orientVectors( cn, ii );
		}
	      return ii;
	    }
	  trace.error() << "[EstimatorHelpers::computeNormals]"
			<< " Error. Cannot estimate II normals if characteristic function is not given."
			<< std::endl;
	}
      else
	{
	  trace.error() << "[EstimatorHelpers::computeNormals]"
			<< " Error. Unknown \"" << estimator << "\". Should be one of True|VCM|II|Trivial."
			<< std::endl;
	}
      return std::vector< RealVector >();
    }


    /// Given a digital shape \a bimage, a sequence of \a surfels,
    /// and some parameters \a vm, returns the mean curvature Integral
    /// Invariant (VCM) estimation at the specified surfels, in the
    /// same order.
    ///
    /// @param[in] vm the options sets in the variable map (arguments
    /// given to the program). Recognized parameters are given in \ref
    /// optionsNormalEstimators.
    /// @param[in] K the digital space where the shape lives.
    /// @param[in] bimage the characteristic function of the shape as a binary image (inside is true, outside is false).
    /// @param[in] surfels the sequence of surfels at which we compute the normals
    ///
    /// @return the vector containing the estimated mean curvatures, in the
    /// same order as \a surfels.
    ///
    /// @note It is better to have surfels in a specific order, as
    /// given for instance by computeDepthFirstSurfelRange.
    static std::vector< Scalar >
    computeIIMeanCurvatures( const po::variables_map&     vm,
			     const KSpace&                K,
			     CountedPtr<BinaryImage>      bimage,
			     const std::vector< Surfel >& surfels )
    {
      typedef functors::IIMeanCurvature3DFunctor<Space> IIMeanCurvFunctor;
      typedef IntegralInvariantVolumeEstimator
	<KSpace, BinaryImage, IIMeanCurvFunctor>        IIMeanCurvEstimator;
      std::vector< Scalar > mc_estimations;
      Scalar h     = vm[ "gridstep" ].as<Scalar>();
      Scalar r     = vm[ "r-radius" ].as<Scalar>();
      Scalar alpha = vm[ "alpha"    ].as<Scalar>();
      if ( alpha != 1.0 ) r *= pow( h, alpha-1.0 );
      trace.info() << " II Mean curvature alpha=" << alpha << std::endl;
      trace.info() << " II Mean curvature r=" << (r*h)  << " (continuous) "
		   << r << " (discrete)" << std::endl;
      IIMeanCurvFunctor   functor;
      functor.init( h, r*h );
      IIMeanCurvEstimator ii_estimator( functor );
      ii_estimator.attach( K, *bimage );
      ii_estimator.setParams( r );
      ii_estimator.init( h, surfels.begin(), surfels.end() );
      ii_estimator.eval( surfels.begin(), surfels.end(),
			 std::back_inserter( mc_estimations ) );
      return mc_estimations;
    }
    
    /// Given a digital shape \a bimage, a sequence of \a surfels,
    /// and some parameters \a vm, returns the gaussian curvature Integral
    /// Invariant (VCM) estimation at the specified surfels, in the
    /// same order.
    ///
    /// @param[in] vm the options sets in the variable map (arguments
    /// given to the program). Recognized parameters are given in \ref
    /// optionsNormalEstimators.
    /// @param[in] K the digital space where the shape lives.
    /// @param[in] bimage the characteristic function of the shape as a binary image (inside is true, outside is false).
    /// @param[in] surfels the sequence of surfels at which we compute the normals
    ///
    /// @return the vector containing the estimated gaussian curvatures, in the
    /// same order as \a surfels.
    ///
    /// @note It is better to have surfels in a specific order, as
    /// given for instance by computeDepthFirstSurfelRange.
    static std::vector< Scalar >
    computeIIGaussianCurvatures( const po::variables_map&     vm,
				 const KSpace&                K,
				 CountedPtr<BinaryImage>      bimage,
				 const std::vector< Surfel >& surfels )
    {
      typedef functors::IIGaussianCurvature3DFunctor<Space> IIGaussianCurvFunctor;
      typedef IntegralInvariantCovarianceEstimator
	<KSpace, BinaryImage, IIGaussianCurvFunctor>        IIGaussianCurvEstimator;
      std::vector< Scalar > mc_estimations;
      Scalar h     = vm[ "gridstep" ].as<Scalar>();
      Scalar r     = vm[ "r-radius" ].as<Scalar>();
      Scalar alpha = vm[ "alpha"    ].as<Scalar>();
      if ( alpha != 1.0 ) r *= pow( h, alpha-1.0 );
      trace.info() << " II Gauss curvature alpha=" << alpha << std::endl;
      trace.info() << " II Gauss curvature r=" << (r*h) << " (continuous) "
		   << r << " (discrete)" << std::endl;
      IIGaussianCurvFunctor   functor;
      functor.init( h, r*h );
      IIGaussianCurvEstimator ii_estimator( functor );
      ii_estimator.attach( K, *bimage );
      ii_estimator.setParams( r );
      ii_estimator.init( h, surfels.begin(), surfels.end() );
      ii_estimator.eval( surfels.begin(), surfels.end(),
			 std::back_inserter( mc_estimations ) );
      return mc_estimations;
    }
    
    
    /// Orient \a v so that it points in the same direction as \a
    /// ref_v (scalar product is then non-negative afterwards).
    ///
    /// @param[in]    ref_v the vectors having the reference orientation.
    /// @param[inout] v the vectors to reorient.
    static void
    orientVectors( const std::vector< RealVector > & ref_v,
		   std::vector< RealVector > &       v )
    {
      std::transform( ref_v.cbegin(), ref_v.cend(), v.cbegin(), v.begin(), 
		      [] ( RealVector rw, RealVector w )
		      { return rw.dot( w ) >= 0.0 ? w : -w; } );
    }

    /// Computes the statistic that measures the angle differences
    /// between the two arrays of unit vectors.
    ///
    /// @param[in] v1 the first array of unit vectors (normals)
    /// @param[in] v2 the second array of unit vectors (normals)
    /// @return their angle difference as a statistic.
    static AngleDevStatistic
    measureAngleDeviation( const std::vector< RealVector > & v1,
			   const std::vector< RealVector > & v2 )
    {
      AngleDevStatistic stat;
      for ( auto it1 = v1.cbegin(), it2 = v2.cbegin(), itE1 = v1.cend();
	    it1 != itE1; ++it1, ++it2 )
	{
          Scalar angle_error = acos( (*it1).dot( *it2 ) );
          stat.addValue( angle_error );
	}
      stat.terminate();
      return stat;
    }

    /// Computes the absolute difference between each element of the two vectors.
    /// @param[in] v1 any vector of values.
    /// @param[in] v2 any vector of values.
    /// @return the vector composed of elemenst |v1[i]-v2[i]|.
    static std::vector< Scalar >
    absoluteDifference( const std::vector< Scalar > & v1,
			const std::vector< Scalar > & v2 )
    {
      std::vector< Scalar > result( v1.size() );
      std::transform( v2.cbegin(), v2.cend(), v1.cbegin(), result.begin(), 
		      [] ( Scalar val1, Scalar val2 )
		      { return fabs( val1 - val2 ); } );
      return result;
    }

    /// Computes the l2-norm of v1-v2, ie the square root of the
    /// mean-squared error of the two vectors.
    ///
    /// @param[in] v1 any vector of values.
    /// @param[in] v2 any vector of values.
    /// @return the normL2 of v1-v2, ie. sqrt( 1/n sum_i (v1[i]-v2[i])^2 ).
    static Scalar
    normL2( const std::vector< Scalar > & v1,
	    const std::vector< Scalar > & v2 )
    {
      Scalar sum = 0;
      for ( unsigned int i = 0; i < v1.size(); i++ )
	sum += ( v1[ i ] - v2[ i ] ) * ( v1[ i ] - v2[ i ] );
      return sqrt( sum / v1.size() );
    }

    /// Computes the l1-norm of v1-v2, ie the average of the absolute
    /// differences of the two vectors.
    ///
    /// @param[in] v1 any vector of values.
    /// @param[in] v2 any vector of values.
    /// @return the normL1 of v1-v2, ie. 1/n sum_i |v1[i]-v2[i]|.
    static Scalar
    normL1( const std::vector< Scalar > & v1,
	    const std::vector< Scalar > & v2 )
    {
      Scalar sum = 0;
      for ( unsigned int i = 0; i < v1.size(); i++ )
	sum += fabs( v1[ i ] - v2[ i ] );
      return sum / v1.size();
    }

    /// Computes the loo-norm of v1-v2, ie the maximum of the absolute
    /// differences of the two vectors.
    ///
    /// @param[in] v1 any vector of values.
    /// @param[in] v2 any vector of values.
    /// @return the normLoo of v1-v2, ie. max_i |v1[i]-v2[i]|.
    static Scalar
    normLoo( const std::vector< Scalar > & v1,
	     const std::vector< Scalar > & v2 )
    {
      Scalar loo = 0;
      for ( unsigned int i = 0; i < v1.size(); i++ )
	loo = std::max( loo, fabs( v1[ i ] - v2[ i ] ) );
      return loo;
    }
    
#ifdef WITH_VISU3D_QGLVIEWER
    static ColorMap
    makeColorMap( const po::variables_map& vm,
		  Scalar                   min,
		  Scalar                   max )
    {
      std::string cmap = vm[ "colormap" ].as<std::string>();
      if      ( cmap == "Cool" )   return ColorMap( min, max, CMAP_COOL );
      else if ( cmap == "Copper" ) return ColorMap( min, max, CMAP_COPPER );
      else if ( cmap == "Hot" )    return ColorMap( min, max, CMAP_HOT );
      else if ( cmap == "Jet" )    return ColorMap( min, max, CMAP_JET );
      else if ( cmap == "Spring" ) return ColorMap( min, max, CMAP_SPRING );
      else if ( cmap == "Summer" ) return ColorMap( min, max, CMAP_SUMMER );
      else if ( cmap == "Autumn" ) return ColorMap( min, max, CMAP_AUTUMN );
      else if ( cmap == "Winter" ) return ColorMap( min, max, CMAP_WINTER );
      // Custom
      ColorMap gradcmap( min, max );
      gradcmap.addColor( Color( 0, 0, 255 ) );
      gradcmap.addColor( Color( 0, 255, 255 ) );
      gradcmap.addColor( Color( 255, 255, 255 ) );
      gradcmap.addColor( Color( 255, 255, 0 ) );
      gradcmap.addColor( Color( 255, 0, 0 ) );
      return gradcmap;
    }

    static ColorMap
    makeTicsColorMap( Scalar                   min,
		      Scalar                   max )
    {
      ColorMap gradcmap( min, max );
      gradcmap.addColor( Color( 0, 0, 155 ) );
      gradcmap.addColor( Color( 0, 155, 155 ) );
      gradcmap.addColor( Color( 155, 155, 155 ) );
      gradcmap.addColor( Color( 155, 155, 0 ) );
      gradcmap.addColor( Color( 155, 0, 0 ) );
      return gradcmap;
    }

    static void
    viewSurfelValues( Viewer& viewer,
		      const po::variables_map&         vm,
		      const std::vector< Surfel >&     surfels,
		      const std::vector< Scalar >&     values,
		      const std::vector< RealVector >& normals )
    {
      Scalar      m = vm.count( "minValue" )
	? vm[ "minValue" ].as<double>()
	: * std::min_element( values.begin(), values.end() );
      Scalar      M = vm.count( "maxValue" )
	? vm[ "maxValue" ].as<double>()
	: * std::max_element( values.begin(), values.end() );
      // std::cout << "m=" << m << " M=" << M << std::endl;
      ColorMap cmap = makeColorMap( vm, m, M );
      // ColorMap cmap = ColorMap( m, M, CMAP_COOL ); 
      // std::cout << "colormap" << std::endl;
      auto      itV  = values.cbegin(); 
      auto      itN  = normals.cbegin(); 
      auto      itNE = normals.cend(); 
      // std::cout << "first value=" << *itV << std::endl;
      Surfel  dummy;
      viewer << SetMode3D( dummy.className(), "Basic" );
      for ( auto it = surfels.cbegin(), itE = surfels.cend(); it != itE; ++it )
	{
	  Scalar          v = *itV++;
	  // std::cout << " " << v << std::endl;
	  v = std::min( M, std::max( m, v ) );
          viewer.setFillColor( cmap( v ) );
	  // Display surfel with or without normal estimation.
	  if ( itN != itNE )
	    DisplayFactory::drawOrientedSurfelWithNormal( viewer, *it, *itN++, false );
	  else
	    viewer << *it;
	}
    }
    
    static void
    viewSurfaceIsolines( Viewer& viewer,
			 const po::variables_map&     vm,
			 CountedPtr<Surface>          surface,
			 const std::vector< Surfel >& surfels,
			 const std::vector< Scalar >& values )
    {
      Scalar      m = vm.count( "minValue" )
	? vm[ "minValue" ].as<double>()
	: * std::min_element( values.begin(), values.end() );
      Scalar      M = vm.count( "maxValue" )
	? vm[ "maxValue" ].as<double>()
	: * std::max_element( values.begin(), values.end() );
      const KSpace& K = surface->container().space();
      // Create map Surfel -> Value
      std::map<Surfel,Scalar> map_values;
      auto      itV = values.begin(); 
      for ( auto it = surfels.begin(), itE = surfels.end(); it != itE; ++it )
	map_values[ *it ] = *itV++;
      
      // Display isolines.
      const Scalar l_zero = vm[ "zero" ].as<double>();
      const Scalar l_tics = vm[ "tics" ].as<double>();
      const ColorMap cmap = makeTicsColorMap( m, M );

      const RealPoint  st = RealPoint::diagonal( -0.5 ); 
      for ( auto it = surfels.begin(), itE = surfels.end(); it != itE; ++it )
	{
	  Vertex       s = *it;
	  Scalar   val_s = map_values[ s ] - l_zero;
	  Integer band_s = (Integer) floor( val_s / l_tics );
	  ArcRange  arcs = surface->outArcs( s );
	  for ( auto a : arcs )
	    {
	      Vertex       t = surface->head( a ); 
	      Scalar   val_t = map_values[ t ] - l_zero;
	      Integer band_t = (Integer) floor( val_t / l_tics );
	      if ( band_s < band_t )
		{
		  Scalar  avg = 0.5 * ( val_s + val_t );
		  avg         = std::min( M, std::max( m, avg ) );
		  Color color = (band_s < 0 && band_t >= 0) ? Color::Black : cmap( avg );
		  viewer.setLineColor( color );
		  viewer.setFillColor( color );
		  Cell   linel = K.unsigns( surface->separator(a) );
		  Dimension  k = * ( K.uDirs( linel ) );
		  Point     p1 = K.uCoords( K.uIncident( linel, k, false ) );
		  Point     p2 = K.uCoords( K.uIncident( linel, k, true ) );
		  viewer.addCylinder( RealPoint( p1[ 0 ], p1[ 1 ], p1[ 2 ] ) + st,
				      RealPoint( p2[ 0 ], p2[ 1 ], p2[ 2 ] ) + st,
				      0.1 );
		}
	    }
	}
    }
#endif
    
  }; // END of class EstimatorHelpers

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined EstimatorHelpers_h

#undef EstimatorHelpers_RECURSES
#endif // else defined(EstimatorHelpers_RECURSES)

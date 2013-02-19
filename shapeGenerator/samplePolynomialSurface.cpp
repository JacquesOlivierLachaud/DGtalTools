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
 * @file samplePolynomialSurface.cpp
 * @ingroup Tools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2012/02/06
 *
 * A simple marching cube algorithm that is used to sample the
 * zero-level of a polynomial surface.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
//! [samplePolynomialSurface-basicIncludes]
#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/SetOfSurfels.h"
#include "DGtal/math/MPolynomial.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"
#include "DGtal/shapes/implicit/ImplicitFunctionDiff1LinearCellEmbedder.h"
#include "DGtal/io/readers/MPolynomialReader.h"
//! [samplePolynomialSurface-basicIncludes]

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace Z3i;

///////////////////////////////////////////////////////////////////////////////



/** 
 * Missing parameter error message.
 * 
 * @param param 
 */
void missingParam(std::string param)
{
  trace.error() << " Parameter: "<<param<<" is required..";
  trace.info() <<std::endl;
  exit(1);
}

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
    ("help,h", "display this message")
    ("polynomial,p", po::value<std::string>(), "The implicit 3d polynomial (like x^2+2*y^2-z^3)")
    ("lower_x,lx",  po::value<double>()->default_value(-1.0), "x-coordinate of lower bound" )
    ("lower_y,ly",  po::value<double>()->default_value(-1.0), "y-coordinate of lower bound" )
    ("lower_z,lz",  po::value<double>()->default_value(-1.0), "z-coordinate of lower bound" )
    ("upper_x,lx",  po::value<double>()->default_value(1.0), "x-coordinate of upper bound" )
    ("upper_y,ly",  po::value<double>()->default_value(1.0), "y-coordinate of upper bound" )
    ("upper_z,lz",  po::value<double>()->default_value(1.0), "z-coordinate of upper bound" )
    ("grid,g",  po::value<double>()->default_value(0.01), "grid step" )
    // ("grid_x,gx",  po::value<double>()->default_value(0.01), "grid step for x-coordinate" )
    // ("grid_y,gy",  po::value<double>()->default_value(0.01), "grid step for y-coordinate" )
    // ("grid_z,gz",  po::value<double>()->default_value(0.01), "grid step for z-coordinate" )
    ("output,o",   po::value<string>()->default_value("surface.xyz"), "Base name of the file containing the surface sampling (format per line is x y z nx ny nz )" )
    ("surface,s",   po::value<string>()->default_value("surface.off"), "Base name of the file containing the surface triangulation (off format)" );

  bool parseOK = true;
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  } catch(const std::exception& ex){
    parseOK = false;
    trace.error() << "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);    
  if ( ! parseOK || vm.count("help") || ( argc <= 1 ) )
    {
      trace.info() << "Generate a cloud of points (a sampling) approximating a implicitly defined polynomial surface. The sampling is obtained with a marching cube algorithm." <<std::endl 
                   << "Basic usage: "<<std::endl
                   << "\t samplePolynomialSurface -p <polynomial> [otherOptions]"<<std::endl
                   << general_opt << "\n";
      return 0;
    }

  //Parse options
  if (!(vm.count("polynomial"))) missingParam("--polynomial");
  std::string polynomialName = vm["polynomial"].as<std::string>();

  //! [samplePolynomialSurface-makeSurface]
  trace.beginBlock( "Making polynomial surface." );
  typedef Space::RealPoint RealPoint;
  typedef Space::RealVector RealVector;
  typedef RealPoint::Coordinate Scalar;
  typedef MPolynomial<3, Scalar> Polynomial3;
  typedef MPolynomialReader<3, Scalar> Polynomial3Reader;
  typedef ImplicitPolynomial3Shape<Space> ImplicitShape;
  typedef GaussDigitizer<Space,ImplicitShape> DigitalShape; 
  typedef DigitalShape::PointEmbedder DigitalEmbedder;

  Polynomial3 P;
  Polynomial3Reader reader;
  std::string::const_iterator iter 
    = reader.read( P, polynomialName.begin(), polynomialName.end() );
  if ( iter != polynomialName.end() )
    {
      trace.error() << "ERROR when reading polynomial: read only <" 
                    << polynomialName.substr( 0, iter - polynomialName.begin() )
                    << ">, created polynomial is P=" << P << std::endl;
      return 2;
    }
  trace.info() << "P( X_0, X_1, X_2 ) = " << P << std::endl;
  ImplicitShape ishape( P );
  DigitalShape dshape;
  dshape.attach( ishape );
  RealPoint p1( vm["lower_x"].as<double>(),
                vm["lower_y"].as<double>(),
                vm["lower_z"].as<double>() );
  RealPoint p2( vm["upper_x"].as<double>(),
                vm["upper_y"].as<double>(),
                vm["upper_z"].as<double>() );
  Scalar step = vm["grid"].as<double>();
  dshape.init( RealPoint( p1 ), RealPoint( p2 ), step );
  Domain domain = dshape.getDomain();
  trace.endBlock();
  //! [samplePolynomialSurface-makeSurface]

  //! [samplePolynomialSurface-KSpace]
  // Construct the Khalimsky space from the image domain
  KSpace K;
  // NB: it is \b necessary to work with a \b closed cellular space
  // since umbrellas use separators and pivots, which must exist for
  // arbitrary surfels.
  bool space_ok = K.init( domain.lowerBound(), 
                          domain.upperBound(), true // necessary
                          );
  if (!space_ok)
    {
      trace.error() << "ERROR in the Khamisky space construction." << std::endl;
      return 3;
    }
  //! [samplePolynomialSurface-KSpace]

  //! [samplePolynomialSurface-SurfelAdjacency]
  typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
  MySurfelAdjacency surfAdj( true ); // interior in all directions.
  //! [samplePolynomialSurface-SurfelAdjacency]


  //! [samplePolynomialSurface-ExtractingSurface]
  trace.beginBlock( "Extracting boundary by scanning the space. " );
  typedef KSpace::SurfelSet SurfelSet;
  typedef SetOfSurfels< KSpace, SurfelSet > MySetOfSurfels;
  typedef DigitalSurface< MySetOfSurfels > MyDigitalSurface;
  typedef MyDigitalSurface::ConstIterator ConstIterator;
  MySetOfSurfels theSetOfSurfels( K, surfAdj );
  Surfaces<KSpace>::sMakeBoundary( theSetOfSurfels.surfelSet(),
                                   K, dshape,
                                   domain.lowerBound(),
                                   domain.upperBound() );
  MyDigitalSurface digSurf( theSetOfSurfels );
  
  trace.info() << "Digital surface has " << digSurf.size() << " surfels."
               << std::endl;
  trace.endBlock();
  //! [samplePolynomialSurface-ExtractingSurface]

  // The cell embedder is used to place vertices closer to the set
  // P(x,y,z)=0
  typedef 
    ImplicitFunctionDiff1LinearCellEmbedder< KSpace, 
                                             ImplicitShape, 
                                             DigitalEmbedder >
    CellEmbedder;
  typedef CellEmbedder::GradientMap GradientMap;
  CellEmbedder cellEmbedder;
  cellEmbedder.init( K, ishape, dshape.pointEmbedder() );
  GradientMap gradient = cellEmbedder.gradientMap();
  //! [samplePolynomialSurface-makingOFF]
  if ( vm.count("surface") )
    {
      trace.beginBlock( "Making OFF surface. " );
      std::string surfaceName = vm["surface"].as<std::string>();
      trace.info() << "- creating OFF surface file <" << surfaceName << ">."
                   << std::endl;
      ofstream out( surfaceName.c_str() );
      if ( out.good() )
        digSurf.exportEmbeddedSurfaceAs3DOFF( out, cellEmbedder );
      out.close();
      trace.endBlock();
    }
  //! [samplePolynomialSurface-makingOFF]
  
  //! [samplePolynomialSurface-makingXYZ]
  if ( vm.count("output") )
    {
      trace.beginBlock( "Making X,Y,Z / NX,NY,NZ sampled surface." );
      std::string outputName = vm["output"].as<std::string>();
      trace.info() << "- creating sample surface file <" << outputName << ">."
                   << std::endl;
      ofstream out( outputName.c_str() );
      if ( out.good() )
        {
          typedef CellEmbedder::Cell Cell;
          typedef DGtal::uint64_t Number;
          Number nbv = digSurf.size();
          // Outputs XYZ header.
          out << "# X Y Z NX NY NZ" << std::endl
              << "# " << nbv << " points." << std::endl
              << "# Generated by DGtal::samplePolynomialSurface." << std::endl;
          // Outputs vertex coordinates (the 3 first ones).
          RealPoint p;
          RealVector v;
          for ( ConstIterator it = digSurf.begin(), it_end = digSurf.end();
                it != it_end; ++it )
            {
              Cell c = K.unsigns( *it );
              p = cellEmbedder( c );
              v = gradient( c );
              double norm = v.norm();
              if ( norm != 0.0 ) v /= norm;
              out << p[ 0 ] << " " << p[ 1 ] << " " << p[ 2 ] << " "
                  << v[ 0 ] << " " << v[ 1 ] << " " << v[ 2 ] << std::endl;
            }
        }
      out.close();
      trace.endBlock();
   }
  //! [samplePolynomialSurface-makingXYZ]

  return 0;
}

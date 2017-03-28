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
 * @file volAntiAliasing.cpp
 * @ingroup Tools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2017/02/28
 *
 * A tool file named volAntiAliasing.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <DGtal/io/readers/VolReader.h>
#include <DGtal/io/writers/GenericWriter.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/IntervalForegroundPredicate.h>
#include <DGtal/topology/SetOfSurfels.h>
#include <DGtal/topology/DigitalSurface.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include <DGtal/math/linalg/EigenSupport.h>
#include <DGtal/dec/DiscreteExteriorCalculus.h>
#include <DGtal/dec/DiscreteExteriorCalculusSolver.h>
#include "NormalEstimation.h"

using namespace std;
using namespace DGtal;
using namespace DGtal::rt;
namespace po = boost::program_options;


namespace DGtal {

  template <typename TKSpace, typename TBooleanImage>
  struct ImplicitDigitalVolume {
    
    typedef TBooleanImage                        BooleanImage;
    typedef TKSpace                              KSpace;
    typedef typename KSpace::SCell               SCell;
    typedef typename KSpace::Cell                Cell;
    typedef typename KSpace::SCellSet            SCellSet;
    typedef typename KSpace::Space               Space;
    typedef typename Space::Point                Point;
    typedef typename Space::Vector               Vector;
    typedef typename Space::RealVector           RealVector;
    typedef typename RealVector::Component       Scalar;
    typedef typename KSpace::Surfel              Surfel;
    typedef SetOfSurfels<KSpace, SCellSet>       SurfaceStorage;
    typedef DigitalSurface< SurfaceStorage >     Surface;
    typedef std::map<SCell,RealVector>           VectorField;
    typedef VectorField                          NormalMap;
    
    typedef EigenLinearAlgebraBackend             LinearAlgebra;
    typedef HyperRectDomain<Space>                Domain;
    typedef DiscreteExteriorCalculus<3,3, LinearAlgebra>   Calculus;
    typedef Calculus::Index                       Index;
    typedef Calculus::DualForm0                   DualForm0;
    typedef Calculus::DualForm1                   DualForm1;
    typedef Calculus::DualIdentity0               DualIdentity0;
    typedef Calculus::PrimalDerivative0           PrimalDerivative0;
    typedef Calculus::PrimalDerivative1           PrimalDerivative1;
    typedef Calculus::DualDerivative0             DualDerivative0;
    typedef Calculus::DualDerivative1             DualDerivative1;
    typedef Calculus::PrimalAntiderivative1       PrimalAntiderivative1;
    typedef Calculus::PrimalAntiderivative2       PrimalAntiderivative2;
    typedef Calculus::DualAntiderivative1         DualAntiderivative1;
    typedef Calculus::DualAntiderivative2         DualAntiderivative2;
    typedef Calculus::PrimalHodge0                PrimalHodge0;
    typedef Calculus::PrimalHodge1                PrimalHodge1;
    typedef Calculus::PrimalHodge2                PrimalHodge2;
    typedef Calculus::DualHodge0                  DualHodge0;
    typedef Calculus::DualHodge1                  DualHodge1;
    typedef Calculus::DualHodge2                  DualHodge2;
    typedef LinearAlgebra::SolverSimplicialLLT    LinearAlgebraSolver;
    typedef DiscreteExteriorCalculusSolver<Calculus, LinearAlgebraSolver,
                                           0, DUAL, 0, DUAL> SolverU;

    BOOST_STATIC_ASSERT(( KSpace::dimension == 3 ));

    /// Creates a parallelogram of vertices \a a, \a b and \a c, and
    /// last vertex is computed as \f$ a + b-a + c-a \f$.
    inline
    ImplicitDigitalVolume( Clone<BooleanImage> anImage )
      : calculus(), u( calculus ), bimage( anImage ), ptrSurface( 0 )
    {
      K.init( bimage.domain().lowerBound(),
              bimage.domain().upperBound(), true );
      SCellSet boundary;
      Surfaces<KSpace>::sMakeBoundary( boundary, K, bimage,
                                       K.lowerBound(), K.upperBound() );
      ptrSurface = new Surface( new SurfaceStorage( K, true, boundary ) );
      // Build trivial normal field.
      for ( auto surfel : boundary )
        {
          Vector N = trivialNormal( surfel );
          normals[ surfel ] = N;
        }
    }
    
    /// Virtual destructor since object contains virtual methods.
    ~ImplicitDigitalVolume() {}

    BooleanImage&  image()   { return bimage; }
    KSpace&        space()   { return K; }
    const Surface& surface() { return *ptrSurface; }

    /// computes the implicit function
    void computeImplicitFunction()
    {
      const Scalar alpha = 0.001;
      DGtal::trace.beginBlock( "Computing DEC" );
      DGtal::trace.info() << "- add cells" << std::endl;
      calculus.myKSpace = space();
      // Adds all the cell
      for ( auto surfel : surface() )
        {
          calculus.insertSCell( surfel );
          Dimension k = space().sOrthDir( surfel );
          calculus.insertSCell( space().sIncident( surfel, k, true ) );
          calculus.insertSCell( space().sIncident( surfel, k, false ) );
        }
      calculus.updateIndexes();
      DGtal::trace.info() << "- dual_D0" << std::endl;
      DualDerivative0 dual_D0  = calculus.template derivative<0,DUAL>();
      DualAntiderivative1 dual_AD1  = calculus.template antiderivative<1,DUAL>();
      DualIdentity0   dual_Id0 = calculus.template identity  <0,DUAL>();
      // DualIdentity0   M        = dual_D0.transpose() * dual_D0 + alpha * dual_Id0;
      DualIdentity0   M        = dual_AD1 * dual_D0 + alpha * dual_Id0;
      DualForm1       n( calculus );
      DualForm0       f( calculus );
      DGtal::trace.info() << "- dual 1-form n and dual 0-form f" << std::endl;
      for ( auto it = normals.cbegin(), itE = normals.cend(); it != itE; ++it )
        {
          SCell  surfel = it->first;
          Dimension   k = space().sOrthDir( surfel );
          SCell  in_vox = space().sDirectIncident( surfel, k );
          SCell out_vox = space().sIndirectIncident( surfel, k );
          Index in_idx  = calculus.getCellIndex( space().unsigns( in_vox ) );
          Index out_idx = calculus.getCellIndex( space().unsigns( out_vox ) );
          f.myContainer( in_idx )  = 1; // 0.5;
          f.myContainer( out_idx ) = 0; // -0.5;
        }
      DualForm1       t( calculus );
      t = dual_D0 * f;
      for ( auto it = normals.cbegin(), itE = normals.cend(); it != itE; ++it )
        {
          SCell  surfel = it->first;
          Dimension   k = space().sOrthDir( surfel );
          RealVector Ne = it->second;
          Index    idx  = calculus.getCellIndex( space().unsigns( surfel ) );
          std::cout << " " << idx << ":" << t.myContainer( idx );
          n.myContainer( idx )     = t.myContainer( idx ) > 0.0
            ? Ne[ k ] : -Ne[ k ];
        }
      DGtal::trace.info() << "- prefactoring matrix M := A'^t A' + a Id" << std::endl;
      SolverU solver;
      solver.compute( M );
      DGtal::trace.info() << "- solving M u = A'^t n + a f" << std::endl;
      // DualForm0 v = dual_D0.transpose() * n + alpha * f;
      DualForm0 v = dual_AD1 * n + alpha * f;
      u = solver.solve( v );
      DGtal::trace.info() << ( solver.isValid() ? "=> OK" : "ERROR" )
                          << " " << solver.myLinearAlgebraSolver.info() << std::endl;
      // TODO
      DGtal::trace.endBlock();
      }

    template <typename GrayLevelImage>
    void getImplicitFunctionImage( GrayLevelImage& output, int max_value = 255 )
    {
      output = GrayLevelImage( bimage.domain() );
      for ( auto p : output.domain() )
        output.setValue( p, bimage( p ) ? max_value : 0 );
      const Index dual_nb0 = u.myContainer.rows();
      for ( Index i = 0; i < dual_nb0; i++ )
        {
          Scalar val = u.myContainer( i );
          SCell  vox = u.getSCell( i );
          Point  p   = space().sCoords( vox );
          Scalar bdv = std::min( 1.0, std::max( 0.0, val ) ); //cut-off values
          int    glv = (int) round( bdv * max_value );
          // std::cout << " " << val;
          output.setValue( p, glv );
        }
    }

    /// @return the trivial normal to the given surfel.
    Vector trivialNormal( SCell surfel ) const
    {
      Dimension k  = K.sOrthDir( surfel );
      bool ext_dir = ! K.sDirect( surfel, k );
      Vector t; // zero vector
      t[ k ] = ext_dir ? 1.0 : -1.0;
      return t;
    }

    // ----------------------------------------------------------------------

    /// The discrete exterior calculus instance.
    Calculus calculus;
    /// The implicit function per voxel
    DualForm0 u;
    /// A reference to the boolean image.
    BooleanImage bimage;
    /// the Khalimsky space
    KSpace       K;
    /// The digital surface corresponding to the boundary of bimage.
    Surface*     ptrSurface;
    /// Map cell -> input normal vector
    VectorField  normals;
  };

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  using namespace DGtal;
  namespace po = boost::program_options;
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "input binary (or binarized) vol file (.vol) , pgm3d (.p3d or .pgm3d, pgm (with 3 dims)) file or sdp (sequence of discrete points)" )
    ("output,o", po::value<std::string>(), "output gray-level vol file (.vol) , pgm3d (.p3d or .pgm3d, pgm (with 3 dims)) file or sdp (sequence of discrete points)" )
    ("thresholdMin,m",  po::value<int>()->default_value(0), "threshold min (excluded) to define binary shape" )
    ("thresholdMax,M",  po::value<int>()->default_value(255), "threshold max (included) to define binary shape" )
    ("gridstep,g", po::value< double >()->default_value( 1.0 ), "the gridstep that defines the digitization (often called h). " )
    ("estimator,e", po::value<string>()->default_value( "VCM" ), "the chosen normal estimator: VCM | II " )
    ("R-radius,R", po::value<double>()->default_value( 5 ), "the constant for parameter R in R(h)=R h^alpha (VCM)." )
    ("r-radius,r", po::value<double>()->default_value( 3 ), "the constant for parameter r in r(h)=r h^alpha (VCM,II,Trivial)." )
    ("kernel,k", po::value<string>()->default_value( "hat" ), "the function chi_r, either hat or ball." )
    ("alpha", po::value<double>()->default_value( 0.0 ), "the parameter alpha in r(h)=r h^alpha (VCM)." )
    ("trivial-radius,t", po::value<double>()->default_value( 3 ), "the parameter t defining the radius for the Trivial estimator. Also used for reorienting the VCM." )
    ("embedding,E", po::value<int>()->default_value( 0 ), "the surfel -> point embedding for VCM estimator: 0: Pointels, 1: InnerSpel, 2: OuterSpel." )
    ;
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }catch(const exception& ex){
    parseOK=false;
    cerr << "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify(vm);
  if( ! parseOK || vm.count("help") || ! vm.count( "input" ) || ! vm.count( "output" ) )
    {
      cerr << "Usage: " << argv[0] << " -i <volume.vol> -o <volume.vol> -e <estimator> [options]\n"
           << "Uses the given digital normal estimator to smooth the binary volume such that its gradient is aligned with normal estimations. You may choose between several normal estimators (II|VCM), specified with -e."
           << endl
           << general_opt << "\n";
      cerr << "Example of use:\n"
           << "./volAntiAliasing -i \"fandisk-128.vol\" -o \"fandisk-gl.vol -e VCM -R 3 -r 3 -t 2 -E 0" << endl << endl;
        return 0;
    }
  
  string inputFilename = vm["input"].as<string>();
  string outputFilename= vm["output"].as<string>();
  int thresholdMin     = vm["thresholdMin"].as<int>();
  int thresholdMax     = vm["thresholdMax"].as<int>();
  string estimator     = vm["estimator"].as<string>();

  trace.beginBlock( "Reading vol file into an image." );
  typedef KhalimskySpaceND< 3, int >                    KSpace;
  typedef KSpace::Space                                 Space;
  typedef HyperRectDomain< Space >                      Domain;
  typedef ImageContainerBySTLVector< Domain, int >      Image;
  typedef ImageContainerBySTLVector< Domain, unsigned char > GrayLevelImage;
  typedef ImageContainerBySTLVector< Domain, bool >     BooleanImage;
  typedef functors::IntervalForegroundPredicate<Image>  ThresholdedImage;
  typedef ImplicitDigitalVolume< KSpace, BooleanImage > Volume;
  Image image = VolReader<Image>::importVol(inputFilename);
  ThresholdedImage thresholdedImage( image, thresholdMin, thresholdMax );
  trace.endBlock();
  trace.beginBlock( "Making binary image and building digital volume." );
  BooleanImage bimage( image.domain() );
  for ( auto p : bimage.domain() )
    bimage.setValue( p, thresholdedImage( p ) );
  // Build volume.
  Volume vol( bimage );
  trace.endBlock();
  // Setting normals
  rt::chooseKernel( vm, vol.space(), vol.surface(), vol.image(), vol.normals );
  // TODO ...
  trace.beginBlock( "Computing implicit function." );
  vol.computeImplicitFunction();
  trace.endBlock();
  trace.beginBlock( "Outputing vol image." );
  GrayLevelImage output( image.domain() );
  vol.getImplicitFunctionImage( output, 255 );
  output >> outputFilename;
  trace.endBlock();
  return 0;
}

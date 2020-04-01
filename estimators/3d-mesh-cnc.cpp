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
 * @file 3d-mesh-cnc.cpp
 * @ingroup surfaceTools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2020/02/14
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
#include "DGtal/math/Statistic.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/helpers/Shortcuts.h"
#include "DGtal/helpers/ShortcutsGeometry.h"
#include "DGtal/shapes/PolygonalSurface.h"
#include "EstimatorHelpers.h"
#include "SimplifiedMesh.h"
#include "CorrectedNormalCurrentComputer.h"
#include "RusinkiewiczCurvatureFormula.h"
#include "NormalCycleComputer.h"
//#include "FastCorrectedNormalCurrent.h"


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

namespace DGtal {
  template < typename TPoint = Z3i::RealPoint,
	     typename TVector =  Z3i::RealVector >
  struct PolygonalSurfaceReader {
    typedef TPoint Point;
    typedef TVector Vector;
    typedef PolygonalSurfaceReader< Point, Vector > Self;
    typedef DGtal::PolygonalSurface< Point > PolygonalSurface;
    typedef typename PolygonalSurface::template IndexedPropertyMap< Vector > NormalMap ;
    typedef typename PolygonalSurface::Index Index;

    static bool verifyIndices( const std::vector< Index > indices )
    {
      std::set<Index> sindices( indices.begin(), indices.end() );
      return sindices.size() == indices.size();
    }
    static
    std::vector< std::string > split( const std::string& str, char delim = ' ')
    {
      std::stringstream ss(str);
      std::string token;
      std::vector< std::string > cont;
      while ( std::getline( ss, token, delim ) ) cont.push_back(token);
      return cont;
    }    
    static
    bool read( std::istream& input,
	       PolygonalSurface& psurf,
	       NormalMap& vtx_normal_map  = NormalMap(),
	       NormalMap& face_normal_map = NormalMap() )
    {
      psurf.clear();
      std::vector<Point>  vertices;
      std::vector<Vector> normals;
      std::vector< std::vector< Index > > faces;
      std::vector< std::vector< Index > > faces_normals;
      std::string linestr;
      std::string keyword;
      std::string indices;
      Point  p;
      Vector n;
      std::getline( input, linestr );
      Index l = 0;
      for ( ; input.good() && ! input.eof(); std::getline( input, linestr ), l++ )
	{
	  if ( linestr.empty() ) continue; // skip empty line
	  if ( linestr[0] == '#' ) continue; // skip comment line
	  istringstream lineinput( linestr );
	  std::operator>>( lineinput,  keyword );
	  if ( keyword == "v" ) {
	    lineinput >> p[ 0 ] >> p[ 1 ] >> p[ 2 ];
	    // std::cout << "[" << l << "] v " << p << std::endl;
	    vertices.push_back( p );
	  } else if ( keyword == "vn" ) {
	    lineinput >> n[ 0 ] >> n[ 1 ] >> n[ 2 ];
	    normals.push_back( n );
	  } else if ( keyword == "f" ) {
	    std::vector< Index > face, face_normals;
	    while ( ! lineinput.eof() ) {
	      std::operator>>( lineinput, indices);
	      if ( indices.empty() ) break;
	      auto vtxinfo = split( indices, '/' );
	      if ( vtxinfo.size() == 0 ) break;
	      Index v  = std::stoi( vtxinfo[ 0 ] );
	      Index vn = vtxinfo.size() >= 3 ? std::stoi( vtxinfo[ 2 ] ) : v;
	      face.push_back( v - 1 );
	      face_normals.push_back( vn - 1 );
	      indices = "";
	    }
	    if ( ! face.empty() && verifyIndices( face ) ) {
	      faces.push_back( face );
	      faces_normals.push_back( face );
	    }
	  }
	  // Weird: necessary to clear them.
	  keyword = ""; linestr = "";
	} // while ( ! input.eof() )
      // Creating Polygonal Surface
      std::cout << "#V=" << vertices.size()
		<< " #F=" << faces.size() << std::endl;
      for ( auto v : vertices ) psurf.addVertex( v );
      for ( auto f : faces )    psurf.addPolygonalFace( f );
      bool ok = psurf.build();
      if ( ! ok ) return false;
      if ( ( ! normals.empty() ) && ( normals.size() == vertices.size() ) )
	{ // Build vertex normal map
	  vtx_normal_map = psurf.makeVertexMap( Vector() );
	  Index i = 0;
	  for ( auto n : normals ) vtx_normal_map[ i++ ] = n;
	}
      if ( ! normals.empty() )
	{ // Build face normal map
	  face_normal_map = psurf.makeFaceMap( Vector() );
	  Index i = 0;
	  for ( auto face_n_indices : faces_normals )
	    { 
	      Vector n;
	      for ( auto k : face_n_indices ) n += normals[ k ];
	      n /= face_n_indices.size();
	      face_normal_map[ i++ ] = n;
	    }
	}
      return ! input.bad();
    }
  };
} // namespace DGtal

  // PolySurf  psurf;
  // NormalMap vtx_normal_map;
  // NormalMap face_normal_map;
  // auto mesh_file = vm[ "input"   ].as<std::string>();
  // trace.info() << "Reading file <" << mesh_file << ">" << std::endl;
  // ifstream mesh_input( mesh_file.c_str() );
  // bool ok = PolySurfReader::read( mesh_input, psurf, vtx_normal_map, face_normal_map );
  // if ( ! ok ) {
  //   trace.error() << "Error reading file <" << mesh_file << ">" << std::endl;
  //   return 1;
  // }
  // trace.info() << "PolygonalSurface: #V=" << psurf.nbVertices()
  // 	       << " #E=" << psurf.nbEdges()
  // 	       << " #F=" << psurf.nbFaces()
  // 	       << " #NV=" << vtx_normal_map.size()
  // 	       << " #NF=" << face_normal_map.size()
  // 	       << std::endl;

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

// template <typename Colors>
// Colors getSimplifiedColorMap( const Colors& colors, unsigned char lost_mask = 0x03 )
// {
//   auto    new_colors = colors;
//   unsigned char mask = 0xff ^ lost_mask;
//   for ( Color& c : new_colors )
//     {
//       c.red  ( c.red()   & mask );
//       c.green( c.green() & mask );
//       c.blue ( c.blue()  & mask );
//       c.alpha( c.alpha() & mask );
//     }
//   return new_colors;
// }

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

/**
   Main class pour computing curvatures on mesh/digital surfaces/polynomials
*/
struct CurvatureComputer
{
  typedef Z3i::Space                            Space;
  typedef Z3i::KSpace                           KSpace;
  typedef Space::RealPoint                      RealPoint;
  typedef Space::RealVector                     RealVector;
  typedef PolygonalSurface<RealPoint>           PolySurf;
  typedef PolygonalSurfaceReader<RealPoint,RealVector> PolySurfReader;
  typedef PolySurfReader::NormalMap             NormalMap;
  typedef SimplifiedMesh<RealPoint,RealVector>  SimpleMesh;
  typedef SimplifiedMeshReader<RealPoint,RealVector> SimpleMeshReader;
  typedef SimplifiedMeshWriter<RealPoint,RealVector> SimpleMeshWriter;
  typedef SimplifiedMeshHelper<RealPoint,RealVector> SimpleMeshHelper;
  typedef SimpleMesh::Index                     Index;
  typedef Shortcuts<KSpace>                     SH;
  typedef ShortcutsGeometry<KSpace>             SHG;
  typedef EstimatorHelpers<KSpace>              EH;
  typedef SH::IdxDigitalSurface                 Surface;
  typedef Surface::DigitalSurfaceContainer      Container;
  typedef SH::BinaryImage                       BinaryImage;
  typedef SH::ImplicitShape3D                   ImplicitShape3D;
  typedef SH::DigitalSurface                    DigitalSurface;
  typedef SH::IdxDigitalSurface                 IdxDigitalSurface;
  typedef CorrectedNormalCurrentComputer< RealPoint, RealVector > CNCComputer;
  typedef NormalCycleComputer< RealPoint, RealVector > NCComputer;
  typedef CNCComputer::Scalar                   Scalar;
  typedef CNCComputer::Scalars                  Scalars;
  typedef CNCComputer::RealTensors              RealTensors;
  typedef std::vector< RealPoint >              RealPoints;
  typedef std::vector< RealVector >             RealVectors;
  typedef GradientColorMap<Scalar>              ColorMap;

  /// Creates the object with the meaningful options.
  CurvatureComputer();
  /// Parses command-line options.
  bool parseCommandLine( int argc, char** argv );

  /// Build input from data for further processing.
  bool buildInput();
  /// Build inputs from implicit polynomial shapes or 3D .vol image file.
  bool buildPolynomialOrVolInput();
  /// Build inputs from a set of predefined meshes
  bool buildPredefinedMesh();
  /// Build inputs from a polygonal mesh .obj file
  bool buildObjInput();

  /// Perturbates positions
  void perturbatePositions();
  /// Compute normals
  void computeNormals();

  /// Calls the appropriate curvature method.
  bool computeCurvatures();
  /// Computes interpolated CNC for estimating curvatures
  bool computeInterpolatedCorrectedNormalCurrent();
  /// Compute Rusinkiewicz curvatures
  bool computeRusinkiewiczCurvatures();
  /// Computes normal cycle for estimating curvatures
  bool computeNormalCycle();
  /// Compute all curvatures from curvature measures, i.e. localize
  /// measure computations.
  bool computeCurvaturesFromMeasures();
  /// If face or vertex curvatures are not computed, interpolate them
  /// from the others.
  bool computeVertexFaceCurvatures();

  /// Output curvatures (G, H, etc) as mesh OBJ files
  bool outputCurvaturesAsMeshObj();
  /// Output curvature errors (G, H, etc) as mesh OBJ files
  bool outputCurvatureErrorsAsMeshObj();
  /// Output Curvature errors statistics.
  bool outputCurvatureErrorStatistics();
  /// Dedicated method for outputing the statistics of any scalar
  /// curvature.
  bool outputScalarCurvatureErrorStatistics
  ( std::ostream& output, std::string information,
    Scalars expected_values, Scalars measured_values );

  ColorMap computeGeometryTypeColormap();
  Color computeGeometryColor( Scalar k1, Scalar k2,
			      Scalar zero, Scalar max ) const;
  
  //---------------------------------------------------------------------------
public:
  
  /// Meaningful options for this object.
  po::options_description general_opt;
  /// Parsed command-line
  po::variables_map vm;
  /// command-line as string
  std::string       command_line;
  /// Input filename
  std::string       filename;
  /// Mesh name
  std::string       meshname;
  /// Mesh args
  std::vector< std::string > meshargs;
  /// Possible quantities
  std::string       quantity;
  /// Possible anisotropy
  std::string       anisotropy;
  /// Gridstep
  double            h;
  /// Subdivision parameters for predefined meshes.
  double            sub_m, sub_n;
  /// Parameters object
  Parameters        params;
  /// The simplified mesh on which all computations are done
  SimpleMesh    smesh;
  /// ColorMap used for describing local geometry type (convexe, concave, etc).
  ColorMap      geom_cmap;
  
  bool in_obj;
  bool in_vol;
  bool in_mesh;
  bool in_polynomial;

  KSpace                        K;   // Digital space.
  CountedPtr<BinaryImage>       bimage;
  CountedPtr<ImplicitShape3D>   shape;
  CountedPtr<DigitalSurface>    surface;
  CountedPtr<IdxDigitalSurface> idx_surface;
  bool blocky;
  bool dual;
  
  /// Stores expected normal vectors (if applicable)
  SH::RealVectors   expected_normals;
  /// Stores expected mean curvature values (if applicable)
  SH::Scalars       expected_H_face_values;
  /// Stores expected mean curvature values (if applicable)
  SH::Scalars       expected_H_vertex_values;
  /// Stores expected Gaussian curvature values (if applicable)
  SH::Scalars       expected_G_face_values;
  /// Stores expected Gaussian curvature values (if applicable)
  SH::Scalars       expected_G_vertex_values;
  /// Stores expected first (smallest) principal curvature values (if applicable)
  SH::Scalars       expected_K1_face_values;
  /// Stores expected second (greatest) principal curvature values (if applicable)
  SH::Scalars       expected_K2_face_values;
  /// Stores expected first (smallest) principal direction values (if applicable)
  SH::RealVectors   expected_D1_face_values;
  /// Stores expected second (greatest) principal direction values (if applicable)
  SH::RealVectors   expected_D2_face_values;
  /// Stores expected first (smallest) principal curvature values (if applicable)
  SH::Scalars       expected_K1_vertex_values;
  /// Stores expected second (greatest) principal curvature values (if applicable)
  SH::Scalars       expected_K2_vertex_values;
  /// Stores expected first (smallest) principal direction values (if applicable)
  SH::RealVectors   expected_D1_vertex_values;
  /// Stores expected second (greatest) principal direction values (if applicable)
  SH::RealVectors   expected_D2_vertex_values;

  /// Stores measured mean curvature values (if applicable)
  SH::Scalars       measured_H_face_values;
  /// Stores measured mean curvature values (if applicable)
  SH::Scalars       measured_H_vertex_values;
  /// Stores measured Gaussian curvature values (if applicable)
  SH::Scalars       measured_G_face_values;
  /// Stores measured Gaussian curvature values (if applicable)
  SH::Scalars       measured_G_vertex_values;
  /// Stores measured normal vectors (if applicable)
  SH::RealVectors   measured_normals;
  /// Stores measured first (smallest) principal curvature values (if applicable)
  SH::Scalars       measured_K1_face_values;
  /// Stores measured second (greatest) principal curvature values (if applicable)
  SH::Scalars       measured_K2_face_values;
  /// Stores measured first (smallest) principal direction values (if applicable)
  SH::RealVectors   measured_D1_face_values;
  /// Stores measured second (greatest) principal direction values (if applicable)
  SH::RealVectors   measured_D2_face_values;
  /// Stores measured first (smallest) principal curvature values (if applicable)
  SH::Scalars       measured_K1_vertex_values;
  /// Stores measured second (greatest) principal curvature values (if applicable)
  SH::Scalars       measured_K2_vertex_values;
  /// Stores measured first (smallest) principal direction values (if applicable)
  SH::RealVectors   measured_D1_vertex_values;
  /// Stores measured second (greatest) principal direction values (if applicable)
  SH::RealVectors   measured_D2_vertex_values;

  /// Provides statistics on expected Gaussian curvatures (if applicable)
  Statistic<double> stat_expected_G_curv;
  /// Provides statistics on measured Gaussian curvatures (if applicable)
  Statistic<double> stat_measured_G_curv;
  /// Provides statistics on errors for Gaussian curvatures (if applicable)
  Statistic<double> stat_error_G_curv;
  /// Provides statistics on expected mean curvatures (if applicable)
  Statistic<double> stat_expected_H_curv;
  /// Provides statistics on measured mean curvatures (if applicable)
  Statistic<double> stat_measured_H_curv;
  /// Provides statistics on errors for mean curvatures (if applicable)
  Statistic<double> stat_error_H_curv;

  Scalars       measured_ball_mu0;
  Scalars       measured_ball_mu1;
  Scalars       measured_ball_mu2;
  RealTensors   measured_ball_muXY;
  
  double time_curv_ground_truth;
  double time_curv_estimations;
  double time_mu_estimations;
  double time_normal_estimations;
  
}; // struct CurvatureComputer


int main( int argc, char** argv )
{  

  CurvatureComputer CC;
  bool ok_parse = CC.parseCommandLine( argc, argv );
  if ( ! ok_parse ) return 0;

  CC.buildInput();
  CC.perturbatePositions();
  CC.computeNormals();
  CC.computeCurvatures();
  CC.computeVertexFaceCurvatures();
  CC.outputCurvaturesAsMeshObj();
  CC.outputCurvatureErrorsAsMeshObj();
  CC.outputCurvatureErrorStatistics();

  // const auto minValue         = vm[ "minValue" ].as<double>();
  // const auto maxValue         = vm[ "maxValue" ].as<double>();
  // const auto colormap_name    = vm[ "colormap" ].as<string>();
  // const auto zt               = vm[ "zero-tic" ].as<double>();
  // const double vthickness      = 0.1;
  // const double vlength         = 1.5;
  // auto outputfile    = vm[ "output"   ].as<string>();
  // auto estimator     = vm[ "estimator" ].as<string>();
  // const double h     = vm[ "gridstep"  ].as<double>();
  // const double mcoef = vm[ "m-coef"    ].as<double>();
  // const double mpow  = vm[ "m-pow"     ].as<double>();
  // const double mr    = mcoef * pow( h, mpow );
  // SH::Scalars       displayed_values;
  // SH::Scalars       expected_values;
  // SH::Scalars       measured_values;
  // SH::RealVectors   expected_normals;
  // SH::RealVectors   measured_normals;
  // Statistic<double> stat_expected_curv;
  // Statistic<double> stat_measured_curv;
  // double time_curv_ground_truth  = 0.0;
  // double time_curv_estimations   = 0.0;
  // double time_mu_estimations     = 0.0;
  // double time_normal_estimations = 0.0;
  
  // trace.beginBlock( "Compute surfels" );
  // params( "surfaceTraversal", "Default" );
  // const auto     surfels = SH::getSurfelRange( surface, params );
  // params( "surfaceTraversal", "DepthFirst" );
  // //params( "surfaceTraversal", "Default" ); //< Force default (bug in II ?)
  // const auto dft_surfels = SH::getSurfelRange( surface, params );
  // const auto       match = SH::getRangeMatch ( surfels, dft_surfels );
  // trace.endBlock();

  // // Compute true and/or estimated normal vectors
  // if ( vm.count( "polynomial" ) )
  //   {
  //     trace.beginBlock( "Compute true normals" );
  //     expected_normals = SHG::getNormalVectors( shape, K, dft_surfels, params );
  //     trace.endBlock();
  //   }
  // trace.beginBlock( "Estimate normals" );
  // if ( estimator == "True" && vm.count( "polynomial" ) )
  //   measured_normals = expected_normals;
  // else if ( estimator == "CTrivial" )
  //   measured_normals = SHG::getCTrivialNormalVectors( surface, dft_surfels, params );
  // else if ( estimator == "VCM" )
  //   measured_normals = SHG::getVCMNormalVectors( surface, dft_surfels, params );
  // else if ( estimator == "II" ) {
  //   measured_normals = SHG::getIINormalVectors( bimage, dft_surfels, params );
  //   auto oriented_normals = SHG::getCTrivialNormalVectors( surface, dft_surfels, params );
  //   SHG::orientVectors( measured_normals, oriented_normals );
  // }
  // else // default is "Trivial"
  //   measured_normals = SHG::getTrivialNormalVectors( K, dft_surfels );
  // time_normal_estimations = trace.endBlock();

  // // Compute true or estimated curvatures
  // if ( vm.count( "polynomial" ) )
  //   {
  //     trace.beginBlock( "Compute true curvatures" );
  //     expected_values = ( ( quantity == "H" ) || ( quantity == "Mu1" )
  // 			  || ( quantity == "HII" ) )
  // 	? SHG::getMeanCurvatures    ( shape, K, dft_surfels, params )
  // 	: SHG::getGaussianCurvatures( shape, K, dft_surfels, params );
  //     stat_expected_curv.addValues( expected_values.cbegin(), expected_values.cend() );
  //     stat_expected_curv.terminate();
  //     trace.info() << "- truth curv: avg = " << stat_expected_curv.mean() << std::endl;
  //     trace.info() << "- truth curv: min = " << stat_expected_curv.min() << std::endl;
  //     trace.info() << "- truth curv: max = " << stat_expected_curv.max() << std::endl;
  //     time_curv_ground_truth = trace.endBlock();
  //   }
  // if ( ( quantity == "HII" ) || ( quantity == "GII" ) )
  //   {
  //     trace.beginBlock( "Compute II curvature estimations" );
  //     measured_values = ( quantity == "HII" )
  // 	? SHG::getIIMeanCurvatures    ( bimage, dft_surfels, params )
  // 	: SHG::getIIGaussianCurvatures( bimage, dft_surfels, params );
  //     stat_measured_curv.addValues( measured_values.cbegin(), measured_values.cend() );
  //     stat_measured_curv.terminate();
  //     trace.info() << "- II curv: avg = " << stat_measured_curv.mean() << std::endl;
  //     trace.info() << "- II curv: min = " << stat_measured_curv.min() << std::endl;
  //     trace.info() << "- II curv: max = " << stat_measured_curv.max() << std::endl;
  //     time_curv_estimations = trace.endBlock();
  //   }
  // else 
  //   {
  //     trace.beginBlock( "Compute corrected normal current" );
  //     Current C( *idx_surface, h, vm.count( "crisp" ) );
  //     C.setCorrectedNormals( dft_surfels.begin(), dft_surfels.end(), measured_normals.begin() );
  //     C.setInterpolationMode( vm.count( "interpolate" ) );
  //     trace.info() << C << " m-ball-r = " << mr << "(continuous)"
  // 		   << " " << (mr/h) << " (discrete)" << std::endl;
  //     trace.info() << "Interpolation mode is "
  // 		   << ( vm.count( "interpolate" ) ? "ON" : "OFF" ) << std::endl;
  //     double              area = 0.0;
  //     double              intG = 0.0;
  //     std::vector<double> mu0( dft_surfels.size() );
  //     std::vector<double> mu1( dft_surfels.size() );
  //     std::vector<double> mu2( dft_surfels.size() );
  //     std::vector<double> muOmega( dft_surfels.size() );
  //     SH::RealVectors     muDir0( dft_surfels.size() );
  //     SH::RealVectors     muDir1( dft_surfels.size() );
  //     SH::RealVectors     muDir2( dft_surfels.size() );
  //     SH::Scalars         kappa0( dft_surfels.size() );
  //     SH::Scalars         kappa1( dft_surfels.size() );
  //     SH::Scalars         kappa2( dft_surfels.size() );
  //     bool       mu0_needed = true;
  //     bool       mu1_needed = false;
  //     bool       mu2_needed = false;
  //     bool   muOmega_needed = false;
  //     bool      muXY_needed = false;
  //     if ( quantity == "Mu0" )     mu0_needed = true;
  //     if ( quantity == "Mu1" )     mu1_needed = true;
  //     if ( quantity == "Mu2" )     mu2_needed = true;
  //     if ( quantity == "MuOmega" ) muOmega_needed = true;
  //     if ( quantity == "H" )       mu0_needed = mu1_needed = true;
  //     if ( quantity == "G" )       mu0_needed = mu2_needed = true;
  //     if ( quantity == "Omega" )   mu0_needed = muOmega_needed = true;
  //     if ( quantity == "MuXY" )    muXY_needed = true;
  //     if ( anisotropy == "NMult" || anisotropy == "NAdd" ) mu0_needed = true;
  //     trace.beginBlock( "Compute mu_k everywhere" );
  //     trace.info() << "computeAllMu0" << std::endl;
  //     if ( mu0_needed ) C.computeAllMu0();
  //     trace.info() << "computeAllMu1" << std::endl;
  //     if ( mu1_needed ) C.computeAllMu1();
  //     trace.info() << "computeAllMu2" << std::endl;
  //     if ( mu2_needed ) C.computeAllMu2();
  //     trace.info() << "computeAllMuOmega" << std::endl;
  //     if ( muOmega_needed ) C.computeAllMuOmega();
  //     trace.info() << "computeAllMuXYAnisotropic" << std::endl;
  //     if ( muXY_needed ) C.computeAllAnisotropicMu();
  //     time_mu_estimations = trace.endBlock();
  //     //#pragma omp parallel for schedule(dynamic)
  //     trace.info() << "compute measures" << std::endl;
  //     Vertex              i = 0;
  //     Vertex              j = dft_surfels.size();
  //     const bool normalized = ( anisotropy == "NMult" ) || ( anisotropy == "NAdd" );
  //     const bool M_times_Mt = ( anisotropy == "Mult" )  || ( anisotropy == "NMult" );
  //     for ( auto aSurfel : dft_surfels )
  // 	{
  // 	  // std::cout << i << " / " << j << std::endl;
  // 	  trace.progressBar( i, j );
  // 	  Vertex v = C.getVertex( aSurfel );
  // 	  area    += C.mu0( v );
  // 	  if ( mu0_needed )     mu0[ i ]     = C.mu0Ball( v, mr );
  // 	  if ( mu1_needed )     mu1[ i ]     = C.mu1Ball( v, mr );
  // 	  if ( mu2_needed )     mu2[ i ]     = C.mu2Ball( v, mr );
  // 	  if ( muOmega_needed ) muOmega[ i ] = C.muOmegaBall( v, mr );
  // 	  if ( muXY_needed && ! M_times_Mt ) {
  // 	    auto M = C.anisotropicMuBall( v, mr );
  // 	    auto N = measured_normals[ i ];
  // 	    M += M.transpose();
  // 	    M *= 0.5;
  // 	    const double   coef_N = normalized ? 100.0 * mu0[ i ] : 100.0;
  // 	    for ( int j = 0; j < 3; j++ )
  // 	       for ( int k = 0; k < 3; k++ )
  // 		 M( j, k ) += coef_N * N[ j ] * N[ k ];
  // 	    auto V = M;
  // 	    RealVector L;
  // 	    EigenDecomposition< 3, double>::getEigenDecomposition( M, V, L );
  // 	    muDir0[ i ] = V.column( 0 );
  // 	    muDir1[ i ] = V.column( 1 );
  // 	    muDir2[ i ] = V.column( 2 );
  //           kappa0[ i ] = normalized ? -L[ 0 ] / mu0[ i ] : -L[ 0 ];
  //           kappa1[ i ] = normalized ? -L[ 1 ] / mu0[ i ] : -L[ 1 ];
  //           kappa2[ i ] = normalized ? -L[ 2 ] / mu0[ i ] : -L[ 2 ];
  // 	  } else if ( muXY_needed && M_times_Mt ) {
  // 	    auto M = C.anisotropicMuBall( v, mr );
  // 	    auto MMt = M * M.transpose();
  // 	    auto V = M;
  // 	    RealVector L;
  // 	    EigenDecomposition< 3, double>::getEigenDecomposition( MMt, V, L );
  // 	    muDir0[ i ] = V.column( 0 );
  // 	    muDir1[ i ] = V.column( 1 );
  // 	    muDir2[ i ] = V.column( 2 );
  //           kappa0[ i ] = sqrt( fabs( L[ 0 ] ) ) / ( normalized ? mu0[ i ] : 1.0 );
  //           kappa1[ i ] = sqrt( fabs( L[ 1 ] ) ) / ( normalized ? mu0[ i ] : 1.0 );
  //           kappa2[ i ] = sqrt( fabs( L[ 2 ] ) ) / ( normalized ? mu0[ i ] : 1.0 );
  // 	  }
  // 	  ++i;
  // 	}
  //     // Computing total Gauss curvature.
  //     if ( mu2_needed )
  // 	{
  // 	  trace.info() << "compute total Gauss curvature" << std::endl;
  // 	  for ( auto f : idx_surface->allFaces() ) intG += C.mu2( f );
  // 	}
  //     if ( quantity == "Mu0" )          measured_values = mu0;
  //     else if ( quantity == "Mu1" )     measured_values = mu1;
  //     else if ( quantity == "Mu2" )     measured_values = mu2;
  //     else if ( quantity == "MuOmega" ) measured_values = muOmega;
  //     else if ( quantity == "H" )
  // 	{
  // 	  measured_values.resize( dft_surfels.size() );
  // 	  std::transform( mu0.cbegin(), mu0.cend(),
  // 			  mu1.cbegin(), measured_values.begin(),
  // 			  [] ( double m0, double m1 ) { return m1 / (2.0*m0); } );
  // 	}
  //     else if ( quantity == "G" )
  // 	{
  // 	  measured_values.resize( dft_surfels.size() );
  // 	  std::transform( mu0.cbegin(), mu0.cend(),
  // 			  mu2.cbegin(), measured_values.begin(),
  // 			  [] ( double m0, double m2 ) { return m2 / m0; } );
  // 	}
  //     else if ( quantity == "Omega" )
  // 	{
  // 	  measured_values.resize( dft_surfels.size() );
  // 	  std::transform( mu0.cbegin(), mu0.cend(),
  // 			  muOmega.cbegin(), measured_values.begin(),
  // 			  [] ( double m0, double m2 ) { return m2 / sqrt( m0 ); } );
  // 	}
  //     if ( quantity != "MuXY" )
  // 	{
  // 	  for ( i = 0; i < j; ++i ) stat_measured_curv.addValue( measured_values[ i ] );
  // 	  stat_measured_curv.terminate();
  // 	  trace.info() << "- CNC area      = " << area << std::endl;
  // 	  trace.info() << "- CNC total G   = " << intG << std::endl;
  // 	  trace.info() << "- CNC mu: avg = " << stat_measured_curv.mean() << std::endl;
  // 	  trace.info() << "- CNC mu: min = " << stat_measured_curv.min() << std::endl;
  // 	  trace.info() << "- CNC mu: max = " << stat_measured_curv.max() << std::endl;
  // 	}
  //     time_curv_estimations = trace.endBlock();
  //     if ( quantity == "MuXY" ) {
  //       trace.beginBlock( "Output principal curvatures and directions files." );
  // 	SH::RealPoints positions( dft_surfels.size() );
  // 	std::transform( dft_surfels.cbegin(), dft_surfels.cend(), positions.begin(),
  // 			[&] (const SH::SCell& c) { return sembedder( c ); } );
  // 	auto pos0 = positions,   pos1 = positions,   pos2 = positions;
  //       auto max0 = kappa0[ 0 ], max1 = kappa1[ 0 ], max2 = kappa2[ 0 ];
  //       auto min0 = kappa0[ 0 ], min1 = kappa1[ 0 ], min2 = kappa2[ 0 ];
  // 	for ( SH::Idx i = 0; i < pos0.size(); ++i ) {
  // 	  muDir0[ i ] *= vlength;
  // 	  muDir1[ i ] *= vlength;
  // 	  muDir2[ i ] *= vlength;
  // 	  pos0[ i ] -= 0.5 * muDir0[ i ];
  // 	  pos1[ i ] -= 0.5 * muDir1[ i ];
  // 	  pos2[ i ] -= 0.5 * muDir2[ i ];
  //         max0 = std::max( max0, kappa0[ i ] );
  //         max1 = std::max( max1, kappa1[ i ] );
  //         max2 = std::max( max2, kappa2[ i ] );
  //         min0 = std::min( min0, kappa0[ i ] );
  //         min1 = std::min( min1, kappa1[ i ] );
  //         min2 = std::min( min2, kappa2[ i ] );
  // 	}
  //       trace.info() << min0 << " <= kappa0 <= " << max0 << std::endl;
  //       trace.info() << min1 << " <= kappa1 <= " << max1 << std::endl;
  //       trace.info() << min2 << " <= kappa2 <= " << max2 << std::endl;
  // 	SH::saveOBJ( surface, SH::RealVectors(), SH::Colors(),
  // 		     outputfile+"-primal.obj" );
  //       auto colors = SH::Colors( dft_surfels.size() );
  // 	// TODO 
  // 	const auto cmap_kappa
  // 	  = SH::getZeroTickedColorMap( minValue, maxValue,
  // 				       params( "colormap", colormap_name )
  // 				       ( "zero-tic", zt ) );
  //       for ( SH::Idx i = 0; i < colors.size(); i++ )
  // 	  colors[ i ] = cmap_kappa( kappa0[ i ] );
  // 	SH::saveVectorFieldOBJ( pos0, muDir0, vthickness, getSimplifiedColorMap( colors ), outputfile+"-dir0.obj" );
  //       for ( SH::Idx i = 0; i < colors.size(); i++ )
  // 	  colors[ i ] = cmap_kappa( kappa1[ i ] );
  // 	SH::saveVectorFieldOBJ( pos1, muDir1, vthickness, getSimplifiedColorMap( colors ), outputfile+"-dir1.obj" );
  // 	// normal direction, color has no meaning
  // 	SH::saveVectorFieldOBJ( pos2, muDir2, vthickness, SH::Colors(),
  // 		     outputfile+"-dir2.obj", SH::Color::Black, SH::Color::Blue );
  // 	const auto colormap0 = SH::getColorMap( min0, max0,
  // 						params( "colormap", colormap_name ) );
  //       for ( SH::Idx i = 0; i < colors.size(); i++ )
  // 	  colors[ i ] = colormap0( kappa0[ match[ i ] ] );
  //       SH::saveOBJ( surface, SH::RealVectors(), getSimplifiedColorMap( colors ), outputfile+"-kappa0.obj" );
  //       const auto colormap1 = SH::getColorMap( min1, max1,
  // 						params( "colormap", colormap_name ) );
  //       for ( SH::Idx i = 0; i < colors.size(); i++ )
  // 	  colors[ i ] = colormap1( kappa1[ match[ i ] ] );
  //       SH::saveOBJ( surface, SH::RealVectors(), getSimplifiedColorMap( colors ), outputfile+"-kappa1.obj" );
  //       const auto colormap2 = SH::getColorMap( min2, max2,
  // 						params( "colormap", colormap_name ) );
  //       for ( SH::Idx i = 0; i < colors.size(); i++ )
  // 	  colors[ i ] = colormap2( kappa2[ match[ i ] ] );
  //       SH::saveOBJ( surface, SH::RealVectors(), getSimplifiedColorMap( colors ), outputfile+"-kappa2.obj" );

  // 	SH::saveOBJ( surface, SH::RealVectors(), SH::Colors(),
  // 		     outputfile+"-primal.obj" );
  // 	if ( vm.count( "polynomial" ) ) {
  // 	  SH::Cell2Index c2i;
  // 	  auto pointels = SH::getPointelRange( c2i, surface );
  // 	  SH::RealPoints pos( pointels.size() );
  // 	  std::transform( pointels.cbegin(), pointels.cend(), pos.begin(),
  // 			  [&] (const SH::Cell& c) { return h * embedder( c ); } ); 
  // 	  auto ppos     = SHG::getPositions( shape, pos, params );
  // 	  std::transform( ppos.cbegin(), ppos.cend(), pos.begin(),
  // 			  [&] (const SH::RealPoint& p) { return p / h; } );
  // 	  SH::saveOBJ( surface, [&] (const SH::Cell& c){ return pos[ c2i[ c ] ];},
  // 		       SH::RealVectors(), SH::Colors(),
  // 		       outputfile+"-qproj.obj" );
  // 	  SH::RealPoints pos_surf( dft_surfels.size() );
  // 	  std::transform( dft_surfels.cbegin(), dft_surfels.cend(), pos_surf.begin(),
  // 			  [&] (const SH::SCell& c) { return h * sembedder( c ); } ); 
  // 	  auto ppos_surf = SHG::getPositions( shape, pos_surf, params );
  // 	  std::transform( ppos_surf.cbegin(), ppos_surf.cend(), pos_surf.begin(),
  // 			  [&] (const SH::RealPoint& p) { return p / h; } );
  // 	  // auto real_pos = SH::getMatchedRange( pos_surf, match );
  // 	  pos0 = pos1 = pos2 = pos_surf; //real_pos;
  // 	  for ( SH::Idx i = 0; i < pos0.size(); ++i ) {
  // 	    pos0[ i ] -= 0.5 * muDir0[ i ];
  // 	    pos1[ i ] -= 0.5 * muDir1[ i ];
  // 	    pos2[ i ] -= 0.5 * muDir2[ i ];
  // 	  }
  //       for ( SH::Idx i = 0; i < colors.size(); i++ )
  // 	  colors[ i ] = cmap_kappa( kappa0[ i ] );
  // 	SH::saveVectorFieldOBJ( pos0, muDir0, vthickness, colors, outputfile+"-qproj-dir0.obj" );
  //       for ( SH::Idx i = 0; i < colors.size(); i++ )
  // 	  colors[ i ] = cmap_kappa( kappa1[ i ] );
  // 	SH::saveVectorFieldOBJ( pos1, muDir1, vthickness, colors, outputfile+"-qproj-dir1.obj" );
  // 	SH::saveVectorFieldOBJ( pos2, muDir2, vthickness, SH::Colors(),
  // 		     outputfile+"-qproj-dir2.obj", SH::Color::Black,
  // 		     SH::Color::Blue );
  // 	}
  // 	trace.endBlock();
  //     }
  //   }
  
  // if ( outputfile != "none" )
  //   {
  //     trace.beginBlock( "Save results as OBJ" );
  //     const bool has_ground_truth = expected_values.size() != 0;
  //     const bool has_estimations  = measured_values.size() != 0;
  //     const auto colormap         = SH::getZeroTickedColorMap( minValue, maxValue,
  //                                                              params( "colormap", colormap_name )
  //                                                              ( "zero-tic", zt ) );
      
  //     trace.info() << "#mvalues=" << measured_values.size() << std::endl;
  //     trace.info() << "#evalues=" << expected_values.size() << std::endl;

      
  //     auto colors = SH::Colors( surfels.size() );
  //     if ( has_ground_truth )
  // 	{
  // 	  for ( SH::Idx i = 0; i < colors.size(); i++ )
  // 	    colors[ i ] = colormap( expected_values[ match[ i ] ] ); 
  // 	  SH::saveOBJ( surface, SH::getMatchedRange( expected_normals, match ), colors,
  // 		       outputfile+"-truth.obj" );
  // 	}
  //     if ( has_estimations )
  // 	{
  // 	  for ( SH::Idx i = 0; i < colors.size(); i++ )
  // 	    colors[ i ] = colormap( measured_values[ match[ i ] ] ); 
  // 	  SH::saveOBJ( surface, SH::getMatchedRange( measured_normals, match ), colors,
  // 		       outputfile+"-estimation.obj" );
  // 	}
  //     if ( has_ground_truth && has_estimations )
  // 	{
  // 	  const auto error_values = SHG::getScalarsAbsoluteDifference( measured_values, expected_values );
  // 	  const auto max_error    = vm[ "max-error" ].as<double>();
  // 	  const auto stat_error   = SHG::getStatistic( error_values );
  // 	  const auto error_cmap   = getErrorColorMap( max_error );
  // 	  for ( SH::Idx i = 0; i < colors.size(); i++ )
  // 	    colors[ i ] = error_cmap( error_values[ match[ i ] ] ); 
  // 	  SH::saveOBJ( surface, SH::getMatchedRange( measured_normals, match ), colors,
  // 		       outputfile+"-error.obj" );
  // 	}


  //     trace.endBlock();
  //   }
  
  
  // if ( ! vm.count( "polynomial" ) ) return 0;
  
  // trace.beginBlock( "Output statistics" );
  // auto error_fname = vm[ "error" ].as<std::string>();
  // std::ofstream ferr;
  // ferr.open ( error_fname.c_str(),
  // 	      std::ofstream::out | std::ofstream::app );
  // if ( ! ferr.good() )
  //   trace.warning() << "Unable to open file " << error_fname << std::endl;
  // ferr << "######################################################################"
  //      << std::endl;
  // ferr << "# ";
  // for ( int i = 0; i < argc; ++i ) ferr << " " << argv[ i ];
  // ferr << std::endl;
  // ferr << "#---------------------------------------------------------------------"
  //      << std::endl;
  // ferr << "# Q=" << quantity
  //      << " P="  << vm[ "polynomial" ].as<std::string>()
  //      << " mr= " << mr << "(continuous)"
  //      << " mrd=" << (mr/h) << " (discrete)" << std::endl;
  // ferr << "# time_curv_ground_truth  = " << time_curv_ground_truth << " ms" << std::endl
  //      << "# time_normal_estimations = " << time_normal_estimations << " ms" << std::endl
  //      << "# time_curv_estimations   = " << time_curv_estimations << " ms" << std::endl
  //      << "# time_mu_estimations     = " << time_mu_estimations << " ms" << std::endl;
  // ferr << "# h(1) size(2) l1(3) l2(4) loo(5)"
  //      << " m_mean(6) m_dev(7) m_min(8) m_max(9)"
  //      << " exp_mean(10) exp_dev(11) exp_min(12) exp_max(13) " << std::endl;
  // ferr << "# mr(14) mrd(15) t_curv_gt(16) t_normal_est(17) t_curv_est(18) t_mu_est(19)"
  //      << " r-radius(20) m-coef(21)"
  //      << std::endl;
  // ferr << h << " " << measured_values.size()
  //      << " " << SHG::getScalarsNormL1 ( measured_values, expected_values )
  //      << " " << SHG::getScalarsNormL2 ( measured_values, expected_values )
  //      << " " << SHG::getScalarsNormLoo( measured_values, expected_values )
  //      << " " << stat_measured_curv.mean()
  //      << " " << sqrt( stat_measured_curv.variance() )
  //      << " " << stat_measured_curv.min()
  //      << " " << stat_measured_curv.max()
  //      << " " << stat_expected_curv.mean()
  //      << " " << sqrt( stat_expected_curv.variance() )
  //      << " " << stat_expected_curv.min()
  //      << " " << stat_expected_curv.max();
  // ferr << " " << mr << " " << (mr/h) << " " << time_curv_ground_truth
  //      << " " << time_normal_estimations << " " << time_curv_estimations
  //      << " " << time_mu_estimations
  //      << " " << vm[ "r-radius" ].as<double>()
  //      << " " << vm[ "m-coef" ].as<double>()
  //      << std::endl;
  // ferr << "#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
  //      << std::endl;
  // ferr.close();
  // trace.endBlock();
  
  return 0;
}

///////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
/// Creates the object with the meaningful options.
/////////////////////////////////////////////////////////////////////////////
CurvatureComputer::CurvatureComputer()
  : general_opt( "Allowed options are" ),
    geom_cmap( -1.0, 1.0 ),
    bimage( nullptr ),
    shape ( nullptr ),
    surface( nullptr ),
    idx_surface( nullptr )
{
  // initializes some parameters
  params = SH::defaultParameters() | SHG::defaultParameters();
  time_curv_ground_truth  = 0.0;
  time_curv_estimations   = 0.0;
  time_mu_estimations     = 0.0;
  time_normal_estimations = 0.0;

  // Create command line options
  general_opt.add_options()
    ( "help,h", "display this message" )
    ( "input,i", po::value<std::string>(), "input file: may be a mesh (.OBJ) or a volume image (.vol)" )
    ( "mesh", po::value<std::string>(), "input predefined mesh: {sphere[|-VN|-FN],r,m,n|lantern[|-VN|-FN],r,h,m,n|torus[-VN|-FN],R,r,m,n[,twist]}, where m/n is the number of latitudes/longitudes" )
    ( "output,o", po::value<std::string>()->default_value( "cnc" ), "the basename for output obj files or <none> if no output obj is wanted." )
    ( "average-normals,K", po::value<int>()->default_value( 0 ), "averages normals by performing <n> times vertexNormals -> faceNormals -> vertexNormals." )
    ( "unit-normals,u", "forces the interpolated normals to have unit norm." )
    ( "weights-normals-f2v", po::value<std::string>()->default_value( "EQUAL" ), "specifies how to average face normals when estimating vertex normals: EQUAL|MAX1995" )
    ( "m-coef", po::value<double>()->default_value( 3.0 ), "the coefficient k that defines the radius of the ball used in measures, that is r := k h^b" )
    ( "m-pow", po::value<double>()->default_value( 0.5 ), "the coefficient b that defines the radius of the ball used in measures, that is r := k h^b" )
    ( "uniform-noise", po::value<double>()->default_value( 0.0 ), "perturbates positions with a uniform random noise as a ratio r of average edge length." )
    ( "adaptive-noise", po::value<double>()->default_value( 0.0 ), "perturbates positions with a uniform random noise as a ratio r of local average edge length." );
  EH::optionsImplicitShape   ( general_opt );
  EH::optionsDigitizedShape  ( general_opt );
  general_opt.add_options()
    ("thresholdMin,m",  po::value<int>()->default_value(0), "threshold min (excluded) to define binary shape" )
    ("thresholdMax,M",  po::value<int>()->default_value(255), "threshold max (included) to define binary shape" )
    ("closed",  po::value<int>()->default_value(1), "tells if the cellular space is closed (1:default) or open." );
  //  EH::optionsVolFile         ( general_opt );
  EH::optionsNoisyImage      ( general_opt );
  EH::optionsNormalEstimators( general_opt );
  general_opt.add_options()
    ( "quantity,Q", po::value<std::string>()->default_value( "CNC" ), "the method that computes curvature quantities: CNC|RZ|NC" )
    ( "anisotropy", po::value<std::string>()->default_value( "NAdd" ), "tells how is symmetrized the anisotropic measure mu_XY, in Mult|Add|NMult|NAdd: Mult forces symmetry by M*M^t, Add forces symmetry by 0.5*(M+M^t)+NxN, NMult and NAdd normalized by the area.");
    // ( "quantity,Q", po::value<std::string>()->default_value( "H" ), "the quantity that is evaluated in Mu0|Mu1|Mu2|MuOmega|H|G|Omega|MuXY|HII|GII, with H := Mu1/(2Mu0), G := Mu2/Mu0, Omega := MuOmega/sqrt(Mu0), MuXY is the anisotropic curvature tensor and HII and GII are the mean and gaussian curvatures estimated by II." )
  //   ( "crisp,C", "when specified, when computing measures in a ball, do not approximate the relative intersection of cells with the ball but only consider if the cell centroid is in the ball (faster by 30%, but less accurate)." )
  //   ( "interpolate,I", "when specified, it interpolate the given corrected normal vector field and uses the corresponding measures." );
  EH::optionsDisplayValues   ( general_opt );
  general_opt.add_options()
    ( "zero-tic", po::value<double>()->default_value( 0.0 ), "adds a black band around zero of given thickness in colormaps." );
    
  // //#endif
  general_opt.add_options()
    ( "error", po::value<std::string>()->default_value( "error" ), "the name of the output file that sum up l2 and loo errors in estimation." )
    ( "max-error", po::value<double>()->default_value( 0.2 ), "the error value corresponding to black." );
  // general_opt.add_options()
  //   ( "output,o", po::value<std::string>()->default_value( "none" ), "the basename for output obj files or none if no output obj is wanted." );
  general_opt.add_options()
    ( "digital-surface", po::value<std::string>()->default_value( "DUAL" ), "chooses which kind of digital surface is used for computations in DUAL|PRIMAL|PDUAL|PPRIMAL: DUAL dual marching-cubes surface, PRIMAL blocky quad primal surface, PDUAL same as DUAL but projected onto polynomial true surface (if possible), PPRIMAL same as PRIMAL but projected onto polynomial true surface (if possible).");

  geom_cmap = computeGeometryTypeColormap();
} 

/////////////////////////////////////////////////////////////////////////////
/// Parse command line
/// @param argc the number of parameters
/// @param argv the array of parameters as C-strings.
/////////////////////////////////////////////////////////////////////////////
bool
CurvatureComputer::parseCommandLine( int argc, char** argv )
{
  in_obj        = false;
  in_vol        = false;
  in_mesh       = false;
  in_polynomial = false;

  bool parseOK = EH::args2vm( general_opt, argc, argv, vm );
  bool neededArgsGiven=true;

  if ( vm.count( "polynomial-list" ) )
    {
      trace.info() << "List of predefined polynomials:" << std::endl;
      auto L = SH::getPolynomialList();
      for ( auto p : L ) {
	trace.info() << "  " << p.first << " -> " << p.second << std::endl;
      }
      return false;
    }

  auto mesh     = vm.count( "mesh"  ) ? vm[ "mesh"  ].as<string>() : "";
  filename      = vm.count( "input" ) ? vm[ "input" ].as<string>() : "";
  auto idx      = filename.rfind('.');
  std::string extension = (idx != std::string::npos) ? filename.substr(idx+1) : "";
  trace.info() << "filename=" << filename << " extension=" << extension << std::endl;
  if ( parseOK && ( ! vm.count("polynomial") ) && ( ! vm.count( "input" ) )
       && ( ! vm.count( "mesh" ) ) ) {
    missingParam("--polynomial or --input or --mesh");
    neededArgsGiven = false;
  }
  if ( vm.count( "input" ) && ( extension != "vol" )
       && ( extension != "obj" ) && ( extension != "OBJ" ) ) {
    missingParam("Wrong input file extension (should be .vol, .obj, .OBJ)");
    neededArgsGiven = false;
  }

  meshargs = SimpleMeshReader::split( mesh, ',' );
  meshname = mesh != "" ? meshargs[ 0 ] : "";
  if ( meshname != ""
       && meshname != "sphere" && meshname != "sphere-VN" && meshname != "sphere-FN"
       && meshname != "lantern" && meshname != "lantern-VN" && meshname != "lantern-FN"
       && meshname != "torus" && meshname != "torus-VN" && meshname != "torus-FN" )
    {
      missingParam("Wrong predefined mesh (should be {sphere[|-VN|-FN]|lantern[|-VN|-FN]|torus[-VN|-FN])");
      neededArgsGiven = false;
    }
  
  if ( !neededArgsGiven || !parseOK || vm.count("help") || argc <= 1 )
    {
      trace.info()<< "Builds the 3d corrected normal currents" <<std::endl
                  << general_opt << "\n"
                  << "Basic usage: "<<std::endl
		  << "\t 3d-mesh-fcnc -i mesh.obj" << std::endl
		  << "\t 3d-mesh-fcnc -i image.vol" << std::endl
		  << "\t 3d-mesh-fcnc --mesh sphere,2.5,5,8" << std::endl
		  << "\t 3d-mesh-fcnc -p goursat" << std::endl
		  << "\t 3d-mesh-fcnc -p \"5*x^2-3*x*y*z+2*y^2-z^2-4*x*y^2\""
		  << std::endl << std::endl;
      return false;
    }
  
  command_line = std::string( argv[ 0 ] );
  for ( int i = 1; i < argc; ++i )
    command_line += std::string( " " ) + std::string( argv[ i ] );
  
  /////////////////////////////////////////////////////////////////////////////
  // Checking quantities
  /////////////////////////////////////////////////////////////////////////////
  h          = vm[ "gridstep" ].as<double>();
  quantity   = vm[ "quantity"   ].as<std::string>();
  anisotropy = vm[ "anisotropy" ].as<std::string>();
  std::vector< std::string > quantities = { "CNC", "RZ", "NC" };
    //{ "Mu0", "Mu1", "Mu2", "MuOmega", "H", "G", "Omega", "MuXY", "HII", "GII" };
  if ( std::count( quantities.begin(), quantities.end(), quantity ) == 0 ) {
    trace.error() << "Quantity should be in CNC|RZ";
    // trace.error() << "Quantity should be in Mu0|Mu1|Mu2|MuOmega|H|G|Omega|MuXY|HII|GII.";
    trace.info()  << " I read quantity=" << quantity << std::endl;
    return false;
  }

  int i = 0;
  if ( mesh != "" )
    in_mesh = true, i++;
  if ( vm.count("polynomial") )
    in_polynomial = true, i++;
  if ( vm.count("input") && ( extension == "vol" ) )
    in_vol = true, i++;
  if ( vm.count("input") && ( ( extension == "obj" ) || ( extension == "OBJ" ) ) )
    in_obj = true, i++;
  if ( i != 1 )
    {
      trace.error() << "Only one input at the same time is possible:"
                    << " either --mesh, --input or --polynomial." << std::endl;
      return false;
    }
  return true;
}

/////////////////////////////////////////////////////////////////////////////
/// Build input from data for further processing.
/////////////////////////////////////////////////////////////////////////////
bool
CurvatureComputer::buildInput()
{
  bool build_ok = false;
  if ( in_vol || in_polynomial ) build_ok = buildPolynomialOrVolInput();
  if ( in_mesh )                 build_ok = buildPredefinedMesh();
  if ( in_obj )                  build_ok = buildObjInput();
  return build_ok;
}

/////////////////////////////////////////////////////////////////////////////
/// Taking care of vol image file or implicit polynomial
/////////////////////////////////////////////////////////////////////////////
bool
CurvatureComputer::buildPolynomialOrVolInput()
{
  trace.beginBlock( "Make Shape from vol file or implicit polynomial" );
  // Generic parameters.
  const auto estimator = vm[ "estimator" ].as<string>();
  params( "gridstep",          h ); // gridstep
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

  /////////////////////////////////////////////////////////////////////////////
  // Case where input is polynomial.
  /////////////////////////////////////////////////////////////////////////////
  if ( in_polynomial )
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
      if ( bimage == nullptr ) in_polynomial = false; // failure
    }
  
  /////////////////////////////////////////////////////////////////////////////
  // Case where input is a 3D image vol file.
  /////////////////////////////////////////////////////////////////////////////
  else if ( in_vol )
    {
      // Fill useful parameters
      params( "thresholdMin", vm[ "thresholdMin" ].as<int>() );
      params( "thresholdMax", vm[ "thresholdMax" ].as<int>() );
      params( "closed",       vm[ "closed" ].as<int>() );
      bimage       = SH::makeBinaryImage( filename, params );
      K            = SH::getKSpace( bimage, params );
      if ( bimage == nullptr ) in_vol = false; // failure
    }
  auto size    = K.upperBound() - K.lowerBound();
  trace.info() << "- Domain size is " << ( size[ 0 ] + 1 )
  	       << " x " << ( size[ 1 ] + 1 )
  	       << " x " << ( size[ 2 ] + 1 ) << std::endl;

  if ( bimage != nullptr )
    {
      unsigned int                nb = 0;
      std::for_each( bimage->cbegin(), bimage->cend(),
		     [&nb] ( bool v ) { nb += v ? 1 : 0; } );
      trace.info() << "- digital shape has " << nb << " voxels." << std::endl;
    }
  auto sembedder   = SH::getSCellEmbedder( K );
  auto embedder    = SH::getCellEmbedder( K );
  if ( bimage != nullptr )  surface     = SH::makeDigitalSurface( bimage, K, params );
  if ( surface != nullptr ) idx_surface = SH::makeIdxDigitalSurface( surface );
  if ( bimage != nullptr && surface == nullptr ) {
    trace.error() << "- surface is empty (either empty or full volume). ";
    trace.info()  << std::endl;
    trace.endBlock();
    return false;
  }
  else if ( surface != nullptr )
    trace.info() << "- surface has " << surface->size()<< " surfels." << std::endl;
  trace.endBlock();

  /////////////////////////////////////////////////////////////////////////////
  // Computing ground truth values if possible.
  /////////////////////////////////////////////////////////////////////////////
  SH::Scalars expected_H_values;
  SH::Scalars expected_G_values;
  SH::Scalars expected_K1_values;
  SH::Scalars expected_K2_values;
  SH::RealVectors expected_D1_values;
  SH::RealVectors expected_D2_values;
  trace.beginBlock( "Compute surfels" );
  params( "surfaceTraversal", "Default" );
  const auto     surfels = ( surface != nullptr )
    ? SH::getSurfelRange( surface, params )
    : SH::SurfelRange();
  trace.endBlock();
  if ( in_polynomial )
    {
      trace.beginBlock( "Compute true curvatures" );
      expected_H_values = SHG::getMeanCurvatures( shape, K, surfels, params );
      expected_G_values = SHG::getGaussianCurvatures( shape, K, surfels, params );
      expected_K1_values= SHG::getFirstPrincipalCurvatures( shape, K, surfels, params );
      expected_K2_values= SHG::getSecondPrincipalCurvatures( shape, K, surfels, params );
      expected_D1_values= SHG::getFirstPrincipalDirections( shape, K, surfels, params );
      expected_D2_values= SHG::getSecondPrincipalDirections( shape, K, surfels, params );
      stat_expected_H_curv.addValues( expected_H_values.cbegin(),
				      expected_H_values.cend() );
      stat_expected_G_curv.addValues( expected_G_values.cbegin(),
				      expected_G_values.cend() );
      stat_expected_H_curv.terminate();
      stat_expected_G_curv.terminate();
      trace.info() << "- truth min_H=" << stat_expected_H_curv.min()
		   << " <= avg_H=" << stat_expected_H_curv.mean()
		   << " <= max_H=" << stat_expected_H_curv.max() << std::endl;
      trace.info() << "- truth min_G=" << stat_expected_G_curv.min()
		   << " <= avg_G=" << stat_expected_G_curv.mean()
		   << " <= max_G=" << stat_expected_G_curv.max() << std::endl;
      time_curv_ground_truth = trace.endBlock();
    }

  /////////////////////////////////////////////////////////////////////////////
  // Compute primal/dual surfaces
  /////////////////////////////////////////////////////////////////////////////
  trace.beginBlock( "Compute primal/dual surface" );
  auto digital_surface_mode = vm[ "digital-surface" ].as<std::string>();
  blocky = digital_surface_mode == "PRIMAL" || digital_surface_mode == "DUAL"
    || ! in_polynomial;
  dual   = digital_surface_mode == "DUAL" || digital_surface_mode == "PDUAL";
  SH::RealPoints  pos_surf; 
  SH::RealPoints  ppos_surf;
  SH::RealPoints  vertices;
  SH::RealVectors normals;
  std::vector< std::vector< Index > > faces;
  SH::Cell2Index c2i;
  if ( ! surfels.empty() )
    { // Getting vertices positions.
      trace.info() << "computing vertices" << std::endl;
      if ( dual )
	{ // dual surface
	  pos_surf = SH::RealPoints( surfels.size() );
	  std::transform( surfels.cbegin(), surfels.cend(), pos_surf.begin(),
			  [&] (const SH::SCell& c) { return h * sembedder( c ); } );
	}
      else
	{ // primal surface
 	  auto pointels = SH::getPointelRange( c2i, surface );
	  pos_surf = SH::RealPoints( pointels.size() );
	  std::transform( pointels.cbegin(), pointels.cend(), pos_surf.begin(),
			  [&] (const SH::Cell& c) { return h * embedder( c ); } ); 
	}
      // project onto true surface if asked.
      trace.info() << "projecting vertices" << std::endl;
      if ( ! blocky ) ppos_surf = SHG::getPositions( shape, pos_surf, params );
      vertices = blocky ? pos_surf : ppos_surf;
      // Build faces
      trace.info() << "build faces" << std::endl;
      if ( dual )
        { // dual surface
	  for ( Index f = 0; f < idx_surface->nbFaces(); ++f )
	    {
	      const auto dual_vtcs = idx_surface->verticesAroundFace( f );
	      std::vector< Index > dual_rvtcs( dual_vtcs.rbegin(), dual_vtcs.rend() );
	      faces.push_back( dual_rvtcs );
	    }
	}
      else
        { // primal surface	  
          for ( auto&& surfel : *surface )
            {
              const auto primal_surfel_vtcs = SH::getPointelRange( K, surfel );
	      std::vector< Index > face;	      
	      for ( auto&& primal_vtx : primal_surfel_vtcs )
		face.push_back( c2i[ primal_vtx ] );
	      faces.push_back( face );
	    }
	}
    }
  trace.endBlock();

  /////////////////////////////////////////////////////////////////////////////
  // Compute true and/or estimated normal vectors
  /////////////////////////////////////////////////////////////////////////////
  if ( in_polynomial )
    {
      trace.beginBlock( "Compute true normals" );
      expected_normals = SHG::getNormalVectors( shape, K, surfels, params );
      trace.endBlock();
    }
  if ( in_polynomial || in_vol )
    {
      trace.beginBlock( "Estimate normals" );
      if ( estimator == "True" && in_polynomial )
	measured_normals = expected_normals;
      else if ( estimator == "CTrivial" )
	measured_normals = SHG::getCTrivialNormalVectors( surface, surfels, params );
      else if ( estimator == "VCM" )
	measured_normals = SHG::getVCMNormalVectors( surface, surfels, params );
      else if ( estimator == "II" ) {
	measured_normals = SHG::getIINormalVectors( bimage, surfels, params );
	auto oriented_n = SHG::getCTrivialNormalVectors( surface, surfels, params );
	SHG::orientVectors( measured_normals, oriented_n );
      }
      else if ( estimator == "Trivial" ) 
	measured_normals = SHG::getTrivialNormalVectors( K, surfels );
      // else "Geometric" normals.
      time_normal_estimations = trace.endBlock();
    }
  // We build a mesh
  trace.beginBlock( "Build mesh from primal/dual surface" );
  smesh.init( vertices.cbegin(), vertices.cend(),
	      faces.cbegin(),    faces.cend() );
  if ( ! measured_normals.empty() ) {
    if ( dual )
      smesh.setVertexNormals( measured_normals.cbegin(), measured_normals.cend() );
    else
      smesh.setFaceNormals( measured_normals.cbegin(), measured_normals.cend() );
  }
  trace.info() << smesh << std::endl;
  trace.endBlock();

  // We complete expected values
  trace.beginBlock( "Build expected values on mesh from primal/dual surface" );
  if ( in_polynomial )
    {
      if ( dual )
	{
	  expected_H_vertex_values = expected_H_values;
	  expected_H_face_values   = smesh.computeFaceValuesFromVertexValues( expected_H_values );
	  expected_G_vertex_values = expected_G_values;
	  expected_G_face_values   = smesh.computeFaceValuesFromVertexValues( expected_G_values );
	  expected_K1_vertex_values= expected_K1_values;
	  expected_K1_face_values  = smesh.computeFaceValuesFromVertexValues( expected_K1_values );
	  expected_K2_vertex_values= expected_K2_values;
	  expected_K2_face_values  = smesh.computeFaceValuesFromVertexValues( expected_K2_values );
	  expected_D1_vertex_values= expected_D1_values;
	  expected_D1_face_values  = smesh.computeFaceUnitVectorsFromVertexUnitVectors( expected_D1_values );
	  expected_D2_vertex_values= expected_D2_values;
	  expected_D2_face_values  = smesh.computeFaceUnitVectorsFromVertexUnitVectors( expected_D2_values );
	}
      else
	{
	  expected_H_vertex_values = smesh.computeVertexValuesFromFaceValues( expected_H_values );
	  expected_H_face_values   = expected_H_values;
	  expected_G_vertex_values = smesh.computeVertexValuesFromFaceValues( expected_G_values );
	  expected_G_face_values   = expected_G_values;
	  expected_K1_vertex_values= smesh.computeVertexValuesFromFaceValues( expected_K1_values );
	  expected_K1_face_values  = expected_K1_values;
	  expected_K2_vertex_values= smesh.computeVertexValuesFromFaceValues( expected_K2_values );
	  expected_K2_face_values  = expected_K2_values;
	  expected_D1_vertex_values= smesh.computeVertexUnitVectorsFromFaceUnitVectors( expected_D1_values );
	  expected_D1_face_values  = expected_D1_values;
	  expected_D2_vertex_values= smesh.computeVertexUnitVectorsFromFaceUnitVectors( expected_D2_values );
	  expected_D2_face_values  = expected_D2_values;
	}
    }
  trace.endBlock();
  
  return in_polynomial || in_vol;
}


/////////////////////////////////////////////////////////////////////////////
/// Taking care of predefined parametric meshes
/////////////////////////////////////////////////////////////////////////////
bool
CurvatureComputer::buildPredefinedMesh()
{
  ASSERT( in_mesh );
  Scalars exp_H, exp_G, exp_K1, exp_K2;
  RealVectors exp_D1, exp_D2;
  trace.beginBlock( "Build predefined parametric mesh" );
  if ( meshargs.size() >= 4
       && ( meshname == "sphere" || meshname == "sphere-VN"
	    || meshname == "sphere-FN" ) )
    {
      const double r = std::stof( meshargs[1] );
      const Index  m = std::stoi( meshargs[2] );
      const Index  n = std::stoi( meshargs[3] );
      const auto normals =
        ( meshname == "sphere-VN" ) ? SimpleMeshHelper::Normals::VERTEX_NORMALS :
        ( meshname == "sphere-FN" ) ? SimpleMeshHelper::Normals::FACE_NORMALS :
        SimpleMeshHelper::Normals::NO_NORMALS;
      smesh = SimpleMeshHelper
        ::makeSphere( r, RealPoint(), m, n, normals );
      exp_H = SimpleMeshHelper::sphereMeanCurvatures    ( r, m, n );
      exp_G = SimpleMeshHelper::sphereGaussianCurvatures( r, m, n );
      exp_K1= SimpleMeshHelper::sphereFirstPrincipalCurvatures ( r, m, n );
      exp_K2= SimpleMeshHelper::sphereSecondPrincipalCurvatures( r, m, n );
      exp_D1= SimpleMeshHelper::sphereFirstPrincipalDirections ( r, m, n );
      exp_D2= SimpleMeshHelper::sphereSecondPrincipalDirections( r, m, n );
      sub_m = m;
      sub_n = n;
    }
  else if ( meshargs.size() >= 5 
	    && ( meshname == "lantern" || meshname == "lantern-VN"
		 || meshname == "lantern-FN" ) )
    {
      const double r = std::stof( meshargs[1] );
      const double h = std::stof( meshargs[2] );
      const Index  m = std::stoi( meshargs[3] );
      const Index  n = std::stoi( meshargs[4] );
      const auto normals =
        ( meshname == "lantern-VN" ) ? SimpleMeshHelper::Normals::VERTEX_NORMALS :
        ( meshname == "lantern-FN" ) ? SimpleMeshHelper::Normals::FACE_NORMALS :
        SimpleMeshHelper::Normals::NO_NORMALS;
      smesh = SimpleMeshHelper
        ::makeLantern( r, h, RealPoint(), m, n, normals );
      exp_H = SimpleMeshHelper::lanternMeanCurvatures    ( r, m, n );
      exp_G = SimpleMeshHelper::lanternGaussianCurvatures( r, m, n );
      exp_K1= SimpleMeshHelper::lanternFirstPrincipalCurvatures ( r, m, n );
      exp_K2= SimpleMeshHelper::lanternSecondPrincipalCurvatures( r, m, n );
      exp_D1= SimpleMeshHelper::lanternFirstPrincipalDirections ( r, m, n );
      exp_D2= SimpleMeshHelper::lanternSecondPrincipalDirections( r, m, n );
      sub_m = m;
      sub_n = n;
    }
  else if ( meshargs.size() >= 5 
	    && ( meshname == "torus" || meshname == "torus-VN"
		 || meshname == "torus-FN" ) )
    {
      const double R = std::stof( meshargs[1] );
      const double r = std::stof( meshargs[2] );
      const Index  m = std::stoi( meshargs[3] );
      const Index  n = std::stoi( meshargs[4] );
      const int twist= std::stoi( meshargs.size() >= 6 ? meshargs[5] : "0" );
      const auto normals =
        ( meshname == "torus-VN" ) ? SimpleMeshHelper::Normals::VERTEX_NORMALS :
        ( meshname == "torus-FN" ) ? SimpleMeshHelper::Normals::FACE_NORMALS :
        SimpleMeshHelper::Normals::NO_NORMALS;
      smesh = SimpleMeshHelper
        ::makeTorus( R, r, RealPoint(), m, n, twist, normals );
      exp_H = SimpleMeshHelper::torusMeanCurvatures    ( R, r, m, n, twist );
      exp_G = SimpleMeshHelper::torusGaussianCurvatures( R, r, m, n, twist );
      exp_K1= SimpleMeshHelper::torusFirstPrincipalCurvatures ( R, r, m, n, twist );
      exp_K2= SimpleMeshHelper::torusSecondPrincipalCurvatures( R, r, m, n, twist );
      exp_D1= SimpleMeshHelper::torusFirstPrincipalDirections ( R, r, m, n, twist );
      exp_D2= SimpleMeshHelper::torusSecondPrincipalDirections( R, r, m, n, twist );
      sub_m = m;
      sub_n = n;
    }
  else
    in_mesh = false;
  if ( in_mesh )
    {
      if ( exp_H.size() == smesh.nbVertices() )	{
	  expected_H_face_values   = smesh.computeFaceValuesFromVertexValues( exp_H );
	  expected_H_vertex_values = exp_H;
      } else {
	expected_H_face_values   = exp_H;
	expected_H_vertex_values = smesh.computeVertexValuesFromFaceValues( exp_H );
      }
      if ( exp_G.size() == smesh.nbVertices() )	{
	  expected_G_face_values   = smesh.computeFaceValuesFromVertexValues( exp_G );
	  expected_G_vertex_values = exp_G;
      } else {
	expected_G_face_values   = exp_G;
	expected_G_vertex_values = smesh.computeVertexValuesFromFaceValues( exp_G );
      }
      if ( exp_K1.size() == smesh.nbVertices() )	{
	  expected_K1_face_values   = smesh.computeFaceValuesFromVertexValues( exp_K1 );
	  expected_K1_vertex_values = exp_K1;
      } else {
	expected_K1_face_values   = exp_K1;
	expected_K1_vertex_values = smesh.computeVertexValuesFromFaceValues( exp_K1 );
      }
      if ( exp_K2.size() == smesh.nbVertices() )	{
	  expected_K2_face_values   = smesh.computeFaceValuesFromVertexValues( exp_K2 );
	  expected_K2_vertex_values = exp_K2;
      } else {
	expected_K2_face_values   = exp_K2;
	expected_K2_vertex_values = smesh.computeVertexValuesFromFaceValues( exp_K2 );
      }
      if ( exp_D1.size() == smesh.nbVertices() )	{
	expected_D1_face_values   = smesh.computeFaceUnitVectorsFromVertexUnitVectors( exp_D1 );
	expected_D1_vertex_values = exp_D1;
      } else {
	expected_D1_face_values   = exp_D1;
	expected_D1_vertex_values = smesh.computeFaceUnitVectorsFromVertexUnitVectors( exp_D1 );
      }
      if ( exp_D2.size() == smesh.nbVertices() )	{
	expected_D2_face_values   = smesh.computeFaceUnitVectorsFromVertexUnitVectors( exp_D2 );
	expected_D2_vertex_values = exp_D2;
      } else {
	expected_D2_face_values   = exp_D2;
	expected_D2_vertex_values = smesh.computeFaceUnitVectorsFromVertexUnitVectors( exp_D2 );
      }
    }
  trace.endBlock();
  return in_mesh;
}

/////////////////////////////////////////////////////////////////////////////
/// Taking care of OBJ input files
/////////////////////////////////////////////////////////////////////////////
bool
CurvatureComputer::buildObjInput()
{
  ASSERT( in_obj );
  // Case where input is a mesh obj file.
  trace.beginBlock( "Reading input obj mesh file" );
  trace.info() << "Reading file <" << filename << ">" << std::endl;
  ifstream mesh_input( filename.c_str() );
  bool  ok = SimpleMeshReader::readOBJ( mesh_input, smesh );
  trace.endBlock();
  if ( ! ok ) {
    trace.error() << "Error reading file <" << filename << ">" << std::endl;
    in_obj = false;
  }
  return in_obj;
}


/////////////////////////////////////////////////////////////////////////////
// Perturbate positions
/////////////////////////////////////////////////////////////////////////////
void
CurvatureComputer::perturbatePositions()
{
  auto uniform_noise  = vm[ "uniform-noise" ].as<double>();  
  auto adaptive_noise = vm[ "adaptive-noise" ].as<double>();  
  if ( uniform_noise > 0.0 )
    {
      const auto eps = uniform_noise * smesh.averageEdgeLength();
      smesh.perturbateWithUniformRandomNoise( eps );
    }
  if ( adaptive_noise > 0.0 )
    smesh.perturbateWithAdaptiveUniformRandomNoise( adaptive_noise );
}

/////////////////////////////////////////////////////////////////////////////
/// Compute normals
/////////////////////////////////////////////////////////////////////////////
void
CurvatureComputer::computeNormals()
{
  trace.beginBlock( "Compute normals if necessary" );
  auto weights_normals_f2v = vm[ "weights-normals-f2v" ].as<std::string>();
  if ( smesh.faceNormals().empty() && smesh.vertexNormals().empty() )
    {
      smesh.computeFaceNormalsFromPositions();
      if ( weights_normals_f2v == "MAX1995" )
        smesh.computeVertexNormalsFromFaceNormalsWithMaxWeights();
      else
        smesh.computeVertexNormalsFromFaceNormals();
    }
  else if ( smesh.faceNormals().empty() )
    smesh.computeFaceNormalsFromVertexNormals();
  else if ( smesh.vertexNormals().empty() )
    {
      if ( weights_normals_f2v == "MAX1995" )
        smesh.computeVertexNormalsFromFaceNormalsWithMaxWeights();
      else
        smesh.computeVertexNormalsFromFaceNormals();
    }
  auto nb_avg_normals = vm[ "average-normals"   ].as<int>();
  for ( int i = 0; i < nb_avg_normals; i++ )
    {
      trace.info() << "face normals -> vertex normals" << std::endl;
      smesh.computeFaceNormalsFromVertexNormals();
      trace.info() << "vertex normals -> face normals" << std::endl;
      if ( weights_normals_f2v == "MAX1995" )
        smesh.computeVertexNormalsFromFaceNormalsWithMaxWeights();
      else
        smesh.computeVertexNormalsFromFaceNormals();
    }
  trace.info() << smesh << std::endl;
  time_normal_estimations += trace.endBlock();
}

/////////////////////////////////////////////////////////////////////////////
// Compute curvatures
/////////////////////////////////////////////////////////////////////////////
bool
CurvatureComputer::computeCurvatures()
{
  bool ok = false;
  if ( quantity == "CNC" )
    ok = computeInterpolatedCorrectedNormalCurrent();
  else if ( quantity == "RZ" )
    ok = computeRusinkiewiczCurvatures();
  else if ( quantity == "NC" )
    ok = computeNormalCycle();
  if ( ! ok ) return false;
  stat_measured_H_curv.addValues( measured_H_face_values.cbegin(),
				  measured_H_face_values.cend() );
  stat_measured_G_curv.addValues( measured_G_face_values.cbegin(),
				  measured_G_face_values.cend() );
  trace.info() << "- estimated min_H=" << stat_measured_H_curv.min()
	       << " <= avg_H=" << stat_measured_H_curv.mean()
	       << " <= max_H=" << stat_measured_H_curv.max() << std::endl;
  trace.info() << "- estimated min_G=" << stat_measured_G_curv.min()
	       << " <= avg_G=" << stat_measured_G_curv.mean()
	       << " <= max_G=" << stat_measured_G_curv.max() << std::endl;
  return true;
}

/////////////////////////////////////////////////////////////////////////////
// Compute interpolated CorrectedNormalCurrent measures
/////////////////////////////////////////////////////////////////////////////
bool
CurvatureComputer::computeInterpolatedCorrectedNormalCurrent()
{
  // Compute local CNC measures
  trace.beginBlock( "Compute local CNC measures" );
  bool unit = vm.count( "unit-normals" );
  trace.info() << "Using " << ( unit ? "unit" : "regular" )
	       << " interpolated normals." << std::endl;
  CNCComputer cnc( smesh );
  cnc.computeInterpolatedMeasures( CNCComputer::Measure::ALL_MU, unit );
  double G = 0.0;
  int nb_nan = 0;
  for ( auto g : cnc.mu2 ) {
    if ( ! isnan( g ) ) G += g;
    else nb_nan++;
  }
  trace.info() << nb_nan << " / " << cnc.mu2.size() << " NaN" << std::endl;
  trace.info() << "Total Gauss curvature G=" << G << std::endl;
  time_mu_estimations = trace.endBlock();

  // Compute CNC measures within balls
  trace.beginBlock( "Computes CNC measures within balls" );
  const double mcoef = vm[ "m-coef"    ].as<double>();
  const double mpow  = vm[ "m-pow"     ].as<double>();
  const double mr    = mcoef * pow( h, mpow );
  trace.info() << "measuring ball radius = " << mr << std::endl;
  measured_ball_mu0.resize ( smesh.nbFaces() );
  measured_ball_mu1.resize ( smesh.nbFaces() );
  measured_ball_mu2.resize ( smesh.nbFaces() );
  measured_ball_muXY.resize( smesh.nbFaces() );
  for ( int f = 0; f < smesh.nbFaces(); ++f )
    {
      trace.progressBar( f, smesh.nbFaces() );
      auto wfaces = smesh.computeFacesInclusionsInBall( mr, f );
      measured_ball_mu0 [ f ] = cnc.interpolatedMu0 ( wfaces );
      measured_ball_mu1 [ f ] = cnc.interpolatedMu1 ( wfaces );
      measured_ball_mu2 [ f ] = cnc.interpolatedMu2 ( wfaces );
      measured_ball_muXY[ f ] = cnc.interpolatedMuXY( wfaces );
    }
  time_curv_estimations = trace.endBlock();

  computeCurvaturesFromMeasures();
  return true;
}  

/////////////////////////////////////////////////////////////////////////////
/// Compute Rusinkiewicz curvatures
/////////////////////////////////////////////////////////////////////////////
bool
CurvatureComputer::computeRusinkiewiczCurvatures()
{
  typedef RusinkiewiczCurvatureFormula< RealPoint, RealVector > RZFormula;
  typedef RZFormula::RealTensor2D   RealTensor2D;
  typedef RZFormula::ColumnVector2D RealVector2D;
  if ( smesh.vertexNormals().empty() )
    {
      trace.warning() << "[CurvatureComputer::computeRusinkiewiczCurvatures]"
                      << " Unable to compute curvatures without vertex normals."
		      << std::endl;
      return false;
    }
  time_mu_estimations = 0.0;
  trace.beginBlock( "Compute Rusinkiewicz's curvatures" );
  measured_H_face_values.resize( smesh.nbFaces() );
  measured_G_face_values.resize( smesh.nbFaces() );
  measured_K1_face_values.resize( smesh.nbFaces() );
  measured_K2_face_values.resize( smesh.nbFaces() );
  measured_D1_face_values.resize( smesh.nbFaces() );
  measured_D2_face_values.resize( smesh.nbFaces() );
  Index idx_f = 0;
  for ( auto f : smesh.incidentVertices() )
    {
      RealPoints  p( f.size() );
      RealVectors u( f.size() );
      for ( Index idx_v = 0; idx_v < f.size(); ++idx_v )
        {
          p[ idx_v ] = smesh.positions()    [ f[ idx_v ] ];
          u[ idx_v ] = smesh.vertexNormals()[ f[ idx_v ] ];
        }
      measured_H_face_values[ idx_f ] = RZFormula::meanCurvature( p, u );
      measured_G_face_values[ idx_f ] = RZFormula::gaussianCurvature( p, u );
      const auto    II = RZFormula::secondFundamentalForm( p, u );
      const auto basis = RZFormula::basis( p );
      auto           V = II;
      RealVector2D   L;
      EigenDecomposition<2, double>::getEigenDecomposition( II, V, L );
      measured_K1_face_values[ idx_f ] = L[ 0 ];      
      measured_K2_face_values[ idx_f ] = L[ 1 ];      
      measured_D1_face_values[ idx_f ] =
	V.column( 0 )[ 0 ] * basis.first + V.column( 0 )[ 1 ] * basis.second;      
      measured_D2_face_values[ idx_f ] =
	V.column( 1 )[ 0 ] * basis.first + V.column( 1 )[ 1 ] * basis.second;      
      idx_f++;
    }
  time_curv_estimations = trace.endBlock();

  return true;
}

/////////////////////////////////////////////////////////////////////////////
// Compute normal cycle  measures
/////////////////////////////////////////////////////////////////////////////
bool
CurvatureComputer::computeNormalCycle()
{
  // Compute local normal cycle measures
  trace.beginBlock( "Compute local normal cycle measures" );
  NCComputer nc( smesh );
  nc.computeMeasures( NCComputer::Measure::ALL_MU );
  {
    double G = 0.0;
    int nb_nan = 0;
    for ( auto g : nc.localMu2 ) {
      if ( ! isnan( g ) ) G += g;
      else nb_nan++;
    }
    trace.info() << nb_nan << " / " << nc.localMu2.size() << " NaN mu2" << std::endl;
    trace.info() << "Total Gauss curvature G=" << G << std::endl;
  }
  {
    int nb_nan = 0;
    for ( auto g : nc.localMu1 ) {
      if ( isnan( g ) ) nb_nan++;
    }
    trace.info() << nb_nan << " / " << nc.localMu1.size() << " NaN mu1" << std::endl;
  }
  {
    int nb_nan = 0;
    for ( auto g : nc.localMu0 ) {
      if ( isnan( g ) ) nb_nan++;
    }
    trace.info() << nb_nan << " / " << nc.localMu0.size() << " NaN mu0" << std::endl;
  }
  time_mu_estimations = trace.endBlock();

  // Compute NC measures within balls
  trace.beginBlock( "Computes Normal cycle measures within balls" );
  const double mcoef = vm[ "m-coef"    ].as<double>();
  const double mpow  = vm[ "m-pow"     ].as<double>();
  const double mr    = mcoef * pow( h, mpow );
  trace.info() << "measuring ball radius = " << mr << std::endl;
  measured_ball_mu0.resize ( smesh.nbFaces() );
  measured_ball_mu1.resize ( smesh.nbFaces() );
  measured_ball_mu2.resize ( smesh.nbFaces() );
  measured_ball_muXY.resize( smesh.nbFaces() );
  for ( int f = 0; f < smesh.nbFaces(); ++f )
    {
      trace.progressBar( f, smesh.nbFaces() );
      auto v_we_wf = smesh.computeCellsInclusionsInBall( mr, f );
      measured_ball_mu0 [ f ] = nc.mu0  ( std::get<2>( v_we_wf ) );
      measured_ball_mu1 [ f ] = nc.mu1  ( std::get<1>( v_we_wf ) );
      measured_ball_mu2 [ f ] = nc.mu2  ( std::get<0>( v_we_wf ) );
      measured_ball_muXY[ f ] = nc.muXY2( std::get<1>( v_we_wf ) );
    }
  time_curv_estimations = trace.endBlock();

  computeCurvaturesFromMeasures();
  return true;
}  


/////////////////////////////////////////////////////////////////////////////
/// If vertex or face curvatures are not computed, compute them from
/// the others.
/////////////////////////////////////////////////////////////////////////////
bool
CurvatureComputer::computeVertexFaceCurvatures()
{
  if ( ! measured_H_face_values.empty() )
    {
      if ( measured_H_vertex_values.empty() )
	measured_H_vertex_values
	  = smesh.computeVertexValuesFromFaceValues( measured_H_face_values );
    }
  else if ( ! measured_H_vertex_values.empty() )
    measured_H_face_values
      = smesh.computeFaceValuesFromVertexValues( measured_H_vertex_values );
  if ( ! measured_G_face_values.empty() )
    {
      if ( measured_G_vertex_values.empty() )
	measured_G_vertex_values
	  = smesh.computeVertexValuesFromFaceValues( measured_G_face_values );
    }
  else if ( ! measured_G_vertex_values.empty() )
    measured_G_face_values
      = smesh.computeFaceValuesFromVertexValues( measured_G_vertex_values );
  if ( ! measured_K1_face_values.empty() )
    {
      if ( measured_K1_vertex_values.empty() )
	measured_K1_vertex_values
	  = smesh.computeVertexValuesFromFaceValues( measured_K1_face_values );
    }
  else if ( ! measured_K1_vertex_values.empty() )
    measured_K1_face_values
      = smesh.computeFaceValuesFromVertexValues( measured_K1_vertex_values );
  if ( ! measured_K2_face_values.empty() )
    {
      if ( measured_K2_vertex_values.empty() )
	measured_K2_vertex_values
	  = smesh.computeVertexValuesFromFaceValues( measured_K2_face_values );
    }
  else if ( ! measured_K2_vertex_values.empty() )
    measured_K2_face_values
      = smesh.computeFaceValuesFromVertexValues( measured_K2_vertex_values );
  if ( ! measured_D1_face_values.empty() )
    {
      if ( measured_D1_vertex_values.empty() )
	measured_D1_vertex_values
	  = smesh.computeVertexUnitVectorsFromFaceUnitVectors( measured_D1_face_values );
    }
  else if ( ! measured_D1_vertex_values.empty() )
    measured_D1_face_values
      = smesh.computeFaceUnitVectorsFromVertexUnitVectors( measured_D1_vertex_values );
  if ( ! measured_D2_face_values.empty() )
    {
      if ( measured_D2_vertex_values.empty() )
	measured_D2_vertex_values
	  = smesh.computeVertexUnitVectorsFromFaceUnitVectors( measured_D2_face_values );
    }
  else if ( ! measured_D2_vertex_values.empty() )
    measured_D2_face_values
      = smesh.computeFaceUnitVectorsFromVertexUnitVectors( measured_D2_vertex_values );
  return true;
}

/////////////////////////////////////////////////////////////////////////////
/// Compute all curvatures from curvature measures, i.e. localize
/// measure computations.
/////////////////////////////////////////////////////////////////////////////
bool
CurvatureComputer::computeCurvaturesFromMeasures()
{
  trace.beginBlock( "Computes curvatures from CNC/NC measures" );
  measured_H_face_values.resize( smesh.nbFaces() );
  measured_G_face_values.resize( smesh.nbFaces() );
  measured_K1_face_values.resize( smesh.nbFaces() );
  measured_K2_face_values.resize( smesh.nbFaces() );
  measured_D1_face_values.resize( smesh.nbFaces() );
  measured_D2_face_values.resize( smesh.nbFaces() );
  bool cnc = ( quantity == "CNC" );
  bool nc  = ( quantity == "NC" );
  for ( SH::Idx i = 0; i < smesh.nbFaces(); i++ )
    {
      if ( isnan( measured_ball_mu0[ i ] ) )
        trace.warning() << "At face " << i << " mu0 is NaN." << std::endl;
      Scalar area = measured_ball_mu0[ i ];
      // Computing mean curvature H
      measured_H_face_values[ i ] = area != 0.0	? measured_ball_mu1[ i ] / area	: 0.0;
      // Computing Gaussian curvature H
      measured_G_face_values[ i ] = area != 0.0	? measured_ball_mu2[ i ] / area	: 0.0;
      // Computing principal curvatures k1 and k2
      auto M = measured_ball_muXY [ i ];
      auto N = smesh.faceNormals()[ i ];
      M += M.transpose();
      M *= 0.5;
      const double   coef_N = 1000.0 * area;
      // Adding 1000 area n x n to anisotropic measure
      for ( int j = 0; j < 3; j++ )
	for ( int k = 0; k < 3; k++ )
	  M( j, k ) += coef_N * N[ j ] * N[ k ];
      auto V = M;
      RealVector L;
      EigenDecomposition< 3, double>::getEigenDecomposition( M, V, L );
      if ( cnc )
	{
	  measured_D1_face_values[ i ] = V.column( 1 );
	  measured_D2_face_values[ i ] = V.column( 0 );
	  measured_K1_face_values[ i ] = area != 0.0 ? -L[ 1 ] / area : 0.0;
	  measured_K2_face_values[ i ] = area != 0.0 ? -L[ 0 ] / area : 0.0;
	}
      else if ( nc )
	{
	  measured_D1_face_values[ i ] = V.column( 0 );
	  measured_D2_face_values[ i ] = V.column( 1 );
	  measured_K1_face_values[ i ] = area != 0.0 ? -L[ 1 ] / area : 0.0;
	  measured_K2_face_values[ i ] = area != 0.0 ? -L[ 0 ] / area : 0.0;
	}
    }
  trace.endBlock();
  return true;
}

/////////////////////////////////////////////////////////////////////////////
/// Output curvatures (G, H, etc) as mesh OBJ files
/////////////////////////////////////////////////////////////////////////////
bool
CurvatureComputer::outputCurvaturesAsMeshObj()
{
  const auto output_basefile  = vm[ "output"   ].as<std::string>();
  if ( output_basefile == "" || output_basefile == "none" ) return false;
  const auto minValue         = vm[ "minValue" ].as<double>();
  const auto maxValue         = vm[ "maxValue" ].as<double>();
  const auto colormap_name    = vm[ "colormap" ].as<std::string>();
  const auto zt               = vm[ "zero-tic" ].as<double>();
  const auto colormapH = SH::getColorMap
    ( minValue, maxValue, params );
  const auto colormapG = SH::getColorMap
    ( ( minValue < 0.0 ? -1.0 : 1.0 ) * 0.5 * minValue*minValue,
      0.5 * maxValue*maxValue, params );

  trace.beginBlock( "Output curvatures in mesh OBJ file" );
  auto output_objfile  = output_basefile + "-primal.obj";
  ofstream mesh_output( output_objfile.c_str() );
  bool okw = SimpleMeshWriter::writeOBJ( mesh_output, smesh );
  mesh_output.close();
  if ( ! okw ) {
    trace.error() << "Error writing file <" << output_objfile << ">" << std::endl;
    return false;
  }
  trace.info() << "Writing mean and Gaussian curvature OBJ files." << std::endl;
  {
    auto colorsH = SH::Colors( smesh.nbFaces() );
    auto colorsG = SH::Colors( smesh.nbFaces() );
    for ( SH::Idx i = 0; i < smesh.nbFaces(); i++ )
      {
	colorsH[ i ] = colormapH( measured_H_face_values[ i ] );
	colorsG[ i ] = colormapG( measured_G_face_values[ i ] );
      }
    okw = SimpleMeshWriter::writeOBJ( output_basefile+"-H", smesh, colorsH );
    okw = SimpleMeshWriter::writeOBJ( output_basefile+"-G", smesh, colorsG );
    if ( zt > 0.0 )
      {
	okw = SimpleMeshWriter::writeIsoLinesOBJ
	  ( output_basefile+"-H-zero", smesh,
	    measured_H_face_values, measured_H_vertex_values, 0.0, zt );
	okw = SimpleMeshWriter::writeIsoLinesOBJ
	  ( output_basefile+"-G-zero", smesh,
	    measured_G_face_values, measured_G_vertex_values, 0.0, zt );
      }
  }
  trace.info() << "Writing K1 and K2 curvature OBJ files." << std::endl;
  {
    auto colorsK1 = SH::Colors( smesh.nbFaces() );
    auto colorsK2 = SH::Colors( smesh.nbFaces() );
    for ( SH::Idx i = 0; i < smesh.nbFaces(); i++ )
      {
	colorsK1[ i ] = colormapH( measured_K1_face_values[ i ] );
	colorsK2[ i ] = colormapH( measured_K2_face_values[ i ] );
      }
    okw = SimpleMeshWriter::writeOBJ( output_basefile+"-K1", smesh, colorsK1 );
    okw = SimpleMeshWriter::writeOBJ( output_basefile+"-K2", smesh, colorsK2 );
    if ( zt > 0.0 )
      {
	okw = SimpleMeshWriter::writeIsoLinesOBJ
	  ( output_basefile+"-K1-zero", smesh,
	    measured_K1_face_values, measured_K1_vertex_values, 0.0, zt );
	okw = SimpleMeshWriter::writeIsoLinesOBJ
	  ( output_basefile+"-K2-zero", smesh,
	    measured_K2_face_values, measured_K2_vertex_values, 0.0, zt );
      }
  }
  trace.info() << "Writing D1 and D2 principal directions OBJ files." << std::endl;
  {
    const auto avg_e = smesh.averageEdgeLength();
    SH::RealPoints positions( smesh.nbFaces() );
    auto d1 = measured_D1_face_values;
    for ( Index f = 0; f < positions.size(); ++f )
      {
	d1[ f ] *= smesh.localWindow( f );
	positions[ f ] = smesh.faceCentroid( f ) - 0.5 * d1[ f ];
      }
    bool ok_d1 = SH::saveVectorFieldOBJ
      ( positions, d1, 0.05 * avg_e, SH::Colors(),
	output_basefile+"-D1-green", SH::Color::Black, SH::Color( 0, 128, 0 ) );
    auto d2 = measured_D2_face_values;
    for ( Index f = 0; f < positions.size(); ++f )
      {
	d2[ f ] *= smesh.localWindow( f );
	positions[ f ] = smesh.faceCentroid( f ) - 0.5 * d2[ f ];
      }
    bool ok_d2 = SH::saveVectorFieldOBJ
      ( positions, d2, 0.05 * avg_e, SH::Colors(),
	output_basefile+"-D2-magenta", SH::Color::Black, SH::Color(128, 0,128 ) );
  }
  trace.info() << "Writing geometry type OBJ files." << std::endl;
  {
    auto colors = SH::Colors( smesh.nbFaces() );
    for ( SH::Idx i = 0; i < smesh.nbFaces(); i++ )
      {
	colors[ i ] = computeGeometryColor( measured_K1_face_values[ i ],
					    measured_K2_face_values[ i ],
					    0.01, maxValue );
      }
    okw = SimpleMeshWriter::writeOBJ( output_basefile+"-GT", smesh, colors );
  }
  trace.endBlock();
  return okw;
}


/////////////////////////////////////////////////////////////////////////////
/// Output Curvature errors (G, H, etc) as mesh OBJ files
/////////////////////////////////////////////////////////////////////////////
bool
CurvatureComputer::outputCurvatureErrorsAsMeshObj()
{
  bool okw = false;
  if ( ! ( in_polynomial || in_mesh ) ) return okw;
  const auto output_basefile  = vm[ "output"   ].as<std::string>();
  if ( output_basefile == "" || output_basefile == "none" ) return false;
  const auto max_error    = vm[ "max-error" ].as<double>();
  const auto error_cmap   = getErrorColorMap( max_error );
  trace.beginBlock( "Output errors in mesh OBJ file" );
  // Error for mean curvature
  const auto error_H_values =
    SHG::getScalarsAbsoluteDifference( measured_H_face_values,
				       expected_H_face_values );
  stat_error_H_curv       = SHG::getStatistic( error_H_values );
  std::vector<Color> colorsHE( smesh.nbFaces() );
  for ( SH::Idx i = 0; i < colorsHE.size(); i++ )
    colorsHE[ i ] = error_cmap( error_H_values[ i ] ); 
  okw = SimpleMeshWriter::writeOBJ( output_basefile+"-H-error.obj",
				    smesh, colorsHE );
  const auto error_H_l2 = SHG::getScalarsNormL2( measured_H_face_values,
						 expected_H_face_values );
  trace.info() << "|H-H_CNC|_oo = " << stat_error_H_curv.max() << std::endl;
  trace.info() << "|H-H_CNC|_2  = " << error_H_l2 << std::endl;

  // Error for gaussian curvature
  const auto error_G_values =
    SHG::getScalarsAbsoluteDifference( measured_G_face_values,
				       expected_G_face_values );
  stat_error_G_curv       = SHG::getStatistic( error_G_values );
  std::vector<Color> colorsGE( smesh.nbFaces() );
  for ( SH::Idx i = 0; i < colorsGE.size(); i++ )
    colorsGE[ i ] = error_cmap( error_G_values[ i ] ); 
  okw = SimpleMeshWriter::writeOBJ( output_basefile+"-G-error.obj",
				    smesh, colorsGE );
  const auto error_G_l2 = SHG::getScalarsNormL2( measured_G_face_values,
						 expected_G_face_values );
  trace.info() << "|G-G_CNC|_oo = " << stat_error_G_curv.max() << std::endl;
  trace.info() << "|G-G_CNC|_2  = " << error_G_l2 << std::endl;
  trace.endBlock();
  return okw;
}

/////////////////////////////////////////////////////////////////////////////
/// Output Curvature errors statistics.
/////////////////////////////////////////////////////////////////////////////
bool
CurvatureComputer::outputCurvatureErrorStatistics()
{
  auto bname = vm[ "error" ].as<std::string>();
  trace.beginBlock( "Output statistics" );
  {
    std::ofstream ferr;
    ferr.open ( ( bname+"-stat-H.txt" ).c_str(),
		std::ofstream::out | std::ofstream::app );
    if ( ! ferr.good() )
      trace.warning() << "Unable to open file: " << bname << "-stat-H.txt" << std::endl;
    outputScalarCurvatureErrorStatistics( ferr, command_line,
					  measured_H_face_values,
					  expected_H_face_values );
  }
  {
    std::ofstream ferr;
    ferr.open ( ( bname+"-stat-G.txt" ).c_str(),
		std::ofstream::out | std::ofstream::app );
    if ( ! ferr.good() )
      trace.warning() << "Unable to open file: " << bname << "-stat-G.txt" << std::endl;
    outputScalarCurvatureErrorStatistics( ferr, command_line,
					  measured_G_face_values,
					  expected_G_face_values );
  }
  {
    std::ofstream ferr;
    ferr.open ( ( bname+"-stat-K1.txt" ).c_str(),
		std::ofstream::out | std::ofstream::app );
    if ( ! ferr.good() )
      trace.warning() << "Unable to open file: " << bname << "-stat-K1.txt" << std::endl;
    outputScalarCurvatureErrorStatistics( ferr, command_line,
					  measured_K1_face_values,
					  expected_K1_face_values );
  }
  {
    std::ofstream ferr;
    ferr.open ( ( bname+"-stat-K2.txt" ).c_str(),
		std::ofstream::out | std::ofstream::app );
    if ( ! ferr.good() )
      trace.warning() << "Unable to open file: " << bname << "-stat-K2.txt" << std::endl;
    outputScalarCurvatureErrorStatistics( ferr, command_line,
					  measured_K2_face_values,
					  expected_K2_face_values );
  }
  trace.endBlock();
  return true;
}

bool
CurvatureComputer::outputScalarCurvatureErrorStatistics
( std::ostream& output, std::string information,
  Scalars measured_values, Scalars expected_values )
{
  if ( measured_values.size() == 0 ) measured_values = expected_values;
  if ( expected_values.size() == 0 ) expected_values = measured_values;
  output << "######################################################################"
	 << std::endl;
  output << "# " << information << std::endl;
  output << "# time_normal_estimations = "
	 << time_normal_estimations << " ms" << std::endl
	 << "# time_curv_estimations   = "
	 << time_curv_estimations + time_mu_estimations << " ms" << std::endl
	 << "# time_mu_estimations     = "
	 << time_mu_estimations << " ms" << std::endl;
  output << "# grid_h(1) subdi_m(2) subdi_n(3)"
	 << " #vertic(4)  #edges(5)  #faces(6)"
	 << "  err_L1(7)  err_L2(8) err_Loo(9)"
	 << " ex_avg(10) ex_dev(11) ex_min(12) ex_max(13)"
	 << " ms_avg(14) ms_dev(15) ms_min(16) ms_max(17)"
	 << " er_avg(18) er_dev(19) er_min(20) er_max(21)";
  output << " rho_rl(22) rho_dg(23) t_norm(24) t_curv(25)"
	 << "   t_mu(26)" << std::endl;
  auto err_values = SHG::getScalarsAbsoluteDifference( measured_values, expected_values );
  Statistic<Scalar> stats_meas;
  Statistic<Scalar> stats_exps;
  Statistic<Scalar> stats_errs;
  stats_meas.addValues( measured_values.cbegin(), measured_values.cend() );
  stats_exps.addValues( expected_values.cbegin(), expected_values.cend() );
  stats_errs.addValues( err_values.cbegin(), err_values.cend() );

  output << std::setprecision(6) << std::fixed; 
  output << std::setw(10) << h
	 << " " << std::setw(10) << sub_m << " " << std::setw(10) << sub_n;
  output << " " << std::setw(10) << smesh.nbVertices()
	 << " " << std::setw(10) << smesh.nbEdges()
	 << " " << std::setw(10) << smesh.nbFaces();
  output << " " << std::setw(10) << SHG::getScalarsNormL1 ( measured_values, expected_values )
	 << " " << std::setw(10) << SHG::getScalarsNormL2 ( measured_values, expected_values )
	 << " " << std::setw(10) << SHG::getScalarsNormLoo( measured_values, expected_values );
  output << " " << std::setw(10) << stats_exps.mean()
	 << " " << std::setw(10) << sqrt( stats_exps.variance() )
	 << " " << std::setw(10) << stats_exps.min()
	 << " " << std::setw(10) << stats_exps.max();
  output << " " << std::setw(10) << stats_meas.mean()
	 << " " << std::setw(10) << sqrt( stats_meas.variance() )
	 << " " << std::setw(10) << stats_meas.min()
	 << " " << std::setw(10) << stats_meas.max();
  output << " " << std::setw(10) << stats_errs.mean()
	 << " " << std::setw(10) << sqrt( stats_errs.variance() )
	 << " " << std::setw(10) << stats_errs.min()
	 << " " << std::setw(10) << stats_errs.max();
  const double rho_coef = vm[ "m-coef" ].as<double>();
  const double rho_pow  = vm[ "m-pow" ].as<double>();
  const double rho_real = rho_coef * pow( h, rho_pow );
  const double rho_dig  = rho_real / h;
  output << " " << std::setw(10) << rho_real
	 << " " << std::setw(10) << rho_dig 
	 << " " << std::setw(10) << time_normal_estimations
	 << " " << std::setw(10) << (time_curv_estimations + time_mu_estimations)
	 << " " << std::setw(10) << time_mu_estimations << std::endl;
  output << "#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
	 << std::endl;
  return true;
}

CurvatureComputer::ColorMap
CurvatureComputer::computeGeometryTypeColormap()
{
  ColorMap gradcmap( -1.0, 1.0 );
  gradcmap.addColor( Color( 0, 0, 255 ) ); // concave -1
  // BLUE to BLUE   1.0 <= |k1/k2| <= 0.2, k1 <= k2  < 0.0
  gradcmap.addColor( Color(   0,   0, 255 ) ); // concave -0.9
  gradcmap.addColor( Color(   0,   0, 255 ) ); // concave -0.8
  gradcmap.addColor( Color(   0,   0, 255 ) ); // concave -0.7
  gradcmap.addColor( Color(   0,   0, 255 ) ); // concave -0.6
  // BLUE to CYAN   0.2 <= |k1/k2| <= 0.0, k1 <= k2  < 0.0
  gradcmap.addColor( Color(   0, 255, 255 ) ); // cyl ccv -0.5, k2 approx 0
  // CYAN to GREEN  0.0 <= |k1/k2| <= 0.2, k1 <  0.0 < k2
  gradcmap.addColor( Color(   0, 255,   0 ) ); // hyperbo -0.4
  // GREEN to GREEN  0.2 <= min(|k2/k1|,|k1/k2|) <= 1, k1 <  0.0 < k2
  gradcmap.addColor( Color(   0, 255,   0 ) ); // hyperbo -0.4
  gradcmap.addColor( Color(   0, 255,   0 ) ); // hyperbo -0.3
  gradcmap.addColor( Color(   0, 255,   0 ) ); // hyperbo -0.2
  gradcmap.addColor( Color(   0, 255,   0 ) ); // hyperbo -0.1
  gradcmap.addColor( Color(   0, 128,   0 ) ); // hyperbo  0.0
  gradcmap.addColor( Color(   0, 255,   0 ) ); // hyperbo  0.1
  gradcmap.addColor( Color(   0, 255,   0 ) ); // hyperbo  0.2
  gradcmap.addColor( Color(   0, 255,   0 ) ); // hyperbo  0.3
  gradcmap.addColor( Color(   0, 255,   0 ) ); // hyperbo  0.4
  // GREEN TO YELLOW
  gradcmap.addColor( Color( 255, 255,   0 ) ); // cyl cvx  0.5, k1 approx 0
  // YELLOW TO RED   0.2 <= |k1/k2| <= 1.0, 0.0 < k1 <= k2
  gradcmap.addColor( Color( 255,   0,   0 ) ); // convex   0.6
  gradcmap.addColor( Color( 255,   0,   0 ) ); // convex   0.7
  gradcmap.addColor( Color( 255,   0,   0 ) ); // convex   0.8
  gradcmap.addColor( Color( 255,   0,   0 ) ); // convex   0.9
  gradcmap.addColor( Color( 255,   0,   0 ) ); // convex   1.0
  return gradcmap;
}

/// Given curvature \a k1 and \a k2, `k1 <= k2`, a near zero value \a
/// zero, and a maximum value \a max, returns an associated color
/// corresponding to its geometry type (convex, concave, hyperbolic,
/// cylindric or flat).
Color
CurvatureComputer::computeGeometryColor( Scalar k1, Scalar k2,
					 Scalar zero, Scalar max ) const
{
  if ( k2 < k1 ) std::swap( k1, k2 );
  Scalar M = std::max( fabs( k1 ), fabs( k2 ) );
  Scalar m = std::min( fabs( k1 ), fabs( k2 ) );
  Scalar r = ( M <= zero ) ? 0.0 : m / M;
  Color flat ( Color::White );
  Color target( Color::Black );
  if ( fabs( k1 ) <= zero && fabs( k2 ) <= zero ) return flat; // flat
  if ( k2 < -zero ) // k1 <= k2 < 0.0, r = |k2| / |k1|
    target = geom_cmap( std::max( -0.5 - r/2.0, -1.0 ) );
  else if ( k2 <= zero ) // k1 < 0.0, k2 = 0
    target = geom_cmap( -0.5 );
  else if ( k1 < -zero ) 
    target = ( fabs(k2) <= fabs(k1) )
      ? geom_cmap( -0.5 + r/2.0 )
      : geom_cmap(  0.5 - r/2.0 );
  else if ( k1 <= zero ) // k1 = 0, k2 > 0.0 
    target = geom_cmap( 0.5 );
  else
    target = geom_cmap( std::max( 0.5 + r/2.0, 1.0 ) );
  Scalar s = std::min( M / max, 1.0 );
  s = round( s * 8.0 ) / 8.0;
  return (1.0 - s) * flat + s * target;
}


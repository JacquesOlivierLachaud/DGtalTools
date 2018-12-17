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
 * @file at-u2-v0.cpp
 * @ingroup Tools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 * @author Marion Foare (\c marion.foare@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2016/10/12
 *
 * A tool file named at-u2-v0.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <string>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
// always include EigenSupport.h before any other Eigen headers
#include "DGtal/math/linalg/EigenSupport.h"
#include "DGtal/dec/DiscreteExteriorCalculus.h"
#include "DGtal/dec/DiscreteExteriorCalculusSolver.h"
#include "DGtal/dec/DiscreteExteriorCalculusFactory.h"

#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/readers/GenericReader.h"

using namespace DGtal;
using namespace std;

template <typename Calculus>
typename Calculus::PrimalIdentity1
proj1( const Calculus& dec, Dimension k )
{
  auto px = dec.template identity<1,PRIMAL>();
  typedef typename Calculus::Index Index;
  typename Calculus::PrimalForm1 f( dec );
  for ( int j = 0; j < px.myContainer.outerSize(); ++j )
    for ( typename Calculus::SparseMatrix::InnerIterator it( px.myContainer, j ); it; ++it )
      {
        Index i = it.row();
        auto slinel = f.getSCell( i );
        if ( ! dec.myKSpace.sIsOpen( slinel, k ) )
          it.valueRef() = 0.0;
      }
  px.myContainer.makeCompressed();
  return px;
}

template <typename Calculus>
typename Calculus::DualIdentity1
proj1_( const Calculus& dec, Dimension k )
{
  auto px = dec.template identity<1,DUAL>();
  typedef typename Calculus::Index Index;
  typename Calculus::PrimalForm1 f( dec );
  for ( int j = 0; j < px.myContainer.outerSize(); ++j )
    for ( typename Calculus::SparseMatrix::InnerIterator it( px.myContainer, j ); it; ++it )
      {
        Index i = it.row();
        auto slinel = f.getSCell( i );
        if ( dec.myKSpace.sIsOpen( slinel, k ) ) // inversion due to duality
          it.valueRef() = 0.0;
      }
  px.myContainer.makeCompressed();
  return px;
}


template <typename Calculus, typename Form1>
Form1 dx( const Calculus& dec, const Form1& alpha1 )
{
  Form1 result = alpha1;
  typedef typename Calculus::Index Index;
  for ( Index i = 0; i < alpha1.length(); ++i )
    {
      auto slinel = alpha1.getSCell( i );
      if ( dec.myKSpace.sIsOpen( slinel, 1 ) ) // vertical
        result.myContainer[ i ] = 0.0;
    }
  return result;
}

template <typename Calculus, typename Form1>
Form1 dy( const Calculus& dec, const Form1& alpha1 )
{
  Form1 result = alpha1;
  typedef typename Calculus::Index Index;
  for ( Index i = 0; i < alpha1.length(); ++i )
    {
      auto slinel = alpha1.getSCell( i );
      if ( dec.myKSpace.sIsOpen( slinel, 0 ) ) // horizontal
        result.myContainer[ i ] = 0.0;
    }
  return result;
}

// template <typename Calculus>
// typename Calculus::PrimalAntiderivative1
// Kxx( const Calculus& dec, const typename Calculus::DualForm1& alpha )


int main( int argc, char* argv[] )
{
  const Z2i::Domain domain(Z2i::Point(0,0), Z2i::Point(9,9));
  Z2i::DigitalSet aSet( domain );
  for ( auto p : domain ) aSet.insert( p );
  // create discrete exterior calculus from set
  //! [calculus_creation]
  typedef DiscreteExteriorCalculus<2, 2, EigenLinearAlgebraBackend> Calculus;
  typedef DiscreteExteriorCalculusFactory<EigenLinearAlgebraBackend> CalculusFactory;
  Calculus calculus = CalculusFactory::createFromDigitalSet( aSet, true );
  //! [calculus_creation]
  trace.info() << calculus << endl;
  
  typedef typename Calculus::Index                       Index;
  typedef typename Calculus::PrimalForm0                 PrimalForm0;
  typedef typename Calculus::PrimalForm1                 PrimalForm1;
  typedef typename Calculus::PrimalForm2                 PrimalForm2;
  typedef typename Calculus::PrimalIdentity0             PrimalIdentity0;
  typedef typename Calculus::PrimalIdentity1             PrimalIdentity1;
  typedef typename Calculus::PrimalIdentity2             PrimalIdentity2;
  typedef typename Calculus::DualIdentity0             DualIdentity0;
  typedef typename Calculus::DualIdentity1             DualIdentity1;
  typedef typename Calculus::DualIdentity2             DualIdentity2;
  typedef typename Calculus::PrimalDerivative0           PrimalDerivative0;
  typedef typename Calculus::PrimalDerivative1           PrimalDerivative1;
  typedef typename Calculus::DualDerivative0             DualDerivative0;
  typedef typename Calculus::DualDerivative1             DualDerivative1;
  typedef typename Calculus::PrimalAntiderivative1       PrimalAntiderivative1;
  typedef typename Calculus::PrimalAntiderivative2       PrimalAntiderivative2;
  typedef typename Calculus::DualAntiderivative1         DualAntiderivative1;
  typedef typename Calculus::DualAntiderivative2         DualAntiderivative2;
  typedef typename Calculus::PrimalHodge0                PrimalHodge0;
  typedef typename Calculus::PrimalHodge1                PrimalHodge1;
  typedef typename Calculus::PrimalHodge2                PrimalHodge2;
  typedef typename Calculus::DualHodge0                  DualHodge0;
  typedef typename Calculus::DualHodge1                  DualHodge1;
  typedef typename Calculus::DualHodge2                  DualHodge2;
  PrimalDerivative0 d0  = calculus.template derivative<0,PRIMAL>();
  PrimalDerivative1 d1  = calculus.template derivative<1,PRIMAL>();
  DualDerivative0   d0_ = calculus.template derivative<0,DUAL>();
  DualDerivative1   d1_ = calculus.template derivative<1,DUAL>();
  PrimalHodge0      h0  = calculus.template hodge<0,PRIMAL>();
  PrimalHodge1      h1  = calculus.template hodge<1,PRIMAL>();
  PrimalHodge2      h2  = calculus.template hodge<2,PRIMAL>();
  DualHodge0        h0_ = calculus.template hodge<0,DUAL>();
  DualHodge1        h1_ = calculus.template hodge<1,DUAL>();
  DualHodge2        h2_ = calculus.template hodge<2,DUAL>();
  
  Calculus::DualForm0 dirac(calculus);
  dirac.myContainer(calculus.getCellIndex( calculus.myKSpace.uSpel(Z2i::Point(4,4))) ) = 1;
  Calculus::DualIdentity0 laplace = calculus.laplace<DUAL>() + 0.01 * calculus.identity<0, DUAL>();
  
  typedef EigenLinearAlgebraBackend::SolverSimplicialLDLT LinearAlgebraSolver;
  typedef DiscreteExteriorCalculusSolver<Calculus, LinearAlgebraSolver, 0, DUAL, 0, DUAL> Solver;
  
  Solver solver;
  solver.compute(laplace);
  Calculus::DualForm0 solution = solver.solve(dirac);
  trace.info() << solver.isValid() << " " << solver.myLinearAlgebraSolver.info() << endl;
  trace.info() << solution << endl;

  auto X = proj1_( calculus, 0 );
  auto Y = proj1_( calculus, 1 );

  Calculus::DualVectorField vfu(calculus);
  for (Index ii=0; ii<vfu.length(); ii++)
    {
      const Z2i::RealPoint cell_center = Z2i::RealPoint(vfu.getSCell(ii).preCell().coordinates)/2.;
      vfu.myCoordinates(ii, 0) = cos(-.5*cell_center[0]+ .3*cell_center[1]);
      vfu.myCoordinates(ii, 1) = cos(.4*cell_center[0]+ .8*cell_center[1]);
    }
  trace.info() << vfu << endl;
  Calculus::DualForm0 dux_dx(calculus);
  Calculus::DualForm0 duy_dy(calculus);
  for (Index ii=0; ii<dux_dx.length(); ii++)
    {
      const Z2i::RealPoint cell_center = Z2i::RealPoint(dux_dx.getSCell(ii).preCell().coordinates)/2.;
      dux_dx.myContainer[ ii ] =  0.5 * sin( -.5*cell_center[0] + .3*cell_center[1] );
      duy_dy.myContainer[ ii ] = -0.8 * sin(  .4*cell_center[0] + .8*cell_center[1] );
    }
  trace.info() << dux_dx << duy_dy << endl;
  Calculus::PrimalForm0 duy_dx_dux_dy(calculus);
  for (Index ii=0; ii<duy_dx_dux_dy.length(); ii++)
    {
      const Z2i::RealPoint cell_center = Z2i::RealPoint(duy_dx_dux_dy.getSCell(ii).preCell().coordinates)/2.;
      duy_dx_dux_dy.myContainer[ ii ] =  -0.3 * sin( -.5*cell_center[0] + .3*cell_center[1] )
        - 0.4 * sin(  .4*cell_center[0] + .8*cell_center[1] );
    }
  trace.info() << duy_dx_dux_dy << endl;

  
  Calculus::DualForm1 u = calculus.flat( vfu ); // d0_ * solution;
  Calculus::DualForm1 ux = X * u;
  Calculus::DualForm1 uy = Y * u;
  {
    Board2D board;
    board << domain;
    board << CustomStyle("VectorField", new VectorFieldStyle2D(1.0));
    board << vfu << u;
    board.saveSVG("solve_tensor_laplace.svg");
  }
  DualAntiderivative1 Kxx = h2 * d1 * h1_ * X;
  DualAntiderivative1 Kyy = h2 * d1 * h1_ * Y;
  LinearOperator<Calculus, 1, DUAL, 0, PRIMAL> Kxy = 0.5 * h2_ * ( d1_ * Y - 1.0 * d1_ * X );
  LinearOperator<Calculus, 1, DUAL, 2, DUAL>   Kxy_ = 0.5 * ( d1_ * Y - 1.0 * d1_ * X );
  {
    Board2D board;
    board << domain;
    board << ( Kxx * u );
    board.saveSVG("solve_tensor_Kxx.svg");
  }
  {
    Board2D board;
    board << domain;
    board << dux_dx;
    board.saveSVG("solve_tensor_Kxx_truth.svg");
  }
  {
    Board2D board;
    board << domain;
    board << ( Kyy * u );
    board.saveSVG("solve_tensor_Kyy.svg");
  }
  {
    Board2D board;
    board << domain;
    board << duy_dy;
    board.saveSVG("solve_tensor_Kyy_truth.svg");
  }
  {
    Board2D board;
    board << domain;
    board << ( Kxy_ * u );
    board.saveSVG("solve_tensor_Kxy.svg");
  }
  {
    Board2D board;
    board << domain;
    board << duy_dx_dux_dy;
    board.saveSVG("solve_tensor_Kxy_truth.svg");
  }

  // Compute trace( K )^2
  DualAntiderivative1 TrK   = Kxx + Kyy;
  auto                TrK_t = TrK.transpose();
  DualIdentity1       TrK2  = TrK_t * TrK;
  trace.info() << "Tr(K)^2 = " << u.myContainer.dot( ( TrK2 * u ).myContainer ) << endl;
  DualIdentity1       Kxx2 = Kxx.transpose() * Kxx;
  DualIdentity1       Kyy2 = Kyy.transpose() * Kyy;
  DualIdentity1       Kxy2 = Kxy_.transpose() * Kxy_;
  DualIdentity1       KdK  = Kxx2 + Kyy2 + 2.0 * Kxy2;
  trace.info() << "K.K     = " << u.myContainer.dot( ( KdK * u ).myContainer ) << endl;
  
  return 0;
}

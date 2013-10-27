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
 * @file 2dContourAnalysis.cpp
 * @ingroup Tools
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 *
 * @date 2013/10/23
 *
 * DGtal tool for analyzing 2d contours.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/base/Common.h"
#include "DGtal/base/Clock.h"

//space / domain
#include "DGtal/kernel/SpaceND.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/helpers/StdDefs.h"

//Grid curve
#include "DGtal/geometry/curves/FreemanChain.h"
#include "DGtal/geometry/curves/GridCurve.h"
#include "DGtal/geometry/curves/ArithmeticalDSS.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include "DGtal/math/Statistic.h"
#include "DGtal/math/Histogram.h"

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

using namespace DGtal;

template <typename KSpace, typename Iterator>
void analyseLengthMS( Statistic<double> & statD, Statistic<double> & statE, 
                      Iterator itb, Iterator ite )
{
  typedef typename KSpace::Space Space;
  typedef typename Space::Point Point;
  typedef typename Space::Vector Vector;
  typedef ArithmeticalDSS<Iterator,int,4> SegmentComputer;
  typedef SaturatedSegmentation<SegmentComputer> Decomposition;
  typedef typename Decomposition::SegmentComputerIterator SegmentComputerIterator;
  // Computes the tangential cover
  SegmentComputer algo;
  Decomposition theDecomposition( itb, ite, algo);
  statD.clear();
  statE.clear();
  for ( SegmentComputerIterator scIt = theDecomposition.begin(), scItEnd = theDecomposition.end();
	scIt != scItEnd; ++scIt )
    {
      const SegmentComputer & sc = *scIt;
      int64_t l = 0; 
      for ( Iterator ptIt = sc.begin(), ptItEnd = sc.end(); ptIt != ptItEnd; ++ptIt )
	++l;
      statD.addValue( (double) l );
      Vector v = *( sc.end() - 1 ) - *( sc.begin() );
      statE.addValue( v.norm() );
    }
}

double kappaGridStep( double length )
{
  return 5.0*5.0*5.0/(length*length*length);
}

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
    ("help,h", "display this message")
    ("freeman-chain,f", po::value<std::string>(), "specifies the input contour as a freeman chain file name")
    ("histogram-discrete-length,d", po::value<unsigned int>(), "outputs an histogram with <arg> bins of the discrete length of maximal segments.")  
    ("histogram-euclidean-length,e", po::value<unsigned int>(), "outputs an histogram with <arg> bins of the Euclidean length of maximal segments.")  
    ("histogram-around-mean,m", po::value<double>(), "indicates that the histogram should be in interval [mu-<arg>*dev,mu+<arg>*dev] instead of [min,max], where mu is the mean and dev is the standard deviation.")
    ("stats,s", "outputs statistics on the length of maximal segments along the contour.")
    ;
  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
  }
  po::notify(vm);    
  if ( ( ! parseOK ) || vm.count("help") || ( argc <= 1 ) )
    {
      trace.info() << "Analyse a 2d input contour with maximal segments so as to extract some characteristics, like digitization gridstep, range of curvatures." << std::endl 
                   << "Basic usage: " <<std::endl
                   << "\t2dContourAnalysis [options] -s --freeman-chain <file.fc>"<<std::endl
                   << general_opt << std::endl;
      trace.info() << "Outputed statistics are:" << std::endl
                   << "\t - i       : the index of the freeman chain curve." << std::endl
                   << "\t - Nb      : the number of maximal segments." << std::endl
                   << "\t - AvgL    : the average of the discrete length of maximal segments." << std::endl
                   << "\t - UVarL   : the unbiased variance of the discrete length of maximal segments." << std::endl
                   << "\t - MinL    : the minimal discrete length of maximal segments." << std::endl
                   << "\t - MaxL    : the maximal discrete length of maximal segments." << std::endl
                   << "\t - MedL    : the median of the discrete length of maximal segments." << std::endl
                   << "\t - AvgEL   : the average of the Euclidean length of maximal segments." << std::endl
                   << "\t - UVarEL  : the unbiased variance of the Euclidean length of maximal segments." << std::endl
                   << "\t - MinEL   : the minimal Euclidean length of maximal segments." << std::endl
                   << "\t - MaxEL   : the maximal Euclidean length of maximal segments." << std::endl
                   << "\t - MedEL   : the median of the Euclidean length of maximal segments." << std::endl
                   << "\t - kh      : the estimated kappa * h obtained with MedL, kappa the average curvature, h the gridstep. " << std::endl
                   << "\t - 1/kh    : the estimated 1/(kappa * h) obtained with MedL " << std::endl
                   << "\t - khE     : the estimated kappa * h obtained with MedEL, kappa the average curvature, h the gridstep. " << std::endl
                   << "\t - 1/khE   : the estimated 1/(kappa * h) obtained with MedEL " << std::endl
                   << "\t - Avg(1/kh): the estimated 1/(kappa * h) obtained with averaging over all discrete lengths." << std::endl
                   << "\t - Avg(1/khE): the estimated 1/(kappa * h) obtained with averaging over all euclidean lengths." << std::endl;
      return 0;
    }
  
  std::string fileName = vm.count("freeman-chain") 
    ? vm["freeman-chain"].as<std::string>()
    : "/dev/stdin";

  typedef Z2i::Space Space; 
  typedef Space::Point Point; 
  typedef Space::Integer Integer;  
  typedef Z2i::KSpace KSpace; 
  typedef FreemanChain<Integer> FreemanChain; 

  std::vector< FreemanChain > vectFcs = 
    PointListReader< Point >::getFreemanChainsFromFile<Integer> (fileName);  

  for(unsigned int i=0; i<vectFcs.size(); i++){
    
    // Freeman chain
    FreemanChain fc = vectFcs.at(i); 
    // Create GridCurve
    GridCurve< KSpace > gridcurve;
    gridcurve.initFromPointsRange( fc.begin(), fc.end() );
    typedef GridCurve<KSpace>::PointsRange Range;
    Range r = gridcurve.getPointsRange();//building range
   
    trace.info() << "# grid curve " << i+1 << "/" 
		 << gridcurve.size() << " "
		 << ( (gridcurve.isClosed())?"closed":"open" ) << std::endl;
    Statistic<double> statMSL( true );
    Statistic<double> statMSEL( true );
    if ( gridcurve.isClosed() )
      analyseLengthMS<KSpace>( statMSL, statMSEL, r.c(), r.c() );
    else 
      analyseLengthMS<KSpace>( statMSL, statMSEL, r.begin(), r.end() );
    statMSL.terminate();
    Statistic<double> resolution( false );
    for ( Statistic<double>::ConstIterator 
            it = statMSL.begin(), itE = statMSL.end(); it != itE; ++it )
      resolution.addValue( 1.0 / kappaGridStep( *it ) );
    Statistic<double> resolutionE( false );
    for ( Statistic<double>::ConstIterator 
            it = statMSEL.begin(), itE = statMSEL.end(); it != itE; ++it )
      resolutionE.addValue( 1.0 / kappaGridStep( *it ) );
    if ( vm.count( "stats" ) )
      {
        std::cout << "# i-th Nb AvgL UVarL MinL MaxL MedL";
        std::cout << " AvgEL UVarEL MinEL MaxEL MedEL";
        std::cout << " kh 1/kh khE 1/khE Avg(1/kh) Avg(1/khE)" << std::endl;
        std::cout << i 
                  << " " << statMSL.samples() 
                  << " " << statMSL.mean() 
                  << " " << statMSL.unbiasedVariance() 
                  << " " << statMSL.min() << " " << statMSL.max()
                  << " " << statMSL.median() 
                  << " " << statMSEL.mean() 
                  << " " << statMSEL.unbiasedVariance() 
                  << " " << statMSEL.min() << " " << statMSEL.max()
                  << " " << statMSEL.median() 
                  << " " << kappaGridStep( statMSL.median() )
                  << " " << 1.0 / kappaGridStep( statMSL.median() )
                  << " " << kappaGridStep( statMSEL.median() )
                  << " " << 1.0 / kappaGridStep( statMSEL.median() )
                  << " " << resolution.mean() 
                  << " " << resolutionE.mean() 
                  << std::endl;
      }
    if ( vm.count( "histogram-discrete-length" ) )
      {
        unsigned int bins = vm[ "histogram-discrete-length" ].as<unsigned int>();
        Histogram<double> hist;
        if ( vm.count( "histogram-around-mean" ) )
          {
            double coef = vm[ "histogram-around-mean" ].as<double>();
            RegularBinner<double> rbinner
              ( statMSL.mean() - coef * sqrt( statMSL.unbiasedVariance() ),
                statMSL.mean() + coef * sqrt( statMSL.unbiasedVariance() ),
                bins );
            hist.init( rbinner );
          }
        else hist.init( bins, statMSL );
        hist.addValues( statMSL.begin(), statMSL.end() );
        hist.terminate();
        std::cout << "# Histogram of discrete length of maximal segments" << std::endl;
        std::cout << "# bins=" << bins << std::endl;
        for ( unsigned int i = 0; i < bins; ++i )
          std::cout << i << " " << hist.pdf( i ) << " " << hist.cdf( i ) << std::endl;
      }
    if ( vm.count( "histogram-euclidean-length" ) )
      {
        unsigned int bins = vm[ "histogram-euclidean-length" ].as<unsigned int>();
        Histogram<double> hist;
        if ( vm.count( "histogram-around-mean" ) )
          {
            double coef = vm[ "histogram-around-mean" ].as<double>();
            RegularBinner<double> rbinner
              ( statMSEL.mean() - coef * sqrt( statMSEL.unbiasedVariance() ),
                statMSEL.mean() + coef * sqrt( statMSEL.unbiasedVariance() ),
                bins );
            hist.init( rbinner );
          }
        else hist.init( bins, statMSEL );
        hist.addValues( statMSEL.begin(), statMSEL.end() );
        hist.terminate();
        std::cout << "# Histogram of euclidean length of maximal segments" << std::endl;
        std::cout << "# bins=" << bins << std::endl;
        for ( unsigned int i = 0; i < bins; ++i )
          std::cout << i << " " << hist.pdf( i ) << " " << hist.cdf( i ) << std::endl;
      }
  }
  return 0;
}
  

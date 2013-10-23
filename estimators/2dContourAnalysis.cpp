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

///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

using namespace DGtal;

template <typename KSpace, typename Iterator>
void analyseLengthMS( Statistic<double> & stat, Iterator itb, Iterator ite )
{
  typedef typename KSpace::Space Space;
  typedef typename Space::Point Point;
  typedef ArithmeticalDSS<Iterator,int,4> SegmentComputer;
  typedef SaturatedSegmentation<SegmentComputer> Decomposition;
  typedef typename Decomposition::SegmentComputerIterator SegmentComputerIterator;
  // Computes the tangential cover
  SegmentComputer algo;
  Decomposition theDecomposition( itb, ite, algo);
  stat.clear();
  for ( SegmentComputerIterator scIt = theDecomposition.begin(), scItEnd = theDecomposition.end();
	scIt != scItEnd; ++scIt )
    {
      const SegmentComputer & sc = *scIt;
      int64_t l = 0; 
      for ( Iterator ptIt = sc.begin(), ptItEnd = sc.end(); ptIt != ptItEnd; ++ptIt )
	++l;
      stat.addValue( (double) l );
    }
}

int main( int argc, char** argv )
{
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("freeman-chain,f", po::value<std::string>(), "Input contour as a freeman chain file name");
  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
  }
  po::notify(vm);    
  if( ( ! parseOK ) || vm.count("help") || ( argc <= 1 ) || ( ! vm.count("freeman-chain") ) )
    {
      trace.info()<< "Analyse a 2d input contour with maximal segments so as to extract some characteristics, like digitization gridstep, range of curvatures." <<std::endl << "Basic usage: "<<std::endl
      << "\t2dContourAnalysis [options] --freeman-chain <file.fc>"<<std::endl
      << general_opt << "\n";
      return 0;
    }
  
  std::string fileName = vm["freeman-chain"].as<std::string>();

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
  }
  return 0;
}
  

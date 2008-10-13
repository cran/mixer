/* GraphReader.cc
 *
 * Copyright (C) 2006 Laboratoire Statistique & Génome
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <GraphReader.h>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>

using namespace std;

namespace ermg {

  long GraphReader::DestinationAdder::_linksExamined = 0;

  GraphReader::GraphReader(const bool symmetrize, const bool enable_loop)
  {
    _enable_loop=enable_loop;
    _graph=new Graph(symmetrize);
  }

  GraphReader::~GraphReader() 
  {
  }

  void GraphReader::loadFromFile(const string &fname)
  {
    LineReader lr(this);

    ifstream gstream;
    string line;
    gstream.open(fname.c_str());
    if (!gstream.good())
      throw GraphReaderException("Unable to open file "+fname);

    while (!gstream.eof()) {
      getline(gstream,line);
      lr(line);
    }
    gstream.close();
  }

  Graph *GraphReader::getGraph()
  {
    return _graph;
  }

  GraphReader::LineReader::LineReader(GraphReader *gr) : _gr(gr), _curLine(0)
  {
  }

  void GraphReader::LineReader::operator()(const string &line) 
  {
    _curLine++;
    istringstream linestream(line);
    string fromLabel;
    linestream >> fromLabel;
    if (fromLabel.length()) {
      DestinationAdder adder(_gr,fromLabel);
      for_each(istream_iterator<string>(linestream),istream_iterator<string>(),
	       adder);
      if (adder.linksExamined()==0) {
	string msg;
	ostringstream msgstream(msg);

	msgstream << "Invalid data file : line " << _curLine << " contains no destination vertices.";
      	throw GraphReaderException(msgstream.str());
      }
	   
    }
  }


  GraphReader::DestinationAdder::DestinationAdder(GraphReader *gr, const string &fromLabel) : 
    _gr(gr), _fromLabel(fromLabel)
    {
    	_linksExamined=0;
    }
 

  long GraphReader::DestinationAdder::linksExamined()
  {
    return _linksExamined;
  }

  void GraphReader::DestinationAdder::operator()(const string &destinationLabel)
  {
 	++_linksExamined;
 	if (_gr->_enable_loop == true || 
 		(_gr->_enable_loop == false && _fromLabel != destinationLabel))
		_gr->getGraph()->addLink(_fromLabel,destinationLabel);
  }

}


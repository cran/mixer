/* GraphReader.h
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

/*!
 *   \file
 *     \author Mark Hoebeke
 *       \brief Graph Reader
 *
 *          19/05/2006
 *          
*/

#ifndef ERMG_GRAPHREADER_H
#define ERMG_GRAPHREADER_H
/*!
  \class GraphReader libermg/GraphReader.h
  \brief GraphReader class
*/


#include <exception>
#include <list>
#include <map>
#include <string>
#include <vector>


#include <Graph.h>

namespace ermg {

  class GraphReader {

  public:
    GraphReader(const bool symmetrize = false, const bool enable_loop = true);
    virtual ~GraphReader();
    void loadFromFile(const std::string &);
    Graph *getGraph();

  private:
    GraphReader(const GraphReader &);
    GraphReader &operator=(const GraphReader &);

    Graph *_graph;
    bool _enable_loop;

    class LineReader {
    public:
      LineReader(GraphReader *);
      void operator()(const std::string &);
    private:
      GraphReader *_gr;
      long _curLine;
      
    };

    class DestinationAdder {
    public:
      DestinationAdder(GraphReader *, const std::string &);
      void operator()(const std::string &);
      long linksExamined();
    private:
      GraphReader *_gr;
      std::string _fromLabel;
      static long _linksExamined;
    };
  };

  class GraphReaderException : public std::exception {
    public :
      GraphReaderException(const std::string &s) : 
      std::exception(), _msg(s) {};
    virtual ~GraphReaderException() throw() {} ;
    virtual const char * what() throw() { return _msg.c_str(); } ;
    private :
      std::string _msg;
  }; 
  

}


#endif

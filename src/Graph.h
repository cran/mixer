/* Graph.h
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
 *       \brief Graph
 *
 *          19/05/2006
 *          
*/

#ifndef ERMG_GRAPH_H
#define ERMG_GRAPH_H
/*!
  \class Graph libermg/Graph.h
  \brief Graph class
*/


#include <exception>
#include <list>
#include <map>
#include <string>
#include <vector>
#include <set>


namespace ermg {

  class Graph {

  public:
    Graph(const bool symmetrize);

    void addLink(const std::string &, const std::string &);
    bool hasLink(const std::string &, const std::string &) const;

    std::vector< std::vector <int> > getMatrix();
    std::vector< std::vector <int> > getTransposedMatrix();

    std::string getLabel(const int) const;
    int getId(const std::string &) const;

    int nbVertices() const{
      return _labelmap.size();
    }

    //! shuffle the vertices order
    void shuffle();
    //! reorder the vertices from classified vertices and return their classes 
    std::vector <int> reorder(const std::string &);

    virtual ~Graph();

  private:

    mutable std::map<std::string,int> _labelmap;
    mutable std::map<int,std::string> _idmap;
    mutable std::map< int, std::list <int> > _linkmap;

    int _nextid;
    bool _symmetrize;
  };
}


#endif

/* ModelDesigner.h
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
 *     \author Vincent Miele
 *       \brief Erdös Reni Mixture for Graphs design
 *
 *          13/10/2006
 *          
*/

#ifndef ERMG_MODELDESIGNER_H
#define ERMG_MODELDESIGNER_H
/*!
  \class ModelDesigner.h libermg/ModelDesigner.h
  \brief Erdös Reni Mixture for Graphs design  class
*/


#include <cmath>
#include <Emd.h>
#include <Ermdg.h>

namespace ermg {
  
  class ModelDesigner
  {   
  protected:     
    //! number of classes
    int _q;
    //! silent mode
    bool _silent;
    
    //! curent Ermg
    Ermg& _curr_ermg;
  public:
    ModelDesigner( Ermg& ermg, int q, bool silent );
    void load(const std::vector<int>& inputclasses);
    void initialRandom(int n_init);
    void em(EmCore& em);
    void outFile(const std::string& ofile, const Graph *g, const std::string& odir );
  };
  
}
#endif

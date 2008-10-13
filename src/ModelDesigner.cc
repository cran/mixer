/* ModelDesigner.cc
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

#include <ModelDesigner.h>
using namespace std;


namespace ermg {

  ModelDesigner::ModelDesigner( Ermg& ermg, 
				int q, bool silent )
    : _curr_ermg(ermg)
  {
    _q = q;
    _silent = silent;
  }

  void ModelDesigner::load(const vector<int>& inputclasses)
  {    	      
    if (!_silent)
      cout<<"Initializing classes...\t"<< flush; 
    _curr_ermg.replaceCahForSubXadjByLoad(inputclasses);      
    if (!_silent)
      cout<<"done"<<endl;
  }
  
  void ModelDesigner::initialRandom(int n_init)
  {   	      
    if (!_silent)
      cout<<"Random "<< flush;
    vector<int> randomclasses(_curr_ermg.nbVertices(), -1);
    for (int i=0; i<n_init; i++)
      randomclasses[i] = int(double(rand())/double(RAND_MAX)*_q);
      
    load(randomclasses);
  }

  void ModelDesigner::em(EmCore& em)
  {      
    if (!_silent)
      cout<<"\rEM for Q="<<_q<<"...\t"<< flush;
    _curr_ermg.em(em, _q);      
    if (!_silent)
      cout<<"done"<<endl;
  }
  
  void ModelDesigner::outFile(const string& ofile, const Graph *g, const string& odir )
  {
    if (!_silent)
      cout<<"Saving...\t"<< flush; 
      _curr_ermg.outFile( ofile, g, odir );
      
    stringstream strout;
    if (odir.length()>0){
      string::size_type islash = ofile.find_last_of("/");
      strout<<odir<<"/"<<ofile.substr(islash+1, ofile.size()-4-islash-1)<<".likelihoods";
    }
    else{
      strout<<ofile.substr(0, ofile.size()-4)<<".likelihoods";
    }
    ofstream fout;
    fout.open(strout.str().c_str(), ios::out);
    fout<<"# Q="<<_q
	<<" Complete Likelihood:\n"<<_curr_ermg.completeLikelihood()<<endl;  
    fout<<endl<<endl;
    fout<<"# Q="<<_q
	<<" Entropy:\n"<<_curr_ermg.entropy()<<endl;    
    fout<<endl<<endl;
    fout<<"# Q="<<_q
	<<" ICL:\n"<<_curr_ermg.ICL()<<endl;    
    fout.close(); 
    if (!_silent)
      cout<<"done"<<endl; 
  }


}

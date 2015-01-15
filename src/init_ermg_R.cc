/* ermg_R.cc
 *
 * Copyright (C) 2008 Laboratoire Statistique & G�nome
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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>
using namespace std;


#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <ModelImprover.h>
#include <Emd.h>
#include <GraphReader.h>

#include <R.h>
#include <Rmath.h>



using namespace ermg;

extern "C" {
  void init_ermg(int*    symmetrize, // param for ERMG initialzation
		 int*    loop,
		 int*    undirected,
		 int*    silent,
		 int*    kmeansnbclass,
		 int*    kmeansnbiter,
		 int*    nokmeans,
                 int*    q_p,
		 int*    nbrEdges,	// nbr of edges, in the given graph
		 int*    nbrNodes,
		 int*    m,	// the list of edges
		 double* tau) {

    GetRNGstate(); // Initialize random number generator
    
    int                 q = *q_p; 
    int     improvenbiter = 0;
    int     emnbiter      = 0 ;
    double  emeps         = 1e5;
    int     fpnbiter      = 0 ;
    double  fpeps         = 1e5 ;
    // int     improve       = 1;
    int     classif       = 1;
    int     qmin          = q;
    int     qmax          = q;

    // GRAPH
    bool tosym = *symmetrize;
    if (*undirected){
      tosym = true;
    }
    GraphReader gr(tosym, *loop);
    

    // loading graph from *m
    // need to convert int to string
    char buf1[50];
    char buf2[50];
    int i,n;
    
    try {
      Graph *graph=gr.getGraph();
      for (i=0; i<*nbrEdges*2; i += 2){
	n=sprintf(buf1, "%d", m[i]);
	n=sprintf(buf2, "%d", m[i+1]);
	graph->addLink(buf1, buf2);
      }

      
      if (graph->nbVertices() < qmax){

      }    
      
      // MODEL 
      Ermg* model;
      if (*undirected == false){
	model = new Ermdg( graph->getMatrix(), graph->getTransposedMatrix(),
			   *loop,
			   *kmeansnbiter,
			   *kmeansnbclass,
			   qmin, qmax );
      }
      else{	
	model = new Ermg( graph->getMatrix(),
			  *loop,
			  *kmeansnbiter,
			  *kmeansnbclass,
			  qmin, qmax );
      } 
      
      // METHOD 
      EmCore * em_method = NULL;
      if (*undirected ==false){
	// EM UNDIRECTED
	em_method = new Emd(  dynamic_cast<Ermdg&>(*model),
			      emnbiter, emeps,
			      fpnbiter, fpeps );
      }
      else{
	// EM DIRECTED
	em_method = new Em(  *model,
			     emnbiter, emeps,
			     fpnbiter, fpeps, classif);
      }    
      
      
      // ESTIMATION
      ModelImprover modelimp( *model, qmin, qmax, improvenbiter, *silent );
      modelimp.initialLikelihoods( *em_method, *nokmeans );
      // ??? if (*improve == true){
      //	modelimp.improveLikelihoods( *em_method );
      //}
      
      
      // double icl ;
      std::vector<double>::iterator it_Alpha ;
      BandMemMatrix<double>::BandCursor it_Pi ;
      BandMemMatrix<double>::BandCursor it_t_Tau ;

      //int q=0;
      int i=0;
      int j=0, k=0;
      //      for(q=*qmin; q<=*qmax; q++){
	// icl = modelimp.getModelForNbclass(q).ICL();
	// res[i++] = icl;
      
	// int limit = modelimp.getModelForNbclass(q)._Alpha.size();
	// for (j=0; j < limit; j++)
	//   res[i++] = modelimp.getModelForNbclass(q)._Alpha[j];

	// int lines = modelimp.getModelForNbclass(q)._Pi.nLines();
	// int cols = modelimp.getModelForNbclass(q)._Pi.nCols();
	// for (j=0; j < lines; j++)
	//  for(k=0; k < cols; k++)
	//    res[i++] = modelimp.getModelForNbclass(q)._Pi.value(j,k);
	
	int lines = modelimp.getModelForNbclass(q)._t_Tau.nLines(); 
	int cols = modelimp.getModelForNbclass(q)._t_Tau.nCols();
	for (k=0; k < lines; k++) {
	  for (j=0; j < cols; j++) {
	    tau[i++] = modelimp.getModelForNbclass(q)._t_Tau.value(k,j);
	  }  
	}

      delete model;
      delete graph;
      delete em_method;
      
    } catch (GraphReaderException &gre) {

    }
     PutRNGstate(); // send RNG status to R
  }
 
}
 

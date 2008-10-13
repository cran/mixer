/* ermg_R.cc
 *
 * Copyright (C) 2008 Laboratoire Statistique & Génome
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


using namespace ermg;

extern "C" {
  void main_ermg(int*    symmetrize, // param for main_ermg function
		 int*    loop,
		 int*    undirected,
		 int*    silent,
		 int*    improvenbiter,
		 int*    kmeansnbclass,
		 int*    kmeansnbiter,
		 int*    emnbiter,
		 double* emeps,
		 int*    fpnbiter,
		 double* fpeps,
		 int*    improve,
		 int*    classif,
		 int*    nokmeans,
		 int*    qmax,
		 int*    qmin,
		 int*    nbrEdges,	// nbr of edges, in the given graph
		 int*    nodes, // nbr of nodes of the graph
		 int*    m,	// the list of edges
		 double* res){
    
    srand((long)getpid());
    
    // GRAPH
    bool tosym = *symmetrize;
    if (*undirected){
      tosym = true;
    }
    GraphReader gr(tosym, *loop);
    if (!(*silent))
      printf("Loading graph...\t");
    

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
      if (!(*silent))
	printf("done");
      
      if (graph->nbVertices() < *qmax){
	printf("Error: Q max is greater than (then set to) the number of vertices");
       	exit(0);
      }    
      
      
      // MODEL 
      Ermg* model;
      if (*undirected == false){
	model = new Ermdg( graph->getMatrix(), graph->getTransposedMatrix(),
			   *loop,
			   *kmeansnbiter,
			   *kmeansnbclass,
			   *qmin, *qmax );
      }
      else{	
	model = new Ermg( graph->getMatrix(),
			  *loop,
			  *kmeansnbiter,
			  *kmeansnbclass,
			  *qmin, *qmax );
      } 
      
      // METHOD 
      EmCore * em_method = NULL;
      if (*undirected ==false){
	// EM UNDIRECTED
	em_method = new Emd(  dynamic_cast<Ermdg&>(*model),
			      *emnbiter, *emeps,
			      *fpnbiter, *fpeps );
      }
      else{
	// EM DIRECTED
	em_method = new Em(  *model,
			     *emnbiter, *emeps,
			     *fpnbiter, *fpeps, *classif);
      }    
      
      
      // ESTIMATION
      ModelImprover modelimp( *model, *qmin, *qmax, *improvenbiter, *silent );
      modelimp.initialLikelihoods( *em_method, *nokmeans );
      if (*improve == true){
	modelimp.improveLikelihoods( *em_method );
      }
      
      
      double icl ;
      std::vector<double>::iterator it_Alpha ;
      BandMemMatrix<double>::BandCursor it_Pi ;
      BandMemMatrix<double>::BandCursor it_t_Tau ;

      int q=0;
      int i=0;
      int j=0, k=0;
      for(q=*qmin; q<=*qmax; q++){
	icl = modelimp.getModelForNbclass(q).ICL();
	res[i++] = icl;
      
	int limit = modelimp.getModelForNbclass(q)._Alpha.size();
	for (j=0; j < limit; j++)
	  res[i++] = modelimp.getModelForNbclass(q)._Alpha[j];

	int lines = modelimp.getModelForNbclass(q)._Pi.nLines();
	int cols = modelimp.getModelForNbclass(q)._Pi.nCols();
	for (j=0; j < lines; j++)
	  for(k=0; k < cols; k++)
	    res[i++] = modelimp.getModelForNbclass(q)._Pi.value(j,k);
	
	lines = modelimp.getModelForNbclass(q)._t_Tau.nLines(); 
	cols = modelimp.getModelForNbclass(q)._t_Tau.nCols();
	for (k=0; k < lines; k++)
	  for (j=0; j < cols; j++)
	    res[i++] = modelimp.getModelForNbclass(q)._t_Tau.value(k,j);

      }

      delete model;
     delete graph;
           delete em_method;
      
    } catch (GraphReaderException &gre) {
      printf("Exception : %s\n",gre.what());
    }
  }
}
 

/* ermgo_R.cc
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


#include <sys/types.h>
#include <unistd.h>
#include <ModelImprover.h>
#include <OCEm.h>
#include <OEm.h>
#include <SOCEm.h>
#include <GraphReader.h>


using namespace ermg;

 
extern "C" {
  void main_ermgo(int*    loop,
		  int*    silent,
		  int*    initnbv,
		  int*    improvenbiter,
		  int*    emnbiter,
		  int*    improve,
		  int*    classif,
		  int*    stochastique,
		  int*    qmax,
		  int*    qmin,
		  int*    nbrEdges,	// nbr of edges, in the given graph
		  int*    nodes, // nbr of nodes of the graph
		  int*    m,	// the list of edges
		  double* res){

    srand((long)getpid());
    
    // GRAPH
    bool tosym = true;
    GraphReader gr(tosym, *loop);
    if (!*silent)
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
      
      //graph->shuffle();
      
      
      // MODEL 
      Ermg* model;
      model = new Ermg( graph->getMatrix(),
			*loop,
			0,
			0,
			*qmin, *qmax );
      
      
      // METHOD 
      EmCore * em_method = NULL;
      int nbv = *initnbv;
      if (nbv==0)
	nbv = model->nbVertices()/3;
      if (nbv<MININITNBV)
	nbv = MININITNBV;
      if (nbv>model->nbVertices())
	nbv = model->nbVertices();

      if (*classif==true) 
	// CLASSIFFICATION ONLINE EM 
	em_method = new OCEm( *model, nbv, *emnbiter );	  
      else 
	if (*stochastique == true) 
	  // STOCHASTIC CLASSIFICATION ONLINE EM 
	  em_method = new SOCEm( *model, nbv, *emnbiter);	  
	else 
	  // ONLINE EM 
	  em_method = new OEm( *model, nbv, *emnbiter );	  
      
      
      // ESTIMATION
      ModelImprover modelimp( *model, *qmin, *qmax, *improvenbiter, *silent );
      modelimp.initialLikelihoods( *em_method );
      
      if ( *improve == true){
	modelimp.improveLikelihoods( *em_method );
	//modelimp.improveLikelihoodsBySelection( *em_method );
      }
      
      //     double cpl = modelimp.getModelForNbclass(*qmin)._complete_likelihood;
      
      //     std::vector<double>::iterator it_Alpha = modelimp.getModelForNbclass(*qmin)._Alpha.begin();
      //     double Alpha_0 = modelimp.getModelForNbclass(*qmin)._Alpha[0];
      //     modelimp.getModelForNbclass(*qmin)._Alpha.size();
      
      //     BandMemMatrix<double>::BandCursor it_Pi = modelimp.getModelForNbclass(*qmin)._Pi.begin();
      //     double Pi_00 = modelimp.getModelForNbclass(*qmin)._Pi.value(0,0);
      //     modelimp.getModelForNbclass(*qmin)._Pi.nLines();
      //     modelimp.getModelForNbclass(*qmin)._Pi.nCols();
      
      //     BandMemMatrix<double>::BandCursor it_t_Tau = modelimp.getModelForNbclass(*qmin)._t_Tau.begin();
      //     double t_Tau_00 = modelimp.getModelForNbclass(*qmin)._t_Tau.value(0,0);
      //     modelimp.getModelForNbclass(*qmin)._t_Tau.nLines();
      //     modelimp.getModelForNbclass(*qmin)._t_Tau.nCols();
      
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
      printf("Exception : %s\n", gre.what());
    }
  }

}

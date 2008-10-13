/* Graph.cc
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

#include <Graph.h>

#include <cassert>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>

using namespace std;

namespace ermg {

  Graph::Graph(const bool symmetrize) : _nextid(0), _symmetrize(symmetrize) 
  {
  }

  Graph::~Graph() {}

  vector<vector<int> >  Graph::getMatrix()
  {
    vector<vector<int> > res;
    //res.resize(_linkmap.size());
    res.resize(_idmap.size(), vector<int>());
    map<int,list<int> >::iterator idit=_linkmap.begin();
    int previd = 0;
    while (idit != _linkmap.end()) {
      previd = idit->first;
      assert(idit->second.size()>=0);
      /*
	cerr << " Processsing " << idit->first << "\t( " << _idmap[idit->first] 
	<< ")\t size = " << idit->second.size() << endl;
      */
      res[idit->first].reserve(idit->second.size());
      idit->second.sort();
      unique_copy(idit->second.begin(),idit->second.end(),
		  back_inserter(res[idit->first]));
      
      idit++;
    }

    return res;
  }

  vector<vector<int> >  Graph::getTransposedMatrix()
  {
    vector<vector<int> > res;
    //res.resize(_linkmap.size());
    res.resize(_idmap.size(), vector<int>());

    map<int,list<int> > tmplinkmap=_linkmap;
    map<int,list<int> >::iterator idit=tmplinkmap.begin();
    while (idit != tmplinkmap.end()) {
      idit->second.sort();
      idit->second.unique();

      idit++;
    }

    idit=tmplinkmap.begin();
    while (idit != tmplinkmap.end()) {
      assert(idit->second.size()>=0);      
      /*
	cerr << " Processsing " << idit->first << "\t( " << _idmap[idit->first] 
	<< ")\t size = " << idit->second.size() << endl;
      */
      for (list<int>::iterator listit=idit->second.begin(); listit!=idit->second.end(); listit++){      
	res[*listit].push_back(idit->first);
      }      
      idit++;
    }
    return res;
  }

  string Graph::getLabel(const int id) const
  {
    return _idmap[id];
  }

  int Graph::getId(const string &label) const 
  {
    return _labelmap[label];
  }

  void Graph::addLink(const string & fromLabel, const string &destinationLabel)
  {
    int fromId=-1;
    map<string,int>::iterator fit=_labelmap.find(fromLabel);
    if (fit != _labelmap.end()) {
      fromId=fit->second;
    } else {
      fromId=_nextid++;
      _labelmap[fromLabel]=fromId;
      _idmap[fromId]=fromLabel;
    }


    int toId=-1;
    map<string,int>::iterator toit=_labelmap.find(destinationLabel);
    if (toit != _labelmap.end()) {
      toId=toit->second;
    } else {
      toId=_nextid++;
      _labelmap[destinationLabel]=toId;
      _idmap[toId]=destinationLabel;
    }

    _linkmap[fromId].push_back(toId);
    
    if ((_symmetrize)&&(fromId!=toId))
      _linkmap[toId].push_back(fromId);

  }

  bool Graph::hasLink(const string &fromLabel, const string &toLabel) const
  {
    bool res=false;

    int fromId=_labelmap[fromLabel];
    int toId=_labelmap[toLabel];

    list<int> fromLinks=_linkmap[fromId];
    
    list<int>::iterator it=find(fromLinks.begin(),fromLinks.end(),toId);
    if (it != fromLinks.end())
      res=true;

    return res;
    
  }

  void Graph::shuffle()
  {
    vector<int> shuffle_vec(_idmap.size());
    for (int index=0;index<int(shuffle_vec.size());index++)
      shuffle_vec[index]=index;

    random_shuffle(shuffle_vec.begin(),shuffle_vec.end());

    map<string,int> shuffled_labelmap;
    for (map<string,int>::iterator it=_labelmap.begin();it!=_labelmap.end(); it++) {
      int shuffled_id=shuffle_vec[it->second];
      shuffled_labelmap[it->first]=shuffled_id;
    }
    _labelmap=shuffled_labelmap;

    map<int,string> shuffled_idmap;
    for (map<int,string>::iterator it=_idmap.begin();it!=_idmap.end(); it++) {
      int shuffled_id=shuffle_vec[it->first];
      shuffled_idmap[shuffled_id]=it->second;
    }
    _idmap=shuffled_idmap;

    
    map<int, list<int> > shuffled_linkmap;
    for (int i=0;i<int(shuffle_vec.size());i++) {
      int shuffled_i=shuffle_vec[i];
      shuffled_linkmap[shuffled_i]=_linkmap[i];
      for (list<int>::iterator it=shuffled_linkmap[shuffled_i].begin();
	   it!=shuffled_linkmap[shuffled_i].end();
	   it++) {
	int id=*it;
	*it=shuffle_vec[id];
      }
	
    }
    _linkmap=shuffled_linkmap;
  } 


  vector<int> Graph::reorder(const string& vfile)
  {
    map<string, int> classlabelmap;
    ifstream gstream;
    string line;
    gstream.open(vfile.c_str());
    if (!gstream.good()){
      cerr<<"Unable to open file "<<vfile<<endl;
      exit(1);
    }
    while (!gstream.eof()) {
      getline(gstream,line); 
      istringstream linestream(line);
      string tmplabel;
      linestream>>tmplabel; 
      if (tmplabel.length()) {
	string cl;
	linestream>>cl;
	if (!cl.length()){
	  cerr<<"Error: the line format \"id class\" is not respected in "<<vfile<<endl;
	  exit(1);
	}
	int icl=atoi(cl.c_str());
	classlabelmap[tmplabel]=icl;
      }
    }
    gstream.close();
    vector<int> vclass;
    vclass.reserve(classlabelmap.size());
    for (map<string, int>::iterator it=classlabelmap.begin(); 
	 it!=classlabelmap.end(); it++)
      vclass.push_back(it->second);


    vector<int> rearranged_vec(_idmap.size());
    int nbs = classlabelmap.size(), i=0, j=nbs;
    map<string,int> rearranged_labelmap;
    for (map<string,int>::iterator it=_labelmap.begin();it!=_labelmap.end(); it++) {
      if (classlabelmap.find(it->first)!=classlabelmap.end()){
	rearranged_labelmap[it->first]=i;
	rearranged_vec[it->second]=i;
// 	cerr<<"rearranged_labelmap["<<it->first<<"]="<<i<<endl;
// 	cerr<<"rearranged_vec["<<it->second<<"]="<<i<<endl;
	i++;
      }
      else{
	rearranged_labelmap[it->first]=j;
	rearranged_vec[it->second]=j;
// 	cerr<<"rearranged_labelmap["<<it->first<<"]="<<j<<endl;
// 	cerr<<"rearranged_vec["<<it->second<<"]="<<j<<endl;
	j++;
      }
    }
    _labelmap=rearranged_labelmap;


    map<int,string> rearranged_idmap;
    for (map<int,string>::iterator it=_idmap.begin();it!=_idmap.end(); it++) {
      int rearranged_id=rearranged_vec[it->first];
      rearranged_idmap[rearranged_id]=it->second;
    }
    _idmap=rearranged_idmap;

    
    map<int, list<int> > rearranged_linkmap;
    for (int i=0;i<int(rearranged_vec.size());i++) {
      int rearranged_i=rearranged_vec[i];
      rearranged_linkmap[rearranged_i]=_linkmap[i];
      for (list<int>::iterator it=rearranged_linkmap[rearranged_i].begin();
	   it!=rearranged_linkmap[rearranged_i].end();
	   it++) {
	int id=*it;
	*it=rearranged_vec[id];
      }
	
    }
    _linkmap=rearranged_linkmap;

//     cerr<<"AFTER"<<endl;
//     for (int i=0;i<int(_linkmap.size());i++){
//       cerr<<i<<" ";
//       for (list<int>::iterator it=_linkmap[i].begin();
// 	   it!=_linkmap[i].end();
// 	   it++){
// 	cerr<<*it<<" ";
//       }
//       cerr<<endl;
//     }
      
      

    return vclass;
  }
 
}


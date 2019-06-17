/****************************************************************************
* JMeshExt                                                                  *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
*                                                                           *
* Copyright(C) 2006: IMATI-GE / CNR                                         *
*                                                                           *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#include "exttrimesh.h"
#include <stdio.h>
#include <stdlib.h>


Point laplacianDisplacement(Vertex *v)
{
 List *vv = v->VV();
 Vertex *w;
 Node *m;
 Point np;
 FOREACHVVVERTEX(vv, w, m) np = np+(*w);
 np = (np)/(vv->numels());
 delete(vv);
 return np;
}

Point laplacianDisplacementOrder2(Vertex *v)
{
 List *vv = v->VV();
 List *vl2 = v->VV();
 Vertex *w, *w2;
 Node *m, *m2;
 Point np;
 double total_weight = 0.0;
 double weight;
 FOREACHVVVERTEX(vv, w, m){
   FOREACHVVVERTEX(w->VV(), w2, m2){
     if(vl2->containsNode(w2) == NULL && w2 != v){
       vl2->appendTail(w2);
     }
   }
 }

 FOREACHVVVERTEX(vl2, w, m){
   weight = 1. / (v->distance(*w) + 1e-5);
   np = np + ((*w) * weight);
   total_weight = total_weight + weight;
   //np = np + (*w);
 }
 //np = (np)/(vl2->numels());
 np = (np)/total_weight;
 delete(vv);
 delete(vl2);
 return np;
}


Point sharpLaplacianDisplacement(Vertex *v)
{
 List *ve = v->VE();
 Vertex *w;
 Node *m;
 Edge *e;
 Point np;
 int nse=0;

 FOREACHVEEDGE(ve, e, m)
  if (IS_SHARPEDGE(e) || e->isOnBoundary())
  {
   if (nse==0) np.setValue(e->oppositeVertex(v));
   else if (nse==1) np = np + (*(e->oppositeVertex(v)));
   else {delete(ve); return (*v);}
   nse++;
  }
  else if (!nse) {w = e->oppositeVertex(v); np = np+(*w);}

 if (!nse) np = (np)/(ve->numels());
 else if (nse == 1) np = (*v);
 else np = np/2;

 delete(ve);
 return np;
}

int ExtTriMesh::laplacianSmooth(int ns, double l)
{
 Triangle *t;
 Edge *e;
 Vertex *v;
 Node *n;
 int i = 0, is_selection = 0, ins = ns;
 double ln = 1.0-l;
 Point np;

 FOREACHTRIANGLE(t, n) if (IS_VISITED(t))
  {MARK_VISIT(t->e1); MARK_VISIT(t->e2); MARK_VISIT(t->e3);}
 FOREACHEDGE(e, n) if (IS_VISITED(e))
  {MARK_VISIT(e->v1); MARK_VISIT(e->v2); is_selection = 1;}

 List vts;
 FOREACHVERTEX(v, n) if (!is_selection || IS_VISITED(v)) vts.appendHead(v);

 coord *xyz = (coord *)malloc(sizeof(coord)*vts.numels()*3);
 if (xyz == NULL) {JMesh::warning("Not enough memory for vertex coordinates.\n"); return 0;}

 JMesh::begin_progress();
 for (; ns>0; ns--)
 {
  i=0;
  FOREACHVVVERTEX((&vts), v, n)
  {
   np = sharpLaplacianDisplacement(v);
   if (!(i%3000)) JMesh::report_progress("%d %% done - %d steps left",((i*33)/(vts.numels()) + (100*(ins-ns)))/ins, ns);
   xyz[i++] = np.x*l+v->x*ln; xyz[i++] = np.y*l+v->y*ln; xyz[i++] = np.z*l+v->z*ln; 
  }

  i=0;
  FOREACHVVVERTEX((&vts), v, n)
  {
   v->x = xyz[i++]; v->y = xyz[i++]; v->z = xyz[i++];
  }
 }
 JMesh::end_progress();
 free(xyz);

 FOREACHEDGE(e, n) UNMARK_VISIT(e);
 FOREACHVERTEX(v, n) UNMARK_VISIT(v);

 return 1;
}

int ExtTriMesh::taubinSmooth(int ns, double l)
{
 Vertex *v;
 Node *n;
 int i = 0, is_selection = 0, ins = ns;
 double ln = 1.0-l;
 double m = -1.02*l;
 double mn = 1.0-m;
 Point np;
 //printf("Lambda: %f", l);
 JMesh::info("Lambda: %f\n", l);
 List vts;
 FOREACHVERTEX(v, n) vts.appendHead(v);
 coord *xyz = (coord *)malloc(sizeof(coord)*vts.numels()*3);
 if (xyz == NULL) {JMesh::warning("Not enough memory for vertex coordinates.\n"); return 0;}

 JMesh::begin_progress();
 for (; ns>0; ns--)
 {
  i=0;
  FOREACHVVVERTEX((&vts), v, n)
  {
   //np = laplacianDisplacement(v);
   np = laplacianDisplacementOrder2(v);
   if (!(i%3000)) JMesh::report_progress("%d %% done - %d steps left",((i*33)/(vts.numels()) + (100*(ins-ns)))/ins, ns);
   xyz[i++] = np.x*l+v->x*ln; xyz[i++] = np.y*l+v->y*ln; xyz[i++] = np.z*l+v->z*ln; 
  }

  i=0;
  FOREACHVVVERTEX((&vts), v, n)
  {
   v->x = xyz[i++]; v->y = xyz[i++]; v->z = xyz[i++];
  }

  i=0;
  FOREACHVVVERTEX((&vts), v, n)
  {
   //np = laplacianDisplacement(v);
   np = laplacianDisplacementOrder2(v);
   xyz[i++] = np.x*m+v->x*mn; xyz[i++] = np.y*m+v->y*mn; xyz[i++] = np.z*m+v->z*mn; 
  }

  i=0;
  FOREACHVVVERTEX((&vts), v, n)
  {
   v->x = xyz[i++]; v->y = xyz[i++]; v->z = xyz[i++];
  }

 }
 JMesh::end_progress();
 free(xyz);

 return 1;
}

int ExtTriMesh::removeSpikes(double threshold)
{
 Triangle *t;
 Edge *e;
 Vertex *v;
 Node *n, *n2;
 List vts;
 FOREACHTRIANGLE(t, n) UNMARK_VISIT(t);
 FOREACHVERTEX(v, n)
 {
      if (v->totalAngle() <= threshold)
      {
          FOREACHVTTRIANGLE(v->VT(), t, n2) MARK_VISIT(t);
      }
 }
 this->removeSelectedTriangles();
 this->unmarkEverything();

 return 1;
}


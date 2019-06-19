#include "exttrimesh.h"
#include "component.h"
#include "detectIntersections.h"

/* Assumes that the Triangulation consists of exactly 2 components, each having no selfintersections.
   If they overlap, they will be joined and the overlapping parts will be deleted. */
int ExtTriMesh::joinOverlappingComponentPair() {
    this->deselectTriangles();
    List *components = this->getComponentsList();
    if( components->numels() > 2 ) JMesh::error("Triangulation consists of more than 2 components.\n");
    if( components->numels() < 2 ) { JMesh::info("Only 1 component, nothing joined.\n"); return 0; }
    // select the intersecting triangles, which form the boundaries
    if(!this->selectIntersectingTriangles()) {
        JMesh::info("Component pair doesn't overlap. Nothing joined.\n");
        return 0;
    }
    // Identify the two most distant triangles (which are assumed to be not part of the overlap),
    // to determine which chunks to keep after removing the intersecting triangles
    ComponentStruct c1((List*) components->popHead()), c2((List*) components->popHead());
    Triangle *t1, *t2;
    Node *n;
    FOREACHVTTRIANGLE(c1.triangles, t1, n) if(IS_VISITED(t1)) break;
    FOREACHVTTRIANGLE(c2.triangles, t2, n) if(IS_VISITED(t2)) break;
    c1.vertices = c1.getVertices(2);
    c2.vertices = c2.getVertices(2);
    Vertex *v1, *v2;
    this->mostDistantPartner(t1->v1(), c2.vertices, &v2);
    this->mostDistantPartner(t2->v1(), c1.vertices, &v1);
    c1.clear();
    c2.clear();
    t1 = v1->e0->t1;
    t2 = v2->e0->t1;
    // now delete the boundaries and have at least 2 shells with boundaries
    this->removeSelectedTriangles();
    if(!t1 || !t2) JMesh::error("Algorithm using most distant points didn't work ...\n");
    this->selectConnectedComponent(t1, false);
    this->selectConnectedComponent(t2, false);
    this->invertSelection();
    this->removeSelectedTriangles();
    // remove more triangles close to the overlap, to make the holes bigger for better joining
    this->selectBoundaryTriangles();
    this->growSelection();
    this->removeSelectedTriangles();
    this->removeSmallestComponents(2);
    c1 = ComponentStruct(t1);
    c2 = ComponentStruct(t2);
    int ret = this->joinComponentsBiggestBoundaryPair(c1.triangles, c2.triangles, DBL_MAX);
    this->eulerUpdate();
    this->fillSmallBoundaries(0, true, false);
    return ret;
}

/* Assumes that the Triangulation consists of exactly 2 components, each having no selfintersections.
   If they overlap, they will be joined and the overlapping parts will be deleted. */
int ExtTriMesh::joinOverlappingComponentPair2(bool firstIn, bool secondIn) {
    this->deselectTriangles();  
    List *components = this->getComponentsList(); // this adds new entries to the head --> i.e. second component is first here
    if( components->numels() > 2 ) JMesh::error("Triangulation consists of more than 2 components.\n");
    if( components->numels() < 2 ) { JMesh::info("Only 1 component, nothing joined.\n"); return 0; }
    // select the intersecting triangles, which form the boundaries
    List *first = (List*) components->head()->data;
    Triangle *t; Node *n;
    // mark4 first component (in fact, this is the 2nd mesh on the command line, as getComponentsList adds to the head of list)
    FOREACHVTTRIANGLE(first, t, n) {
        MARK_BIT(t,4);
        MARK_BIT(t->v1(),4);
        MARK_BIT(t->v2(),4);
        MARK_BIT(t->v3(),4);
    }
    delete(first);
    delete(components);
    // mark5 second component
    FOREACHTRIANGLE(t, n) if(!IS_BIT(t,4)) {
        MARK_BIT(t,5);
        MARK_BIT(t->v1(),5);
        MARK_BIT(t->v2(),5);
        MARK_BIT(t->v3(),5);
    }

    this->markTrianglesInsideComponent(7, 5, 4, firstIn); // second item in component list referts to first item on command line
    this->markTrianglesInsideComponent(6, 4, 5, secondIn);

    JMesh::info("first %d; second %d\n",firstIn,secondIn);

    // unmark mask bits that are used by other subfunctions to prevent side effects
    FOREACHTRIANGLE(t, n) { UNMARK_BIT(t, 1); UNMARK_BIT(t, 2); UNMARK_BIT(t, 3); }
    // unmark bits 4 and 5 when inside
    FOREACHTRIANGLE(t, n) if(t->mask & (1<<6 | 1<<7)) { UNMARK_BIT(t, 4); UNMARK_BIT(t, 5); }

    // mark parts to remove
    if(firstIn) { FOREACHTRIANGLE(t, n) if(t->mask & (1<<5)) MARK_VISIT(t); }
    else {FOREACHTRIANGLE(t, n) if(t->mask & (1<<7)) MARK_VISIT(t); }
    if(secondIn) { FOREACHTRIANGLE(t, n) if(t->mask & (1<<4)) MARK_VISIT(t); }
    else {FOREACHTRIANGLE(t, n) if(t->mask & (1<<6)) MARK_VISIT(t); }

    this->removeSelectedTriangles();
    this->eulerUpdate();
    this->removeSmallestComponents(2);

    if(this->shells() != 2) return 0;

    // invert orientation of second component if necessary
    Edge *e;
    FOREACHEDGE(e, n) UNMARK_VISIT(e) ;
    if(secondIn && (secondIn != firstIn)) {
        JMesh::info("inverting triangle normals and edge directions of second component\n");
        FOREACHTRIANGLE(t, n) if((t->mask & (1<<6))) {
            t->invert();
            if(!IS_VISITED(t->e1)) { MARK_VISIT(t->e1);  p_swap((void **)(&(t->e1->v1)), (void **)(&(t->e1->v2))); };
            if(!IS_VISITED(t->e2)) { MARK_VISIT(t->e2);  p_swap((void **)(&(t->e2->v1)), (void **)(&(t->e2->v2))); };
            if(!IS_VISITED(t->e3)) { MARK_VISIT(t->e3);  p_swap((void **)(&(t->e3->v1)), (void **)(&(t->e3->v2))); };
        }
    }

    FOREACHTRIANGLE(t, n) { if(!(IS_BIT(t,4)||IS_BIT(t,6))) break; }
    ComponentStruct c1(t);
    FOREACHTRIANGLE(t, n) { if(!(IS_BIT(t,5)||IS_BIT(t,7))) break; }
    ComponentStruct c2(t);

    int ret = this->joinComponentsBiggestBoundaryPair(c1.triangles, c2.triangles, 10);
    if(!ret) {
        this->selectBoundaryTriangles();
        this->removeSelectedTriangles();

        FOREACHTRIANGLE(t, n) { if(!(IS_BIT(t,4)||IS_BIT(t,6))) break; }
        ComponentStruct c1(t);
        FOREACHTRIANGLE(t, n) { if(!(IS_BIT(t,5)||IS_BIT(t,7))) break; }
        ComponentStruct c2(t);

        ret = this->joinComponentsBiggestBoundaryPair(c1.triangles, c2.triangles, 10);
        if(!ret) return 0;
    }

    this->eulerUpdate();
    this->fillSmallBoundaries(0, true, false);
    this->unmarkEverything();
    return ret;
}


/* Assumes that the Triangulation consists of exactly 2 components, each having no selfintersections.
   If they overlap, the outer part of the second component will be attached to the first one; the first one is kept unchanged*/
/*
int ExtTriMesh::attachTo() {
    this->deselectTriangles();
    List *components = this->getComponentsList(); // this adds new entries to the head --> i.e. second component is first here
    if( components->numels() > 2 ) JMesh::error("Triangulation consists of more than 2 components.\n");
    if( components->numels() < 2 ) { JMesh::info("Only 1 component, nothing joined.\n"); return 0; }
    // select the intersecting triangles, which form the boundaries
    List *first = (List*) components->head()->data;
    Triangle *t; Node *n; Vertex *v;
    // mark4 first component (in fact, this is the 2nd mesh on the command line, as getComponentsList adds to the head of list)
    FOREACHVTTRIANGLE(first, t, n) {
        MARK_BIT(t,4);
        MARK_BIT(t->v1(),4);
        MARK_BIT(t->v2(),4);
        MARK_BIT(t->v3(),4);
    }
    delete(first);
    delete(components);
    // mark5 second component
    FOREACHTRIANGLE(t, n) if(!IS_BIT(t,4)) {
        MARK_BIT(t,5);
        MARK_BIT(t->v1(),5);
        MARK_BIT(t->v2(),5);
        MARK_BIT(t->v3(),5);
    }

    this->markTrianglesInsideComponent(7, 5, 4, 0); // second item in component list referts to first item on command line
    this->markTrianglesInsideComponent(6, 4, 5, 0);

    // unmark mask bits that are used by other subfunctions to prevent side effects
    FOREACHTRIANGLE(t, n) { UNMARK_BIT(t, 1); UNMARK_BIT(t, 2); UNMARK_BIT(t, 3); }

    // remove the parts of the 2nd mesh that are inside the first one
    FOREACHTRIANGLE(t, n) if(t->mask & (1<<6)) MARK_VISIT(t);
    this->removeSelectedTriangles();
    this->eulerUpdate();
    this->removeSmallestComponents(2);
    if(this->shells() != 2) return 0;

    // clone the parts of the 1st mesh that were inside the 2nd one
    FOREACHTRIANGLE(t, n) if(t->mask & (1<<7)) { MARK_VISIT(t); UNMARK_BIT(t, 7); }
    Triangulation* stin=createSubMeshFromSelection(NULL,true); // to do: replace by createSeparSurfFromSelection
    //JMesh::info("a\n");
    //Triangulation* stin=createSeparSurfFromSelection(); does not yet work; try alternative idea first, delete if successful


    // alternative idea: create 2 complete components with duplicated vertices at separating surface
    // protect vertices of separating surface when cleaning etc by using "protection bit"
    // and save as .stl-file

    // idea: keep pointers from original to duplicated vertex --> use for vertex numbering during saveMSH

//JMesh::info("b\n");
    // copy 1st mesh, then delete it from this
    FOREACHTRIANGLE(t, n) {if (IS_BIT(t,5)) break; }
    ExtTriMesh *shellBckp;
    shellBckp = (ExtTriMesh*) this->extractShell(t);

    FOREACHTRIANGLE(t, n) { UNMARK_VISIT(t); if (IS_BIT(t,5)) MARK_VISIT(t); }
    this->removeSelectedTriangles();
    this->eulerUpdate();
    this->removeSmallestComponents(2);

    // invert orientation of cloned part, mark everything of cloned part with bit 7 (as "protection bit")
    Edge *e;
    FOREACHVEEDGE((&(stin->E)), e, n) UNMARK_VISIT(e);
    FOREACHVTTRIANGLE((&(stin->T)), t, n) {
        MARK_BIT(t,7);
        t->invert();
        if(!IS_VISITED(t->e1)) { MARK_VISIT(t->e1);  p_swap((void **)(&(t->e1->v1)), (void **)(&(t->e1->v2))); };
        if(!IS_VISITED(t->e2)) { MARK_VISIT(t->e2);  p_swap((void **)(&(t->e2->v1)), (void **)(&(t->e2->v2))); };
        if(!IS_VISITED(t->e3)) { MARK_VISIT(t->e3);  p_swap((void **)(&(t->e3->v1)), (void **)(&(t->e3->v2))); };
        MARK_BIT(t->v1(),7); MARK_BIT(t->v2(),7); MARK_BIT(t->v3(),7);
        MARK_BIT(t->e1,7); MARK_BIT(t->e2,7); MARK_BIT(t->e3,7);
    }


    // debug
    JMesh::info("stin T %d   E %d   V %d\n",stin->T.numels(),stin->E.numels(),stin->V.numels());
    JMesh::info("before T %d   E %d   V %d\n",this->T.numels(),this->E.numels(),this->V.numels());

    stin->saveMSH("beforeJoin.msh");
    this->saveMSH("attachtst1.msh");
    //shellBckp->saveMSH("shellbcp.msh");

    // add cloned part of 1st mesh to rest of 2nd mesh, fill boundary
    this->joinTailTriangulation(stin);

    FOREACHTRIANGLE(t, n) { if(IS_BIT(t,4)) break; }
    ComponentStruct c1(t);
    FOREACHTRIANGLE(t, n) { if(IS_BIT(t,7)) break; }
    ComponentStruct c2(t);
    int ret = this->joinComponentsBiggestBoundaryPair(c1.triangles, c2.triangles, 10);
    // done: what are the effects of "fair selection" and "refinement"
    // --> size of joining triangles is more equally distributed; but more intersections

    // debug
    JMesh::info("after T %d   E %d   V %d\n",this->T.numels(),this->E.numels(),this->V.numels());
     int i = 0;
     FOREACHTRIANGLE(t, n) if(t->mask & (1<<7)) i++;
     FOREACHVERTEX(v,n) { v->info = NULL; }
      JMesh::info("nmarked %d\n",i);
     //this->saveMSH("attachtst2.msh",1);


     // cleaning etc (this will still create intersections with original 1st mesh
     // to do: clean and checkRepair have to respect the "protection" bit
     // if cleaning & repairing does not succeed then, maybe cut more of from the second surface??
     //-> so: write cleanP and checkRepairP functions
     // to do: proper handling of  no-clean attribute



     this->unmarkEverything(128); // keep "protection" bit 7

     //this->saveMSH("attachtst2b.msh",1);
     //JMesh::report_progress("Cleaning ...");
     //JMesh::quiet = true;
     this->cleanP();
     this->saveMSH("attachtst3.msh",1);

     //this->checkAndRepair(); // isn't this the last step of clean() anyway?? --> no, not always

     // append 1st mesh again
     this->joinTailTriangulation(shellBckp);

     FOREACHVERTEX(v,n) {
             if ((v->info != NULL)&&IS_BIT( ((Vertex *) v->info) ,7)) {
                MARK_BIT(v,7); UNMARK_BIT(((Vertex *) v->info) ,7);
             }
     }
     this->saveMSH("attachtst4.msh",0,0,0,1);

     // done: update vertex numbering in saveMSH, so that vertices of separation surface do not get a number twice

     // debug


    return 0;
}
*/


int ExtTriMesh::joinComponentsBiggestBoundaryPair(List *nl, List *ml, double joinDistance) {
    ComponentStruct cn(nl), cm(ml);
    cn.initializeBoundaries();
    cm.initializeBoundaries();
    List *loop, *biggestLoop1 = NULL, *biggestLoop2 = NULL;
    // get the pair of the biggest boundary loops
    List *tmp = new List();
    while (loop = (List*) cn.boundaries->popHead()) {
        if(biggestLoop1 && loop->numels() <= biggestLoop1->numels()) tmp->appendHead(loop);
        else biggestLoop1 = loop;
    }
    cn.boundaries->joinTailList(tmp);
    while (loop = (List*) cm.boundaries->popHead()) {
        if(biggestLoop2 && loop->numels() <= biggestLoop2->numels()) tmp->appendHead(loop);
        else biggestLoop2 = loop;
    }
    cm.boundaries->joinTailList(tmp);
    delete tmp;

    bool ret = joinBoundaryPair(biggestLoop1, biggestLoop2);
    this->eulerUpdate();
    this->fillSmallBoundaries(this->E.numels());
    cn.clear();
    cm.clear();
    return ret;
}

bool ExtTriMesh::joinBoundaryPair(List *bl1, List *bl2) {
    Vertex *v, *w; Node *n;
    FOREACHVVVERTEX(bl1, v, n) {
        double d = getClosestPartner(v, bl2, &w);
        //if(joinBoundaryLoops(v, w, false, true, false))
        //    return true;
        if(joinBoundaryLoops2(v, w, false, false)) //org: if(joinBoundaryLoops2(v, w, false, false))
            return true;
    }
    return false;
}

bool ExtTriMesh::loopsHaveAllVerticesCloserThanDistance(List *loop, List *loop2, const double &distance) {
    Node *n, *m;
    Vertex *v, *w;
    const double d2 = distance*distance;
    if (loop->numels() && loop2->numels()) {
        FOREACHVVVERTEX(loop, v, n) {
            bool foundClosePartner = false;
            FOREACHVVVERTEX(loop2, w, m) {
                if( v->squaredDistance(w) < d2 ) {
                    foundClosePartner = true;
                    break;
                }
            }
            if (!foundClosePartner) return false;
        }
        return true;
    }
    return false;
}

double ExtTriMesh::closestPair(List *l1, List *l2, Vertex **closest1, Vertex **closest2)
{
    Node *n, *m;
    Vertex *v,*w;
    double adist, mindist = DBL_MAX;
    FOREACHVVVERTEX(l1, v, n) {
        FOREACHVVVERTEX(l2, w, m) {
            if ((adist = w->squaredDistance(v)) < mindist) {
                mindist=adist;
                *closest1 = v;
                *closest2 = w;
            }
        }
    }
    return mindist;
}

double ExtTriMesh::mostDistantPartner(Vertex *v, List *l, Vertex **distantPartner) {
    Node *m;
    Vertex *w;
    double adist, maxdist = 0;
    FOREACHVVVERTEX(l, w, m) {
        if ((adist = w->squaredDistance(v)) > maxdist) {
            maxdist=adist;
            *distantPartner = w;
        }
    }
    return maxdist;
}

double ExtTriMesh::getClosestPartner(Vertex *v, List *l, Vertex **closestParnter) {
    Node *m;
    Vertex *w;
    double adist, mindist = DBL_MAX;
    FOREACHVVVERTEX(l, w, m) {
        if ((adist = w->squaredDistance(v)) < mindist) {
            mindist=adist;
            *closestParnter = w;
        }
    }
    return mindist;
}

int ExtTriMesh::moveVerticesInwards(Point &componentCenter, std::map<Vertex *, Point> &origin, double stepsize, double distance) {
    List todo(new di_cell(this)), cells;
    di_cell *c, *c2;
    int i = 0;
    int ret = 0;
    double stepsize2 = stepsize*stepsize;
    while (c = (di_cell *)todo.popHead()) {
        if (i > DI_MAX_NUMBER_OF_CELLS || c->triangles.numels() <= 10 || (c->Mp-c->mp).length() < distance) cells.appendHead(c);
        else {
            i++;
            JMesh::report_progress(NULL);
            c2 = c->fork();
            if (!c->containsBothShells(1,2)) delete(c); else todo.appendTail(c);
            if (!c2->containsBothShells(1,2)) delete(c2); else todo.appendTail(c2);
        }
    }
    JMesh::report_progress("");
    Node *n, *m;
    Triangle *t, *t2;
    Vertex *v1, *v2;
    double distance2 = distance*distance;
    while (c = (di_cell *)cells.popHead()) {
        std::set<Vertex*> vertices;
        FOREACHVTTRIANGLE((&c->triangles), t, n) { vertices.insert(t->v1()); vertices.insert(t->v2()); vertices.insert(t->v3()); }
        for(std::set<Vertex*>::const_iterator itv1 = vertices.begin(); itv1 != vertices.end(); ++itv1) {
            v1 = *itv1;
            if(IS_BIT(v1, 2)) {
                FOREACHVTTRIANGLE((&c->triangles), t, n) if(IS_BIT(t, 1)) {
                    Point center = t->getCircleCenter();
                    double radius2 = MAX(center.squaredDistance(t->v1()),distance2);
                    if(center.squaredDistance(v1) < radius2) {
                        Point n = t->getNormal()*distance, p = center + n;
                        v2->intersectionWithPlane(v1, &origin[v1], &p, &n);
                        if(v1->squaredDistance(&origin[v1]) > v2->squaredDistance(&origin[v1]))
                            *v1 += *v2 - *v1;
                        UNMARK_BIT(v1,2);
                    }
                }
            }
        }
    }
    FOREACHVERTEX(v1, n) if(IS_BIT(v1,2)) {
        double dist = origin[v1].distance(v1);
        if(dist > stepsize) *v1 += (origin[v1]-(*v1))*(stepsize/dist);
        else { *v1 = origin[v1]; UNMARK_BIT(v1,2); }
        ret++;
    }
    return ret;
}

bool ExtTriMesh::decoupleFirstFromSecondComponent(double minAllowedDistance, unsigned max_iterations, bool treatFirstAsOuter, bool outwards, double stepWidth) {
    bool quiet = JMesh::quiet;
    int iteration_counter = 0;
    short constantBit = 4, decoupleBit = 5 , markBit = 6;
    Triangle *t;
    Node *n;
    Vertex *v;

    if(stepWidth != 0.5) JMesh::info("decoupling step width set to %f\n",stepWidth);

    ExtTriMesh *shellToDecouple, *constantShell; // temporary triangulation for intermediate cleaning etc.
    if(this->shells() != 2) JMesh::error("Must have exactly 2 components.\n");
    shellToDecouple = (ExtTriMesh*) this->extractFirstShell();
    constantShell = (ExtTriMesh*) this->extractFirstShell();
    this->joinTailTriangulation(shellToDecouple);
    delete(shellToDecouple);
    // dilate the constant component by d, to contrain the minAllowedDistance
    constantShell->dilate((treatFirstAsOuter? 1:-1)*minAllowedDistance);
    JMesh::quiet = true; constantShell->clean(); JMesh::quiet = quiet;
    while(iteration_counter++ < max_iterations) {
        this->selectAllTriangles(decoupleBit); // mark component to decouple
        unsigned toDecoupleTriangleNumber = this->T.numels();
        constantShell->selectAllTriangles(constantBit); // mark constant component
        this->joinHeadTriangulation(constantShell);
        JMesh::info("Iteration %d\n", iteration_counter);

        // A.T.: the markBit should not be 0, as 0 is also used to mark the intersecting triangles
        // in markTrianglesInsideComponent
        unsigned nt = this->markTrianglesInsideComponent(markBit, decoupleBit, constantBit, !treatFirstAsOuter && !outwards);
        // A.T.: the last argument !treatFirstAsOuter && !outwards ensures that intersection triangles are treated as outside when the inner shell
        // moved inside, and treated as inside when the outer shell is moved (both for out-out and out-in case)
        FOREACHTRIANGLE(t, n)  { if (IS_BIT(t,markBit)) MARK_BIT(t,0); else UNMARK_BIT(t,0); }
        // unmark mask bits that are used by other subfunctions to prevent side effects
        FOREACHTRIANGLE(t, n) { UNMARK_BIT(t, 1); UNMARK_BIT(t, 2); UNMARK_BIT(t, 3); }

        constantShell = (ExtTriMesh*) this->extractFirstShell();
        if((treatFirstAsOuter && nt == 0 && outwards) || // outer outwards
           (treatFirstAsOuter && nt == 0 && !outwards) || // outer inwards
           (!treatFirstAsOuter && nt == toDecoupleTriangleNumber && !outwards) // inner inwards
           ) break; // finished
        std::map<Vertex*, Point> shift;
        FOREACHVERTEX(v, n) UNMARK_VISIT(v);

        // we have overlapping triangles
        FOREACHTRIANGLE(t, n) { // IS_BIT(t,0) == triangle is inside the constant component
            if(( IS_BIT(t,0) &&  treatFirstAsOuter &&  outwards) || // outer outwards (move inside  triangles out)
               ( IS_BIT(t,0) &&  treatFirstAsOuter && !outwards) || // outer inwards  (move inside  triangles in)
               (!IS_BIT(t,0) && !treatFirstAsOuter && !outwards)    // inner inwards  (move outside triangles in)
               ) { // compute shift for affected vertices
                MARK_VISIT(t->v1()); MARK_VISIT(t->v2()); MARK_VISIT(t->v3());
            }
            UNMARK_BIT(t,0);
        }

        // compute shift as mean weighted normal of surrounding triangles
        FOREACHVERTEX(v, n) if(IS_VISITED(v)) {
            shift[v] = v->getNormal()*stepWidth*(outwards? 1 : -1);
            UNMARK_VISIT(v);
        }
        for(std::map<Vertex*, Point>::iterator it = shift.begin(); it != shift.end(); ++it) {
            v = it->first;
            *v += it->second;
        }
        this->unmarkEverything();
        JMesh::report_progress("Cleaning ...");
        JMesh::quiet = true;
        this->clean(); // and clean and repair it (because the shift could have produced new intersections ...)
        this->checkAndRepair();
        JMesh::quiet = quiet;
        JMesh::report_progress("");
    }
    if(iteration_counter < max_iterations) return true;
    return false;
}


void ExtTriMesh::cutFirstWithSecondComponent(double minAllowedDistance, bool cutOuter) {
    Triangle *t;
    Node *n;
    bool quiet = JMesh::quiet;
    short constantBit = 4, cutBit = 5, markBit=6;
    if(this->shells() != 2) JMesh::error("Must have exactly 2 components.\n");
    ExtTriMesh *shellToCut, *constantShell; // temporary triangulation for intermediate cleaning etc.
    shellToCut = (ExtTriMesh*) this->extractFirstShell();
    constantShell = (ExtTriMesh*) this->extractFirstShell();
    this->joinTailTriangulation(shellToCut);
    delete(shellToCut);
    // dilate the inner component by d, to contrain the minAllowedDistance
    constantShell->dilate(-1*minAllowedDistance);
    JMesh::quiet = true; constantShell->clean(); JMesh::quiet = quiet;
    this->selectAllTriangles(cutBit);
    constantShell->selectAllTriangles(constantBit);
    this->joinHeadTriangulation(constantShell);

    // A.T.: the markBit should not be 0, as 0 is also used to mark the intersecting triangles
    // in markTrianglesInsideComponent
    unsigned nt = this->markTrianglesInsideComponent(markBit, cutBit, constantBit, cutOuter);
    // A.T.: the last argument cutOuter ensures that intersection triangles are treated as outside when the inner shell is cut
    // and as inside when the outer shell is cut
    FOREACHTRIANGLE(t, n)  { if (IS_BIT(t,markBit)) MARK_BIT(t,0); else UNMARK_BIT(t,0); }
    // unmark mask bits that are used by other subfunctions to prevent side effects
    FOREACHTRIANGLE(t, n) {UNMARK_BIT(t, 1); UNMARK_BIT(t, 2); UNMARK_BIT(t, 3); }

    constantShell = (ExtTriMesh*) this->extractFirstShell(); // extract constant component
    if(cutOuter) this->invertSelection();
    this->removeSelectedTriangles();
    this->unmarkEverything();
    this->checkAndRepair();
    this->clean();
}


void ExtTriMesh::refineFirstWithinSecondComponent(double factor) {
    Triangle *t;
    Node *n;
    ExtTriMesh *shellToRefine, *constantShell; // temporary triangulation for intermediate cleaning etc.

    bool quiet = JMesh::quiet;
    short constantBit = 4, firstCompBit = 5, markBit=6;
    if(this->shells() != 2) JMesh::error("Must have exactly 2 components.\n");

    shellToRefine = (ExtTriMesh*) this->extractFirstShell();
    constantShell = (ExtTriMesh*) this->extractFirstShell();
    this->joinTailTriangulation(shellToRefine);
    delete(shellToRefine);

    this->selectAllTriangles(firstCompBit);
    constantShell->selectAllTriangles(constantBit);
    this->joinHeadTriangulation(constantShell);

    // A.T.: the markBit should not be 0, as 0 is also used to mark the intersecting triangles
    // in markTrianglesInsideComponent
    unsigned nt = this->markTrianglesInsideComponent(markBit, firstCompBit, constantBit, false);
    // A.T.: the last argument cutOuter ensures that intersection triangles are treated as outside when the inner shell is cut
    // and as inside when the outer shell is cut
    FOREACHTRIANGLE(t, n)  { if (IS_BIT(t,markBit)) MARK_BIT(t,0); else UNMARK_BIT(t,0); }
    // unmark mask bits that are used by other subfunctions to prevent side effects
    FOREACHTRIANGLE(t, n) {UNMARK_BIT(t, 1); UNMARK_BIT(t, 2); UNMARK_BIT(t, 3); }

    constantShell = (ExtTriMesh*) this->extractFirstShell(); // extract constant component
    this->refineSelectedTriangles(factor);
    this->unmarkEverything();
    this->checkAndRepair();
    this->clean();
}


int ExtTriMesh::refineSelectedTriangles(double factor) {

     Node *n;
     Triangle *t, *t1, *t2;
     Edge *e, *f;
     Vertex *v;
     List toswap, reg, all_edges, interior_edges, boundary_edges, boundary_vertices, interior_vertices;
     double  l;
     int swaps, totits, ntb, nnt=-1, pnnt, gits=0;
     Point vc;
     double meanTrArea=0.0, minAng1, minAng2, hlpAng;

     FOREACHTRIANGLE(t, n) if (IS_VISITED(t)) reg.appendHead(t);

     // get interior and boundary edges
     FOREACHVTTRIANGLE((&reg), t, n)
     {
      e = t->e1; if (!IS_VISITED2(e)) {MARK_VISIT2(e); all_edges.appendHead(e);} else UNMARK_VISIT2(e);
      e = t->e2; if (!IS_VISITED2(e)) {MARK_VISIT2(e); all_edges.appendHead(e);} else UNMARK_VISIT2(e);
      e = t->e3; if (!IS_VISITED2(e)) {MARK_VISIT2(e); all_edges.appendHead(e);} else UNMARK_VISIT2(e);

      meanTrArea+=t->area();
     }
     meanTrArea=meanTrArea/reg.numels();

     while (all_edges.numels())
     {
      e = (Edge *)all_edges.popHead();
      if (IS_VISITED2(e)) {boundary_edges.appendHead(e); UNMARK_VISIT2(e);}
      else {interior_edges.appendHead(e); MARK_VISIT2(e);}
     }

     // get boundary vertices
     FOREACHVEEDGE((&boundary_edges), e, n)
     {
      v = e->v1; if (!IS_VISITED2(v)) {MARK_VISIT2(v); boundary_vertices.appendHead(v);}
      v = e->v2; if (!IS_VISITED2(v)) {MARK_VISIT2(v); boundary_vertices.appendHead(v);}
     }

     FOREACHVVVERTEX((&boundary_vertices), v, n) UNMARK_VISIT2(v);
     FOREACHVEEDGE((&interior_edges), e, n) UNMARK_VISIT2(e);
     FOREACHVEEDGE((&boundary_edges), e, n) MARK_BIT(e, 3);

     // split triangles as long as area>meanTrArea/factor
     do
     {
      pnnt=nnt;
      nnt=0;
      FOREACHVTTRIANGLE((&reg), t, n)
      {
       vc = t->getCenter();

       if (t->area()>meanTrArea/factor)
       {
        ntb = T.numels();
        v = splitTriangle(t,&vc,1); // only split when angles of new triangles are "proper"
        nnt += (T.numels()-ntb);
        if (T.numels() == ntb+2) // check whether split was really performed
        {
         interior_vertices.appendHead(v);
         interior_edges.appendHead(v->e0);
         interior_edges.appendHead(v->e0->leftTriangle(v)->prevEdge(v->e0));
         interior_edges.appendHead(v->e0->rightTriangle(v)->nextEdge(v->e0));
         t1 = ((Triangle *)T.head()->data);
         t2 = ((Triangle *)T.head()->next()->data);
         t1->mask = t2->mask = t->mask;
         reg.appendHead(t1); reg.appendHead(t2);
        }
       }
      }

      // swap edges to distribute egde lengths more uniformly
      FOREACHVEEDGE((&interior_edges), e, n) {MARK_VISIT2(e); toswap.appendHead(e);}
      totits=0; swaps=1;
      while (swaps && totits++ < 10)
      {
       swaps = 0;
       while ((e=(Edge *)toswap.popHead())!=NULL)
       {
        UNMARK_VISIT2(e);
        l = e->squaredLength();
        if (e->swap())
        {
         if (e->squaredLength() >= l*0.999999) e->swap(1);
         else
         {
          swaps++;
          toswap.appendTail(e);
          f = e->t1->nextEdge(e); if (!IS_VISITED2(f) && !IS_BIT(f, 3)) {MARK_VISIT2(f); toswap.appendTail(f);}
          f = e->t1->prevEdge(e); if (!IS_VISITED2(f) && !IS_BIT(f, 3)) {MARK_VISIT2(f); toswap.appendTail(f);}
          f = e->t2->nextEdge(e); if (!IS_VISITED2(f) && !IS_BIT(f, 3)) {MARK_VISIT2(f); toswap.appendTail(f);}
          f = e->t2->prevEdge(e); if (!IS_VISITED2(f) && !IS_BIT(f, 3)) {MARK_VISIT2(f); toswap.appendTail(f);}
         }
        }
       }
      }

      if (pnnt==nnt) gits++;
     } while (nnt && gits<10);

     FOREACHVEEDGE((&boundary_edges), e, n) UNMARK_BIT(e, 3);

     FOREACHVEEDGE((&boundary_edges), e, n) {
         t1=e->t1; t2=e->t2;
         if(t1 == NULL || t2 == NULL) continue;

         minAng1 = fabs(t1->getAngle(t1->v1()));
         hlpAng = fabs(t1->getAngle(t1->v2())); minAng1 = minAng1<hlpAng ? minAng1 : hlpAng;
         hlpAng = fabs(t1->getAngle(t1->v3())); minAng1 = minAng1<hlpAng ? minAng1 : hlpAng;
         hlpAng = fabs(t2->getAngle(t2->v1())); minAng1 = minAng1<hlpAng ? minAng1 : hlpAng;
         hlpAng = fabs(t2->getAngle(t2->v2())); minAng1 = minAng1<hlpAng ? minAng1 : hlpAng;
         hlpAng = fabs(t2->getAngle(t2->v3())); minAng1 = minAng1<hlpAng ? minAng1 : hlpAng;

         e->swap();

         minAng2 = fabs(t1->getAngle(t1->v1()));
         hlpAng = fabs(t1->getAngle(t1->v2())); minAng2 = minAng2<hlpAng ? minAng2 : hlpAng;
         hlpAng = fabs(t1->getAngle(t1->v3())); minAng2 = minAng2<hlpAng ? minAng2 : hlpAng;
         hlpAng = fabs(t2->getAngle(t2->v1())); minAng2 = minAng2<hlpAng ? minAng2 : hlpAng;
         hlpAng = fabs(t2->getAngle(t2->v2())); minAng2 = minAng2<hlpAng ? minAng2 : hlpAng;
         hlpAng = fabs(t2->getAngle(t2->v3())); minAng2 = minAng2<hlpAng ? minAng2 : hlpAng;

         if(minAng2<minAng1) e->swap();
     }

     if (gits>=10) {JMesh::warning("refineSelectedTriangles: Refinement stage failed to converge. Breaking.\n"); return 1;}

     return 0;
}


int ExtTriMesh::markTrianglesInsideComponent(short insideMarkBit, short componentMarkBit1, short componentMarkBit2, bool treatIntersectionsAsOutside) {

    Triangle *t, *tHlp, *t0, *t1, *t2, *t3;
    Node *n, *nHlp;
    Vertex *vHlp;
    Point BBoxMin, BBoxMax, TCtr;
    int counter;

    if ( insideMarkBit<3 || componentMarkBit1<3 || componentMarkBit2<3) JMesh::warning("markTrianglesInsideComponent: All bits in argument line should be > 3 to prevent side effects!\n");

    FOREACHTRIANGLE(tHlp, nHlp) { UNMARK_BIT(tHlp,2); } // ensure mask bit 2 is zero; otherwise side effects can occur (e.g. when BBoxMarkBit is also 2)
    this->forceNormalConsistence(); // uses mark bit 2, and sets mark bit 2 to 0 as last step
    di_cell *c = new di_cell(this), *c2, *tmp;
    List todo(c), cells, tmptl;

    // keep only triangles of the two components
    while(t = (Triangle*) c->triangles.popHead())
        if(t->mask & (1<<componentMarkBit1 | 1<<componentMarkBit2))
            tmptl.appendHead(t);
    c->triangles.joinTailList(&tmptl);

    int ncells = 0;
    int ret = 0;
    // get smallest cells containing at least both shells
    while (c = (di_cell *)todo.popHead()) {
        if (ncells > 10*DI_MAX_NUMBER_OF_CELLS || c->triangles.numels() <= 100) cells.appendHead(c);
        else {
            JMesh::report_progress(NULL);
            tmp = new di_cell(*c);
            c2 = c->fork();
            if (!c->containsBothShells(componentMarkBit1, componentMarkBit2) ||
                !c2->containsBothShells(componentMarkBit1, componentMarkBit2)) {
                delete(c);
                delete(c2);
                cells.appendHead(tmp);
            } else {
                ncells++;
                todo.appendTail(c);
                todo.appendTail(c2);
                delete(tmp);
            }
        }
    }
    JMesh::report_progress("");
    // if no intersections, then all triangles
    int nintersections = this->selectIntersectingTriangles();

    std::set<Vertex*> vertices1, vertices2;
    short outsideMarkBit = 0;
    // get first unused bit
    unsigned char mask = 1<<componentMarkBit1 | 1<<componentMarkBit2 | 1<<insideMarkBit | 1;
    while(1<<outsideMarkBit & mask) outsideMarkBit++;

    // get second unused bit
    short BBoxMarkBit = 0;
    mask = 1<<componentMarkBit1 | 1<<componentMarkBit2 | 1<<insideMarkBit | 1<<outsideMarkBit | 1;
    while(1<<BBoxMarkBit & mask) BBoxMarkBit++;
    FOREACHTRIANGLE(tHlp, nHlp) { UNMARK_BIT(tHlp,BBoxMarkBit); } // ensure BBox mask bit is zero to prevent side effect (e.g. from forceNormalConsistence above)

    // get third unused bit
    short BBoxVisitBit = 0;
    mask = 1<<componentMarkBit1 | 1<<componentMarkBit2 | 1<<insideMarkBit | 1<<outsideMarkBit | 1<< BBoxMarkBit | 1;
    while(1<<BBoxVisitBit & mask) BBoxVisitBit++;

    unsigned char decidedMask = 1<<insideMarkBit | 1<<outsideMarkBit;

    //JMesh::info("markTrianglesInsideComponent: mask %d; componentMarkBit1 %d; componentMarkBit2 %d; insideMarkBit %d; outsideMarkBit %d; BBoxMarkBit %d; BBoxVisitBit %d; decidedMask %d\n",
    //           mask, componentMarkBit1,componentMarkBit2,insideMarkBit,outsideMarkBit,BBoxMarkBit, BBoxVisitBit, decidedMask);

    while(c = (di_cell*) cells.popHead()) {
        JMesh::report_progress("%d%%", (int)round(((double)(ncells - cells.numels()))/ncells*100));
        // get vertices of triangles of component1 in the cell
        Vertex *vt;
        FOREACHVTTRIANGLE((&c->triangles), t, n) {
            bool comp1 = IS_BIT(t, componentMarkBit1);
            vt = t->v1(); if(!comp1) vertices2.insert(vt); else if(!(vt->mask & decidedMask)) vertices1.insert(vt);
            vt = t->v2(); if(!comp1) vertices2.insert(vt); else if(!(vt->mask & decidedMask)) vertices1.insert(vt);
            vt = t->v3(); if(!comp1) vertices2.insert(vt); else if(!(vt->mask & decidedMask)) vertices1.insert(vt);
        }
        // decide for each vertex whether inside or outside
        for(std::set<Vertex*>::const_iterator i = vertices1.begin(); i != vertices1.end(); ++i) {
            Vertex *v = *i, *w, *closest;
            // search for yet undecided (unmarked) connected regions
            if(!(v->mask & decidedMask)) {
                // get closest vertex of the other component
                double d, dmin = DBL_MAX;
                for(std::set<Vertex*>::const_iterator j = vertices2.begin(); j != vertices2.end(); ++j) {
                    w = *j;
                    if(!IS_VISITED(w)) {
                        MARK_VISIT(w);
                        d = v->squaredDistance(w);
                        if(d < dmin) { dmin = d; closest = w; }
                    }
                }
                for(std::set<Vertex*>::const_iterator j = vertices2.begin(); j != vertices2.end(); ++j)
                    UNMARK_VISIT(*j);
                if(dmin == DBL_MAX) { /*v->printPoint(); */continue; }
                // get the mean normal at the closest vertex
                Point trianglesNormal = ((Triangle *) closest->VT()->head()->data)->getNormal();
                Point meanNormal = closest->getNormal();
                // decide whether it is inside or outside (using the mean normal plane)
                bool isInside = meanNormal.squaredLength() ? meanNormal*(*v - *closest) < 0 : false;
                bool isInside2 = trianglesNormal.squaredLength() ? trianglesNormal*(*v - *closest) < 0 : false;
                // if the normals are too different, the test results differ too
                if(isInside != isInside2) continue; // skip in that case
                // find an unselected triangle of having vertex v
                List *l = v->VT();
                while(t = (Triangle*) l->popHead()) if(!IS_VISITED(t)) { todo.appendHead(t); break; };
                delete(l);
                if(!todo.numels()) continue; // no unselected triangle in neighborhood
                // spread the decision to all connected triangles, stop at selected triangles (== intersections)
                short markbit = isInside ? insideMarkBit : outsideMarkBit;
                while(t = (Triangle*) todo.popHead()) {
                    MARK_BIT(t, markbit);
                    MARK_BIT(t->v1(), markbit); MARK_BIT(t->v2(), markbit); MARK_BIT(t->v3(), markbit);
                    Triangle *t1 = t->t1(), *t2 = t->t2(), *t3 = t->t3();
                    // stop at selected triangles (== markbit 0 == intersections)
                    if(t1 && !(t1->mask & (decidedMask | 1<<0))) todo.appendHead(t1);
                    if(t2 && !(t2->mask & (decidedMask | 1<<0))) todo.appendHead(t2);
                    if(t3 && !(t3->mask & (decidedMask | 1<<0))) todo.appendHead(t3);
                }
            }
        }
        delete(c);
        vertices1.clear();
    }

    // determine the BBox of all triangles with componentMarkBit2
    BBoxMax.x = -DBL_MAX, BBoxMin.x = DBL_MAX;
    BBoxMax.y = -DBL_MAX, BBoxMin.y = DBL_MAX;
    BBoxMax.z = -DBL_MAX, BBoxMin.z = DBL_MAX;
    FOREACHTRIANGLE(tHlp, nHlp) {
        TCtr = tHlp->getCenter();
        if (tHlp->mask & 1<<componentMarkBit2) {
            if (TCtr.x < BBoxMin.x) BBoxMin.x = TCtr.x;
            if (TCtr.y < BBoxMin.y) BBoxMin.y = TCtr.y;
            if (TCtr.z < BBoxMin.z) BBoxMin.z = TCtr.z;

            if (TCtr.x > BBoxMax.x) BBoxMax.x = TCtr.x;
            if (TCtr.y > BBoxMax.y) BBoxMax.y = TCtr.y;
            if (TCtr.z > BBoxMax.z) BBoxMax.z = TCtr.z;
        }
    }
    //JMesh::info("markTrianglesInsideComponent: BBoxMin %f %f %f; BBoxMax %f %f %f\n",
    //                BBoxMin.x,BBoxMin.y,BBoxMin.z,BBoxMax.x,BBoxMax.y,BBoxMax.z);

    // set BBoxMarkBit for all triangles with componentMarkBit1 & being outside the BBox
    FOREACHTRIANGLE(tHlp, nHlp) {
        if (IS_BIT(tHlp, componentMarkBit1)) {
            TCtr = tHlp->getCenter();
            if ((TCtr.x < BBoxMin.x)||(TCtr.y < BBoxMin.y)||(TCtr.z < BBoxMin.z)||
                (TCtr.x > BBoxMax.x)||(TCtr.y > BBoxMax.y)||(TCtr.z > BBoxMax.z)) {
                MARK_BIT(tHlp, BBoxMarkBit);
            }
        }
    }

    // clear BBoxVisitBit
    FOREACHTRIANGLE(tHlp, nHlp) { UNMARK_BIT(tHlp, BBoxVisitBit);};

    List todoBB;
    // mark bit 0 = intersection!
    // find an unvisited triangle (i.e. not on intersection) with BBoxMarkBit set
    nHlp = T.head(); t0 = NULL;
    while(nHlp) {
        tHlp = (Triangle *)nHlp->data;
        if(IS_BIT(tHlp,BBoxMarkBit)&&!IS_BIT(tHlp,BBoxVisitBit)&&!IS_VISITED(tHlp)) { t0 = tHlp; break; };
        nHlp = nHlp->next();
    };

    // spread BBoxVisitBit across all connected triangles
    while(t0) {
        todoBB.appendHead(t0);
        int ns = 0;

        while (todoBB.numels())
        {
        tHlp = (Triangle *) todoBB.popHead();
        if (!IS_BIT(tHlp,BBoxVisitBit))
         {
          t1 = tHlp->t1(); t2 = tHlp->t2(); t3 = tHlp->t3();

          if (t1 != NULL && !IS_BIT(t1,BBoxVisitBit) &&!IS_VISITED(t1) ) todoBB.appendHead(t1);
          if (t2 != NULL && !IS_BIT(t2,BBoxVisitBit) &&!IS_VISITED(t2) ) todoBB.appendHead(t2);
          if (t3 != NULL && !IS_BIT(t3,BBoxVisitBit) &&!IS_VISITED(t3) ) todoBB.appendHead(t3);

          MARK_BIT(tHlp,BBoxVisitBit);
          ns++;
         }
        }

        // test whether triangles exist with BBoxMarkBit set, but not visited yet (BBoxVisitBit not set)
        nHlp = T.head(); t0 = NULL;
        while(nHlp) {
            tHlp = (Triangle *)nHlp->data;
            if(IS_BIT(tHlp,BBoxMarkBit)&&!IS_BIT(tHlp,BBoxVisitBit)&&!IS_VISITED(tHlp)) { t0 = tHlp; break; };
            nHlp = nHlp->next();
        };
    }

    counter = 0;
    FOREACHTRIANGLE(tHlp, nHlp) if (IS_BIT(tHlp, BBoxMarkBit) && !IS_BIT(tHlp, BBoxVisitBit)) counter++;
    if (counter) JMesh::info("BBox algorithm: %d triangles with BBoxMarkBit, but not with BBoxVisitBit set!\n",counter);

    // all triangles with BBoxVisitBit set are outside
    counter = 0;
    FOREACHTRIANGLE(tHlp, nHlp) {
        if (IS_BIT(tHlp, BBoxVisitBit)) {
            if (IS_BIT(tHlp, insideMarkBit)) counter++;
            MARK_BIT(tHlp, outsideMarkBit);
            UNMARK_BIT(tHlp, insideMarkBit);
        }
    }
    if (counter) JMesh::warning("BBox algorithm corrected %d triangles.\n",counter);

    // code intersection triangles as inside triangles (unless treatIntersectionsAsOutside is true)
    JMesh::report_progress("");
    FOREACHTRIANGLE(t, n) {
        if(IS_BIT(t, componentMarkBit1)) {
            if(IS_VISITED(t) || IS_BIT(t, insideMarkBit)) {

                if(treatIntersectionsAsOutside) {
                    if(IS_VISITED(t)&&!IS_BIT(t, insideMarkBit)) MARK_BIT(t, outsideMarkBit);
                    else ret++;
                } else {
                    MARK_BIT(t, insideMarkBit);
                    ret++;
                }
                UNMARK_VISIT(t);

            } else UNMARK_BIT(t, outsideMarkBit);
            t->v1()->mask = 0; t->v2()->mask = 0; t->v3()->mask = 0;
        } else UNMARK_VISIT(t);
    }

    JMesh::info("Number of triangles inside: %d\n", ret);
    return ret;
}

int ExtTriMesh::moveTooCloseVerticesOutwards(double minAllowedDistance, short componentMarkBit1, short componentMarkBit2) {
    di_cell *c = new di_cell(this), *c2;
    Triangle *t, *t2;
    Node *n, *m;
    List todo(c), cells, tmptl;
    // keep only triangles of the two components
    while(t = (Triangle*) c->triangles.popHead())
        if(t->mask & (1<<componentMarkBit1 | 1<<componentMarkBit2))
            tmptl.appendHead(t);
    c->triangles.joinTailList(&tmptl);
    int ncells = 0;
    int ret = 0;
    // cellsize = sqrt(d^2+d^2+d^2) = sqrt(3*d^3)
    double cellsize2 = 4*3*minAllowedDistance*minAllowedDistance;
    // get smallest cells containing at least both shells
    while (c = (di_cell *)todo.popHead()) {
        JMesh::report_progress(NULL);
        if (ncells > DI_MAX_NUMBER_OF_CELLS || c->triangles.numels() <= 10 || (c->Mp-c->mp).squaredLength() < cellsize2 )
            cells.appendHead(c);
        else {
            ncells++;
            JMesh::report_progress(NULL);
            c2 = c->fork();
            if (c->containsBothShells(componentMarkBit1, componentMarkBit2))
                todo.appendTail(c);
            else delete(c);
            if (c2->containsBothShells(componentMarkBit1, componentMarkBit2))
                todo.appendTail(c2);
            else delete(c2);
        }
    }
    double minAllowedDistance2 = minAllowedDistance*minAllowedDistance;
    std::set<Vertex *> vertices;
    std::map<Vertex *, Point> shift;
    std::map<Vertex *, double> minDist2;
    Vertex *v;
    FOREACHVERTEX(v, n) {
        minDist2[v] = minAllowedDistance2;
        shift[v] = Point();
    }
    while(c = (di_cell*) cells.popHead()) {
        vertices.clear();
        JMesh::report_progress(NULL);
        FOREACHVTTRIANGLE((&c->triangles), t, n) if(IS_BIT(t, componentMarkBit1)) {
            vertices.insert(t->v1());
            vertices.insert(t->v2());
            vertices.insert(t->v3());
        }
        while(t = (Triangle*) c->triangles.popHead()) if(IS_BIT(t, componentMarkBit2)) {
            for(std::set<Vertex *>::iterator i = vertices.begin(); i != vertices.end(); ++i) {
                double dist2 = t->pointTriangleSquaredDistance(*i);
                if (dist2 < minDist2[*i]) {
                    minDist2[*i] = dist2;
                    MARK_VISIT(*i);
                }
            }
        }
        delete(c);
    }
    FOREACHVERTEX(v, n) if(IS_VISITED(v)) {
        List *vtl = v->VT();
        FOREACHVTTRIANGLE(vtl, t2, m) shift[v] += t2->getNormal();
        shift[v] /= vtl->numels();
        delete(vtl);
    }
    FOREACHVERTEX(v, n) if(IS_VISITED(v)) {
        ret++;
        *v += shift[v]*MAX((minAllowedDistance-sqrt(minDist2[v])), MAX(0.1*minAllowedDistance, 1));
        UNMARK_VISIT(v);
    }
    JMesh::report_progress("");
    JMesh::info("Number of too close vertices: %d\n", ret);
    return ret;
}

void ExtTriMesh::dilate(double d) {
    if(d == 0.0) return;
    Vertex *v; Node *n, *m; Triangle *t;
    std::map<Vertex *, Point> shift;
    int nsteps = MAX((int) d,1);
    double step = d/(double)nsteps;
    for(int i = 0; i < nsteps; i++) {
        FOREACHVERTEX(v, n)  {
            shift[v] = Point();
            List *vtl = v->VT();
            FOREACHVTTRIANGLE(vtl, t, m) shift[v] += t->getNormal();
            shift[v] /= vtl->numels()/step;
            delete(vtl);
        }
        FOREACHVERTEX(v, n) *v += shift[v];
        this->clean();
    }
}


// ensure a minimal distance between surfaces
// when the same surface is given as to-be-adjusted and constant shell
// then a minimal distance between close-by surface parts is ensured
// A.T.
int ExtTriMesh::fineTune (double dist, int nsteps, bool secondIn) {

    double shiftDist = dist/((double) nsteps); // in [mm]

    ExtTriMesh *shellToAdjust, *constantShell;
    Vertex *v; Node *n, *m; Triangle *t;
    short constantBit = 4, adjustBit = 5, markBit=6;
    int i, j , k, maskHlp, toAdjustTriangleNumber = 0, nTriOverall;
    double direction = 1.0;
    double dirConstSh = 1.0;

    JMesh::info("Minimal distance: %f; substeps: %d\n",dist,nsteps);

    if (secondIn) {
        JMesh::info("Pushing second component inside\n");
        dirConstSh = -1.0;
    }

    bool treatFirstAsInner = true;
    if (treatFirstAsInner) {
        JMesh::info("Pushing first component inside\n");
        direction = -1.0;
    }

    this->deselectTriangles();
    if(this->shells() != 2) JMesh::error("Must have exactly 2 components.\n");
    shellToAdjust = (ExtTriMesh*) this->extractFirstShell();
    constantShell = (ExtTriMesh*) this->extractFirstShell();
    nTriOverall = shellToAdjust->T.numels();

    this->joinTailTriangulation(shellToAdjust);
    this->selectAllTriangles(adjustBit); // mark component to adjust
    FOREACHVERTEX(v, n) MARK_BIT(v, adjustBit);
    delete(shellToAdjust);

    constantShell->selectAllTriangles(constantBit); // mark constant component
    FOREACHVVVERTEX((&(constantShell->V)), v, n) MARK_BIT(v, constantBit);
    this->joinHeadTriangulation(constantShell);
    delete(constantShell);

    char *t_mask = new char [T.numels()];
    j=0; FOREACHTRIANGLE(t, n) t_mask[j++] = 0;
    std::map<Vertex *, Point> shift;
    for(i=0; i<nsteps; i++) {
      FOREACHVERTEX(v, n) if (IS_BIT(v, constantBit)) {
         shift[v] = Point();
         List *vtl = v->VT();
         FOREACHVTTRIANGLE(vtl, t, m) shift[v] += t->getNormal();
         shift[v] /= (dirConstSh/shiftDist)*vtl->numels();
         delete(vtl);
      }
      FOREACHVERTEX(v, n) if (IS_BIT(v, constantBit)) *v += shift[v];

      FOREACHTRIANGLE(t, n) { UNMARK_BIT(t, 0); UNMARK_BIT(t, 1); UNMARK_BIT(t, 2); UNMARK_BIT(t, 3); UNMARK_BIT(t, markBit); }
      unsigned nt = this->markTrianglesInsideComponent(markBit, adjustBit, constantBit, treatFirstAsInner);
      if (treatFirstAsInner) {
        j=0; k=0; FOREACHTRIANGLE(t, n) { if (!IS_BIT(t,markBit) && IS_BIT(t,adjustBit) && t_mask[j] == 0) { k++; }; j++; };
        if (k > 0.9*nTriOverall) {
            JMesh::warning("%d of %d triangles would be marked: probably too many; skipping iteration\n",k,nTriOverall);
            continue;
        }

        j=0; k=0; FOREACHTRIANGLE(t, n) { if (!IS_BIT(t,markBit) && IS_BIT(t,adjustBit) && t_mask[j] == 0) { t_mask[j] = (char) i+1; k++; }; j++; };
        JMesh::info("%d triangles marked\n",k);
      } else {
        j=0; k=0; FOREACHTRIANGLE(t, n) { if (IS_BIT(t,markBit) && IS_BIT(t,adjustBit) && t_mask[j] == 0) { k++; };  j++; };
        if (k > 0.9*nTriOverall) {
            JMesh::warning("%d of %d triangles would be marked: probably too many; skipping iteration\n",k,nTriOverall);
            continue;
        }

        j=0; k=0; FOREACHTRIANGLE(t, n) { if (IS_BIT(t,markBit) && IS_BIT(t,adjustBit) && t_mask[j] == 0) { t_mask[j] = (char) i+1; k++; };  j++; };
        JMesh::info("%d triangles marked\n",k);
      }
    }

    j=0; FOREACHTRIANGLE(t, n) { if (IS_BIT(t,adjustBit) && t_mask[j]>0) { t->mask = t_mask[j]; toAdjustTriangleNumber++; } else { t->mask = 0; }; j++; };
    JMesh::info("%d triangles will be adjusted\n",toAdjustTriangleNumber);

    constantShell = (ExtTriMesh*) this->extractFirstShell();

    FOREACHVERTEX(v, n) {
        shift[v] = Point();
        List *vtl = v->VT();
        maskHlp = nsteps+1;
        FOREACHVTTRIANGLE(vtl, t, m) {
            shift[v] += t->getNormal();
            if (t->mask>0 && t->mask<maskHlp) maskHlp = t->mask;
        }
        shift[v] /= vtl->numels();
        if (maskHlp>nsteps) shift[v] *= 0.0;
        else shift[v] *= direction*((double) (nsteps+1-maskHlp))*shiftDist;
        delete(vtl);
    }
    FOREACHVERTEX(v, n) *v += shift[v];

    delete [] t_mask;
    return toAdjustTriangleNumber;
}


 bool ExtTriMesh::applyTrafo(const char *trafoFile) {
   FILE *fp;
   int i=0;
   float a,b,c,d;
   Matrix4x4 TMat;
   Vertex *v;
   Node *n;

   // read trafo file
   if ((fp = fopen(trafoFile,"r")) == NULL) JMesh::error("Can't open trafo file '%s'.\n", trafoFile);
   while ((fscanf(fp,"%f %f %f %f",&a,&b,&c,&d) == 4) && (i<4)) {
       TMat.matrix[i][0]=a; TMat.matrix[i][1]=b; TMat.matrix[i][2]=c; TMat.matrix[i++][3]=d; }
   fclose(fp);

   printf("matrix:\n");
   for (i=0; i<4; i++) printf("%f %f %f %f\n",TMat.matrix[i][0], TMat.matrix[i][1], TMat.matrix[i][2], TMat.matrix[i][3]);

   // apply trafo
   FOREACHVERTEX(v, n){ TMat.transform(&v->x, &v->y, &v->z); }
   return true;
 }

#ifndef SPRING_MATCH_H
#define SPRING_MATCH_H

class TriMesh;

void smInitMesh(TriMesh *mesh);
void smInitTarget(TriMesh *mesh);
void smRun(int maxIter = -1);

void smDrawInnerSprings();

#endif
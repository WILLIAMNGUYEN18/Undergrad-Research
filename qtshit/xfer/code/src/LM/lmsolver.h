#ifndef LM_SOLVER_H
#define LM_SOLVER_H

class Keeper;

void steihaugSolver(Keeper* B, double* x, const double trustRadius, bool* boundaryHit);
void solve(const int numVar, const double prec, double trustRadius, double* currSol, bool reAllocate, int maxIteration);

#endif
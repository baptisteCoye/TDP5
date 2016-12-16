#include <stdio.h>
#include <stdlib.h>
#include <perf.h>

#include "datatype.h"
#include "util.h"
#include "TDPConfig.h"

#define FILENAME "../tests/data/file_test2.txt"
#define NB_ITER 100
#define d 100
#define distMin 

int main(int argc, char ** argv){
  int N = 8;
  int n = 4;
  particule * total = malloc(sizeof(particule) * N);  
  vecteur * force = malloc(sizeof(vecteur) * N);
  double * distMin = malloc(sizeof(double) * N);
  double dt;

  fill_bloc(&total[0], n, 0, 20, 20, 20);
  fill_bloc(&total[n], n, 120, 20, 20, 20);

  printf("bloc 1 : \n");
  for (int i = 0; i < n; i++){
    print_particule(total[i]);
  }
  printf("bloc 2 : \n");
  for (int i = n; i < N; i++){
    print_particule(total[i]);
  }
 
  // Initialisation des forces a 0 et des distances à -1
  for (int i = 0; i < nbPartPerProc; i++){
    force[i].x = 0.0;
    force[i].y = 0.0;
    distMin[i] = -1.0;
  }

  P2P(force, total, N, distMin);

  accelerate(total,force,nbPartPerProc);

  dt = determine_dt_forall(total, force, N, distMin, 1);

  move_particules(total,force, N, dt);

  printf("bloc 1 : \n");
  for (int i = 0; i < n; i++){
    print_particule(total[i]);
  }
  printf("bloc 2 : \n");
  for (int i = n; i < N; i++){
    print_particule(total[i]);
  }
}

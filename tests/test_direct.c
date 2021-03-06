#include <stdio.h>
#include <stdlib.h>
#include <perf.h>

#include "datatype.h"
#include "util.h"
#include "TDPConfig.h"

#define FILENAME "../tests/data/file_test2.txt"
#define NB_ITER 100
#define d 100

int main(int argc, char ** argv){
  int N = 8;
  int n = 4;
  particule * total = malloc(sizeof(particule) * N);  
  vecteur * force = malloc(sizeof(vecteur) * N);
  double * distMin = malloc(sizeof(double) * N);
  double dt;

  for (int i = 0; i < N; i++){
    total[i].m = 10;
    total[i].px = 5 + 10*(i%2) + 120*(i/4);
    total[i].py = (5 + 10*(i/2))%20 + 120*(i/4);
  }

  printf("bloc 1 : \n");
  for (int i = 0; i < n; i++){
    print_particule(total[i]);
  }
  printf("bloc 2 : \n");
  for (int i = n; i < N; i++){
    print_particule(total[i]);
  }
 
  // Initialisation des forces a 0 et des distances à -1
  for (int i = 0; i < N; i++){
    force[i].x = 0.0;
    force[i].y = 0.0;
    distMin[i] = -1.0;
  }

  P2P(force, total, N, distMin);

  accelerate(total,force,N);

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

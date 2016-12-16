#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "datatype.h"
#include "util.h"

#define MIN_DISTANCE 0.01

int main(int argc, char ** argv){
  
  if (argc != 3){
    fprintf(stderr, "./randFile <name of the file> <number of particles>\n");
    return 0;
  }

  srand(time(NULL));

  FILE * file = fopen(argv[1], "w+");
  int N = atoi(argv[2]);

  fprintf(file, "%d\n", N);
  particule * parts = malloc(sizeof(parts[0]) * N);

  int i, j;

  for (i = 0; i < N; i++){
    parts[i].m = (double) (rand() % (1000*N)) / N;
    parts[i].px = (double) (rand() % (1000*N)) / N -500;
    parts[i].py = (double) (rand() % (1000*N)) / N -500;
    parts[i].vx = (double) (rand() % (1000*N)) / N -500;
    parts[i].vy = (double) (rand() % (1000*N)) / N -500;
  }

  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      if ((i!=j) && (distance(parts[i], parts[j]) < MIN_DISTANCE)){
	parts[i].px = (double) (rand() % (1000*N)) / N;
	parts[i].py = (double) (rand() % (1000*N)) / N;
	i--;
	break;
      }
    }
  }

  for (i = 0; i < N; i++){
    fprintf(file, "%lf %lf %lf %lf %lf\n", parts[i].m, parts[i].px, parts[i].py, parts[i].vx, parts[i].vy);
  }

  fclose(file);
  free(parts);

  return 0;
}

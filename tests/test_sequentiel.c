#include <stdio.h>
#include <stdlib.h>
#include <perf.h>

#include "datatype.h"
#include "util.h"
#include "TDPConfig.h"

#define FILENAME "../tests/data/file_test2.txt"
#define NB_ITER 100

int main(int argc, void** argv){
  perf_t begin, end;
  perf(&begin);

  int N;
  particule * data;
  vecteur * forces;
  double * distMin;
  double dt;
  int err; 

  FILE * readfile = fopen(argv[1], "r");
  if (readfile == NULL){
    fprintf(stderr, "Erreur lors de la lecture du readfile. Fermeture du programme de test.\n");
    return EXIT_FAILURE;
  }

  err = fscanf(readfile, "%d\n", &N);
  if (err < 0)
    fprintf(stderr, "Erreur lors de la lecture dans le fichier. err = :%d:\n", err);

#if VERBOSE >= 1
  printf("N = %d\n", N);
#endif

  data = malloc(sizeof(particule) * N);
  forces = malloc(sizeof(vecteur) * N);
  distMin = malloc(sizeof(double) * N);

  for (int i = 0; i < N; i++){
    err = fscanf(readfile, "%lf %lf %lf %lf %lf\n", &(data[i].m), &(data[i].px), &(data[i].py), &(data[i].vx), &(data[i].vy));
    if (err < 0)
      fprintf(stderr, "Erreur lors de la lecture dans le fichier. err = :%d:\n", err);
  }

  fclose(readfile);

#ifdef SAVE_RESULTS
  char * initfilename = malloc(sizeof(char) * 20);
  snprintf(initfilename, 20, "save_seq_0.csv");    

  FILE * writefile = fopen(initfilename, "w+");
  if (writefile == NULL) {
    fprintf(stderr, "Erreur lors de l'ouverture du fichier %s en ecriture. Fermeture du programme.\n", initfilename);
    return EXIT_FAILURE;
  }
    
  //    fprintf(writefile, "%d\n", N);
  fprintf(writefile, "X,Y\n");
  for (int i = 0; i < N; i++)
    fprintf(writefile, "%lf,%lf\n", data[i].px, data[i].py);
  //      fprintf(writefile, "%lf %lf %lf %lf %lf\n", data[i].m, data[i].px, data[i].py, data[i].vx, data[i].vy);
      
  fclose(writefile);
  free(initfilename);
#endif

  for (int i = 0; i < NB_ITER; i++){
#if VERBOSE >= 1
    printf("iteration %d\n", i);
#endif

    dt = DT_MAX;

#if VERBOSE >= 1
    printf("     calcul des forces ...\n");
#endif

    for (int i = 0; i < N; i++){
      forces[i].x = 0;
      forces[i].y = 0;
      distMin[i] = -1;
    }

    calcul_local(forces, data, N, distMin);

    double dtTmp;

    accelerate(data, forces, N);

#if VERBOSE >= 1
    printf("     determination du dt ... \n");
#endif

#if VERBOSE >= 2
    printf("           distMin :");
    for (int i = 0; i < N; i++)
      printf("  %lf", distMin[i]);
    printf("\n");
#endif

#if VERBOSE >= 2
    printf("           forces :");
    for (int i = 0; i < N; i++)
      printf("  [%lf, %lf]", forces[i].x, forces[i].y);

    printf("\n");
#endif

#if VERBOSE >= 2
    printf("           positions :");
    for (int i = 0; i < N; i++)
      printf("  [%lf, %lf]", data[i].px, data[i].py);

    printf("\n");
#endif
    
#if VERBOSE >= 2
    printf("           vitesse :");
    for (int i = 0; i < N; i++)
      printf("  [%lf, %lf]", data[i].vx, data[i].vy);

    printf("\n");
#endif

#if VERBOSE >= 2
    printf("           accelerations :");
    for (int i = 0; i < N; i++)
      printf("  [%lf, %lf]", data[i].ax, data[i].ay);

    printf("\n");
#endif

    dt = determine_dt_forall(data, forces, N, distMin, 1);

#if VERBOSE >= 1
    printf("     dt = %lf\n", dt);
#endif

#if VERBOSE >= 1
    printf("     deplacement des particules\n");
#endif

    move_particules(data, forces, N, dt);

#ifdef SAVE_RESULTS

#if VERBOSE >= 1
    printf("     sauvegarde des donnees\n");
#endif

    char * filename = malloc(sizeof(char) * 20);
    snprintf(filename, 20, "save_seq_%d.csv", i+1);

    FILE * writefile = fopen(filename, "w+");
    if (writefile == NULL) {
      fprintf(stderr, "Erreur lors de l'ouverture du fichier %s en ecriture. Fermeture du programme.\n", filename);
      return EXIT_FAILURE;
    }
    
    //    fprintf(writefile, "%d\n", N);
    fprintf(writefile, "X,Y\n");
    for (int i = 0; i < N; i++)
      fprintf(writefile, "%lf,%lf\n", data[i].px, data[i].py);
      //      fprintf(writefile, "%lf %lf %lf %lf %lf\n", data[i].m, data[i].px, data[i].py, data[i].vx, data[i].vy);
      
    fclose(writefile);
    free(filename);
#endif
  }

  free(data); free(forces); free(distMin);
  perf(&end);
  perf_diff(&begin, &end);
  printf("Temps d'execution: "); perf_printmicro(&end);
  return 0;
}

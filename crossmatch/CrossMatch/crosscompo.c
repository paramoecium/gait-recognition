/*
 *  Program name --- crosscompo.c
 */

#include        <stdio.h>
#include        <string.h>
#include        <stdlib.h>
#include	<time.h>
#include	<sys/time.h>
#include	"crossmatch.h"

#define	BUFFER 100
#define FUNCTION 1	// 0: L1 distance, 1: euclid distance

extern double			epsilon;
extern struct SeqList		*list;
extern int			listnum;


/*
 * Allocate memory for data sequence
 */

void malloc_double (struct DataElement *array)
{

  int	len;

  len = array->length;

  array->value = (double *)malloc(sizeof(double) * len);

  if (array->value == NULL) {
   printf ("array->value malloc error\n");
   exit(1);
  }

}


/*
 * Read data sequence
 */

void read_file(
  char *filename,
  double *seq,
  int len
) {

  int j;
  char  read_buf[BUFFER];
  FILE  *fp;

  if ( ( fp = fopen ( filename , "r" ) ) == NULL ){
    printf( "  %s can not open\n" , filename);
    exit(1);
  }
  j = 0;
  while (fgets(read_buf,BUFFER,fp) != NULL ){
    seq[j] = atof(read_buf);
    j++;
    if (j >= len) break;
  }
  if (j < len)
    printf("Warning: sequence is short. length=%d\n", j);
  fclose(fp);

}


/*
 * Compute distance between elements
 */
double element_dist (double *x, double *y) {

  double	dist;

#if FUNCTION
  /* square of euclid distance */
  dist = (*x - *y) * (*x - *y);
#else
  double	abs;
  /* L1 distance */
  abs = *x - *y;
  if (abs < 0) {
    abs = -1 * abs;
  }
  dist = abs;
#endif

  return (dist);

}


/*
 * Compute elapsed time (micro seconds)
 */
long long int get_time_usec (struct timeval start, struct timeval end) {

  long long	sec, usec;

  if (start.tv_usec <= end.tv_usec) {
    sec = end.tv_sec - start.tv_sec;
    usec = end.tv_usec - start.tv_usec;
  } else {
    sec = end.tv_sec - start.tv_sec - 1;
    usec = end.tv_usec - start.tv_usec + 1000000;
  }

  return sec * 1000000 + usec;

}


/*
 * Compute maximum value
 */
double max (double *a, double *b, double *c)
{

  double	maximum;

  if (*a >= *b) {
    maximum = *a;
  } else {
    maximum = *b;
  }

  if (*c > maximum) {
    maximum = *c;
  }

  return (maximum);

}


/*
 * Compute minimum value
 */
double min (double *a, double *b, double *c)
{

  double minmum;

  if (*a <= *b) {
    minmum = *a;
  } else {
    minmum = *b;
  }

  if (*c < minmum) {
    minmum = *c;
  }

  return (minmum);

}

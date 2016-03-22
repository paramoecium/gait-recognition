/*
 *  File name :   naive-scband.c
 */

#include	<stdio.h>
#include	<string.h>
#include	<stdlib.h>
#include	<time.h>
#include	<math.h>
#include	<sys/time.h>
#include	"crossmatch.h"

#define	BUFSIZE	256
#define INITIALIZE_VALUE 1e20
#define	LIST 100		// Upper limit value for candidate subsequence pairs in list arrays
#define EXPERIMENT 1    	// 0: Performance evaluation, 1: Pattern discovery


void	malloc_read();
void	simpair_detection();
void	error(char *);

struct DataElement	**sqd;				/* Data sequences */
struct MatrixElement	*pre_column, *cur_column;	/* Score and position matrices */
struct SeqList		*list;				/* Array for candidate subsequences */
int			Lmin;				/* Minimum length of subsquence match */
double			epsilon;			/* Distance threshold */
int			scwidth;			/* Width of Sakoe-Chiba band */
int			listnum;			/* The number of candidate subsequence pairs in list arrays*/ 


int main(
  int	argc,
  char	*argv[]
) {


  char		dbfile[BUFSIZE], dbfile2[BUFSIZE];
  struct timeval	start, end;
  double	totaltime = 0.0;


  if (argc != 6) {
    error (argv[0]);
    exit(1);
  }
  strcpy(dbfile, argv[1]);	// File name of data sequence X
  strcpy(dbfile2, argv[2]);	// File name of data sequence Y
  Lmin = atoi(argv[3]);		// Minimum length of subsequence match
  epsilon = atof(argv[4]);	// Distance threshold
  scwidth = atof(argv[5]);	// Width of Sakoe-Chiba band
  

  /*
   * Allocate memory and read data sequences 
   */
  malloc_read(dbfile, dbfile2);


  /*
   * Detection of similar subsequence pairs
   */
  gettimeofday (&start, NULL);
  simpair_detection(sqd[0]->value, sqd[1]->value);
  gettimeofday (&end, NULL);


  /*
   * Output computation time
   */
  totaltime = (double)get_time_usec(start, end) / (double)1000000;
#if EXPERIMENT == 0
  printf ("%.8f sec. \n", totaltime);
#endif


  /*
   * Release memory space
   */
  free(sqd[0]);
  free(sqd[1]);
  free(sqd);
  free(pre_column);
  free(cur_column);
  free(list);

  return 0;

}



/*****************************************************************************************
 * Function name	: malloc_read
 * Input1		: char *dbfile	- Data sequence X
 * Input2		: char *dbfile2 - Data sequence Y
 * note			: Allocate memory for the computation and store data sequences
*****************************************************************************************/

void malloc_read (
  char	*dbfile,
  char	*dbfile2
) {


  long int	j;
  char		read_buf[BUFSIZE];
  FILE		*fp;


  sqd = (struct DataElement **)malloc(sizeof(struct DataElement *));
  sqd[0] = (struct DataElement *)malloc(sizeof(struct DataElement));
  sqd[1] = (struct DataElement *)malloc(sizeof(struct DataElement));

  if (sqd == NULL || sqd[0] == NULL || sqd[1] == NULL) {
   printf ("sqd or sqd[0] or sqd[1] malloc error\n");
   exit(1);
  }


  /*
   * Check sequence lengths
   */
  if ((fp = fopen (dbfile, "r")) == NULL) {
   printf ("%s can not open\n", dbfile);
   exit(1);
  }
 
  sqd[0]->length = 0;
  while (fgets(read_buf, BUFSIZE, fp) != NULL) {
   sqd[0]->length++;
  }
  fclose(fp);

  if ((fp = fopen (dbfile2, "r")) == NULL) {
   printf ("%s can not open\n", dbfile2);
   exit(1);
  }

  sqd[1]->length = 0;
  while (fgets(read_buf, BUFSIZE, fp) != NULL) {
   sqd[1]->length++;
  }
  fclose(fp);

  /*
   * Allocate memory for data sequences and read each element
   */
  malloc_double (sqd[0]);
  malloc_double (sqd[1]);
  read_file (dbfile, sqd[0]->value, sqd[0]->length);
  read_file (dbfile2, sqd[1]->value, sqd[1]->length);


  /*
   * Allocate memory for score and position matrices
   */
  pre_column = (struct MatrixElement *)malloc((sqd[1]->length+1) * sizeof(struct MatrixElement));
  cur_column = (struct MatrixElement *)malloc((sqd[1]->length+1)* sizeof(struct MatrixElement));

  if (pre_column == NULL || cur_column == NULL) {
    printf ("matrix malloc error\n");
    exit(1);
  }


  /*
   * Inisialize matrices
   */
  for (j = 0; j < sqd[1]->length+1; j++) {
    pre_column[j].dist = 0;
    pre_column[j].score = 0;
    pre_column[j].x_sp = 0;
    pre_column[j].y_sp = 0;
    cur_column[j].dist = 0;
    cur_column[j].score = 0;
    cur_column[j].x_sp = 0;
    cur_column[j].y_sp = 0;
  }


  /*
   * Allocate and inisialize memory of list arrays for candidate subsequences
   */
  list = (struct SeqList *)malloc(sizeof(struct SeqList));
  if (list == NULL) {
    printf ("list malloc error\n");
    exit(1);
  }

  list[0].flag = 1;
  list[0].x_sp = 0;
  list[0].y_sp = 0;
  list[0].x_ep = 0;
  list[0].y_ep = 0;
  list[0].length = 0;
  list[0].dist = 0;
  list[0].score = 0;


}


/**************************************************************************
 * Function name	: simpair_detection
 * Input1		: double *stream	- Data sequence X
 * Input2		: double *stream2	- Data sequence Y
 * note			: Detect similar subsequence pairs 
**************************************************************************/

void simpair_detection (double *stream, double *stream2)
{

  long int	xpos, ypos;
  long int	i, j;
  long int	band_y_sp, band_y_ep;
  long int	band_yy_sp, band_yy_ep;
  long int	loop_sp_y;
  int		Lx, Ly;
  double	L;
  double	dist_cur;
  double	d_best;
  double	d_min;
  double	d_diag,	d_vert,	d_hori;
  double	diagpt, diagpt2;
  struct MatrixElement	*swap;

  /*
   * Initialize each parameter
   */
  listnum = 0;
  Lx = Ly = 0;
  L = 0;
  dist_cur = INITIALIZE_VALUE;
  d_best = INITIALIZE_VALUE;
  d_min = INITIALIZE_VALUE;
  d_diag = INITIALIZE_VALUE;
  d_vert = INITIALIZE_VALUE;
  d_hori = INITIALIZE_VALUE;


  /*
   * Compute DTW distances
   */

  for (xpos = 1; xpos <= sqd[0]->length; xpos++) { 

    // Update Sacoe-Chiba band
    diagpt = xpos * (double)sqd[1]->length/sqd[0]->length;
    band_y_sp = (int)floor(diagpt - scwidth);
    band_y_ep = (int)ceil(diagpt + scwidth);
    if (band_y_sp < 1)  band_y_sp = 1;
    if (band_y_ep > sqd[1]->length)  band_y_ep = sqd[1]->length;

    for (ypos = 1; ypos <= sqd[1]->length; ypos++) {

      if (ypos < band_y_sp) {
	goto JUMP;
      } else if (band_y_ep < ypos) {
	break;
      }

      dist_cur = element_dist (&stream[xpos-1], &stream2[ypos-1]);
      pre_column[ypos].dist = dist_cur;

      pre_column[ypos].x_sp = xpos;
      pre_column[ypos].y_sp = ypos;

      pre_column[ypos].score = epsilon*1 - pre_column[ypos].dist;
      d_min = pre_column[ypos].dist - epsilon*(1-Lmin);


      for (j = ypos+1; j <= band_y_ep; j++) {
        dist_cur = element_dist (&stream[xpos-1], &stream2[j-1]);
        d_best = pre_column[j-1].dist;

        pre_column[j].dist = dist_cur + d_best;
        pre_column[j].x_sp = xpos;
        pre_column[j].y_sp = pre_column[j-1].y_sp;

	Lx = 1;
	Ly = j - pre_column[j].y_sp + 1;
	L = (double)(Lx + Ly) / 2;
	pre_column[j].score = epsilon*L - pre_column[j].dist;
	d_min = pre_column[j].dist - epsilon*(L-Lmin);
      }

      for (i = xpos+1; i <= sqd[0]->length; i++) {

	diagpt2 = i * (double)sqd[1]->length/sqd[0]->length;
	band_yy_sp = (int)floor(diagpt2 - scwidth);
	band_yy_ep = (int)ceil(diagpt2 + scwidth);
	if (band_yy_sp < 1)  band_yy_sp = 1;
	if (band_yy_ep > sqd[1]->length)  band_yy_ep = sqd[1]->length;
	if (ypos < band_yy_sp) {
	  loop_sp_y = band_yy_sp;
	} else {
	  loop_sp_y = ypos;
	}

	dist_cur = element_dist (&stream[i-1], &stream2[loop_sp_y-1]);
	d_diag = pre_column[loop_sp_y-1].dist;
	d_hori = pre_column[loop_sp_y].dist;

	if (d_diag <= d_hori) {
	  cur_column[loop_sp_y].dist = dist_cur + d_diag;
	  cur_column[loop_sp_y].x_sp = pre_column[loop_sp_y-1].x_sp;
	  cur_column[loop_sp_y].y_sp = pre_column[loop_sp_y-1].y_sp;
	} else {
	  cur_column[loop_sp_y].dist = dist_cur + d_hori;
	  cur_column[loop_sp_y].x_sp = pre_column[loop_sp_y].x_sp;
	  cur_column[loop_sp_y].y_sp = pre_column[loop_sp_y].y_sp;
	}

	Lx = i - cur_column[ypos].x_sp + 1;
	Ly = loop_sp_y - cur_column[loop_sp_y].y_sp + 1;
	L = (double)(Lx + Ly) / 2;

	cur_column[loop_sp_y].score = epsilon*L - cur_column[loop_sp_y].dist;
	d_min = cur_column[loop_sp_y].dist - epsilon*(L-Lmin);


        for (j = loop_sp_y+1; j <= band_yy_ep; j++) {

          dist_cur = element_dist (&stream[i-1], &stream2[j-1]);
	  d_diag = pre_column[j-1].dist; 
	  d_vert = cur_column[j-1].dist; 
	  d_hori = pre_column[j].dist; 

          d_best = min (&d_diag, &d_vert, &d_hori);

          if (d_best == d_diag) {
	    cur_column[j].dist = dist_cur + d_diag;
    	    cur_column[j].x_sp = pre_column[j-1].x_sp;
	    cur_column[j].y_sp = pre_column[j-1].y_sp;
          } else if (d_best == d_vert) {
	    cur_column[j].dist = dist_cur + d_vert;
    	    cur_column[j].x_sp = cur_column[j-1].x_sp;
	    cur_column[j].y_sp = cur_column[j-1].y_sp;
          } else {
	    cur_column[j].dist = dist_cur + d_hori;
    	    cur_column[j].x_sp = pre_column[j].x_sp;
	    cur_column[j].y_sp = pre_column[j].y_sp;
          }

	  Lx = i - cur_column[j].x_sp + 1;
	  Ly = j - cur_column[j].y_sp + 1;
	  L = (double)(Lx + Ly) / 2;

	  cur_column[j].score = epsilon*L - cur_column[j].dist;
	  d_min = cur_column[j].dist - epsilon*(L-Lmin);


#if EXPERIMENT
          if (cur_column[j].dist <= epsilon*(L-Lmin)) {
	    if ((list[0].dist-epsilon*(list[0].length-Lmin)) 
		>= (cur_column[j].dist-epsilon*(L-Lmin))) {
	      list[0].flag = 1;
	      list[0].x_sp = cur_column[j].x_sp;
	      list[0].y_sp = cur_column[j].y_sp;
	      list[0].x_ep = i;
	      list[0].y_ep = j;
	      list[0].dist = cur_column[j].dist;
	      list[0].score = cur_column[j].score;
	    }
	  }
#endif

	  /* Initialization for next computation */
	  L = 0;
	  Lx = Ly = 0;
          dist_cur = INITIALIZE_VALUE;
          d_best = INITIALIZE_VALUE;
	  d_diag = INITIALIZE_VALUE;
	  d_vert = INITIALIZE_VALUE;
	  d_hori = INITIALIZE_VALUE;

        }


	/*
	 * Initialization score matrix
	 */
	diagpt2 = (i+1) * (double)sqd[1]->length/sqd[0]->length;
	band_yy_sp = (int)floor(diagpt2 - scwidth);
	if (band_yy_sp < 1)  band_yy_sp = 1;
	for (j = band_yy_sp-1; j >= band_yy_sp-5; j--) {
 	  if (j == 0) {
	    break;
	  }
	  pre_column[j].dist = INITIALIZE_VALUE;
	}

	/*
	 * Update matrix pointers
	 */
        swap = cur_column;
        cur_column = pre_column;
        pre_column = swap;

      }


#if EXPERIMENT
      /*
       * Report optimal subsequence pair
       */
      if (list[0].flag == 1) {
	printf ("X[%5d:%5d] \t Y[%5d:%5d] \t DTWdist: %.1f\n",
	        list[0].x_sp, list[0].x_ep, list[0].y_sp, list[0].y_ep, list[0].dist);
      }
#endif


      for (j = 0; j < sqd[1]->length+1; j++) {
        pre_column[j].dist = INITIALIZE_VALUE;
        pre_column[j].score = 0;
        pre_column[j].x_sp = 0;
        pre_column[j].y_sp = 0;
        cur_column[j].dist = INITIALIZE_VALUE;
        cur_column[j].score = 0;
        cur_column[j].x_sp = 0;
        cur_column[j].y_sp = 0;
      }

      JUMP:

      list[0].flag = 0;

#if EXPERIMENT
      list[0].flag = 0;
      list[0].x_sp = list[0].y_sp = 0;
      list[0].x_ep = list[0].y_ep = 0;
      list[0].dist = INITIALIZE_VALUE;
      list[0].score = 0;
#endif


    }

  }

}



void error (char *toolname) {

  printf ("Usage: %s  'file of Sequence X'  'file of Sequence Y'  'Lmin'  'epsilon'  'scband width'\n", toolname);

}



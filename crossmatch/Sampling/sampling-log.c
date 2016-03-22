/*
 *  File name :   sampling-log.c
 */

#include	<stdio.h>
#include	<string.h>
#include	<stdlib.h>
#include	<time.h>
#include	<math.h>
#include	<sys/time.h>
#include	"sampling.h"

#define	BUFSIZE	256
#define	LIST 100	// Upper limit value for candidate subsequence pairs in list arrays


// Function declaration
void	malloc_read();
void	simpair_detection();
void	error(char *);

// Grobal variables
struct DataElement	**sqd;				/* Data sequences */
struct MatrixElement	*pre_column, *cur_column;	/* Score and position matrices */
struct SeqList		*list;				/* Array for candidate subsequences */
int			Lmin;				/* Minimum length of subsquence match */
double			epsilon;			/* Distance threshold */
int			scwidth;			/* Width of Sakoe-Chiba band */
int			Tx, Ty;				/* Sampling periods */
int			listnum;			/* The number of candidate subsequence pairs in list arrays*/ 
int			x_start, y_start;		/* Starting position of optimal subsequence pair */  
int			x_end, y_end;			/* End position of optimal subsequence pair */


int main(
  int	argc,
  char	*argv[]
) {


  char			dbfile[BUFSIZE], dbfile2[BUFSIZE];


  if (argc != 12) {
    error (argv[0]);
    exit(1);
  }
  strcpy(dbfile, argv[1]);	// File name of data sequence X
  strcpy(dbfile2, argv[2]);	// File name of data sequence Y
  Lmin = atoi(argv[3]);		// Minimum length of subsequence match
  epsilon = atof(argv[4]);	// Distance threshold
  scwidth = atof(argv[5]);	// Width of Sakoe-Chiba band
  Tx = atoi(argv[6]);		// Sampling period of X
  Ty = atoi(argv[7]);		// Sampling period of Y
  x_start = atoi(argv[8]);	// Starting position of X
  y_start = atoi(argv[9]);	// Starting position of Y
  x_end = atoi(argv[10]);	// End position of X
  y_end = atoi(argv[11]);	// End position of Y


  /*
   * Allocate memory and read data sequences 
   */
  malloc_read(dbfile, dbfile2);


  /*
   * Detection of similar subsequence pairs
   */
  simpair_detection(sqd[0]->value, sqd[1]->value);


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


  long int	i, j;
  char		read_buf[BUFSIZE];
  FILE		*fp;


  /*
   * Allocate memory for data sequences
   */
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
    pre_column[j].x_sp_org = 0;
    pre_column[j].y_sp_org = 0;
    cur_column[j].dist = 0;
    cur_column[j].score = 0;
    cur_column[j].x_sp = 0;
    cur_column[j].y_sp = 0;
    cur_column[j].x_sp_org = 0;
    cur_column[j].y_sp_org = 0;
  }


  /*
   * Allocate and inisialize memory of list arrays for candidate subsequences
   */
  list = (struct SeqList *)malloc(LIST * sizeof(struct SeqList));
  if (list == NULL) {
    printf ("list malloc error\n");
    exit(1);
  }

  for (i = 0; i < LIST; i++) {
    list[i].flag = 1;
    list[i].x_sp = 0;
    list[i].y_sp = 0;
    list[i].x_ep = 0;
    list[i].y_ep = 0;
    list[i].x_sp_org = 0;
    list[i].y_sp_org = 0;
    list[i].x_ep_org = 0;
    list[i].y_ep_org = 0;
    list[i].length = 0;
    list[i].dist = 0;
    list[i].score = 0;
  }


}


/**************************************************************************
 * Function name	: simpair_detection
 * Input1		: double *stream	- Data sequence X
 * Input2		: double *stream2	- Data sequence Y
 * note			: Detect similar subsequence pairs 
**************************************************************************/

void simpair_detection (double *stream, double *stream2)
{

  long int	i, j;
  long int	band_sp, band_ep;
  int		lx, ly;	
  int		max_p;
  double	L;
  double	dist;
  double	s_best, s_diag, s_vert, s_hori;
  double	Bv, Bh, Bd;
  double	diagpt;
  struct MatrixElement	*swap;


  /*
   * Initialize each parameter
   */
  listnum = 0;
  lx = ly = 0;
  L = 0;
  dist = 0;
  s_best = 0;
  s_diag = s_vert = s_hori = 0;
  Bv = (double)Ty / 2;
  Bh = (double)Tx / 2;
  Bd = (double)(Tx + Ty) / 2;

  if (Tx >= Ty) max_p = Tx;
  else max_p = Ty;


  /*
   * Compute scores and starting positions
   */

  for (i = 1; i <= sqd[0]->length; i++) {

    // Update Sacoe-Chiba band
    diagpt = i * (double)sqd[1]->length/sqd[0]->length;
    band_sp = (int)floor(diagpt - scwidth);
    band_ep = (int)ceil(diagpt + scwidth);
    if (band_sp < 1) {
      band_sp = 1;
    }
    if (band_ep > sqd[1]->length) {
      band_ep = sqd[1]->length;
    }

    for (j = band_sp; j <= band_ep; j++) {

      dist = element_dist (&stream[i-1], &stream2[j-1]);
      s_diag = Bd*epsilon - max_p*dist + pre_column[j-1].score;
      s_vert = Bv*epsilon - Ty*dist + cur_column[j-1].score;
      s_hori = Bh*epsilon - Tx*dist + pre_column[j].score;
      s_best = max (&s_diag, &s_vert, &s_hori);

      cur_column[j].score = s_best;

      // Reset score
      if (cur_column[j].score <= 0) {
        cur_column[j].score = 0;
    	cur_column[j].x_sp = i;
	cur_column[j].y_sp = j;
    	cur_column[j].x_sp_org = i * Tx;
	cur_column[j].y_sp_org = j * Ty;
      // Keep score of neighboring cells
      } else {
        if (s_best == s_diag && pre_column[j-1].score != 0) {  // Select diagonal cell
    	  cur_column[j].x_sp = pre_column[j-1].x_sp;
	  cur_column[j].y_sp = pre_column[j-1].y_sp;
    	  cur_column[j].x_sp_org = pre_column[j-1].x_sp_org;
	  cur_column[j].y_sp_org = pre_column[j-1].y_sp_org;
        } else if (s_best == s_vert && cur_column[j-1].score != 0) {  // Select vertical cell
	  cur_column[j].x_sp = cur_column[j-1].x_sp;
	  cur_column[j].y_sp = cur_column[j-1].y_sp;
	  cur_column[j].x_sp_org = cur_column[j-1].x_sp_org;
	  cur_column[j].y_sp_org = cur_column[j-1].y_sp_org;
        } else if (s_best == s_hori && pre_column[j].score != 0) {  // Select horizontal cell
	  cur_column[j].x_sp = pre_column[j].x_sp;
	  cur_column[j].y_sp = pre_column[j].y_sp;
	  cur_column[j].x_sp_org = pre_column[j].x_sp_org;
	  cur_column[j].y_sp_org = pre_column[j].y_sp_org;
        } else {  // Set current cell as starting position
    	  cur_column[j].x_sp = i;
	  cur_column[j].y_sp = j;
    	  cur_column[j].x_sp_org = (i-1) * Tx + 1;
	  cur_column[j].y_sp_org = (j-1) * Ty + 1;
	}
      }

      // Compute DTW distance
      lx = i*Tx - cur_column[j].x_sp_org + 1;
      ly = j*Ty - cur_column[j].y_sp_org + 1;
      L = (double)(lx + ly) / 2;
      cur_column[j].dist = epsilon*L - cur_column[j].score;

      // Output information of optimal subsequence pair
      if (cur_column[j].x_sp == x_start && cur_column[j].y_sp == y_start) {
            printf ("%5ld  %5ld  %5d  %5d  %.10f\n", i, j, cur_column[j].x_sp,
		    cur_column[j].y_sp, cur_column[j].score);
      }


      /* Initialization for next computation */
      L = 0;
      lx = ly = 0;
      dist = 0;
      s_best = 0;
      s_diag = s_vert = s_hori = 0;

    }

    /* Judgement the end of processing */
    if (i > x_end) {
      exit(0);
    }


    /*
     * Initialization score matrix
     */
    diagpt = (i+1) * (double)sqd[1]->length/sqd[0]->length;
    band_sp = (int)floor(diagpt - scwidth);
    if (band_sp < 1)  band_sp = 1;

    for (j = band_sp-1; j >= band_sp-5; j--) {
      if (j == 0) {
	break;
      }
      pre_column[j].score = 0;
    }

    /*
     * Update matrix pointers
     */
    swap = cur_column;
    cur_column = pre_column;
    pre_column = swap;

  }

}


void error (char *toolname) {

  printf ("Usage: %s  'file of Sequence X'  'file of Sequence Y'  'Lmin'  'epsilon'  'scband width'  'sampling period of X'  'sampling period of Y'  'starting position of X'  'starting position of Y'  'end position of X'  'end position of Y'\n", toolname);

}



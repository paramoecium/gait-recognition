/*
 *  File name :   crossmatch.c
 */

#include	<stdio.h>
#include	<string.h>
#include	<stdlib.h>
#include	<time.h>
#include	<math.h>
#include	<sys/time.h>
#include	"crossmatch.h"

#define	BUFSIZE	256
#define Bd 1		// Weight for diagonal direction
#define Bv 0.5		// Weight for vertical direction
#define Bh 0.5		// Weight for horizontal direction
#define	LIST 100	// Upper limit value for candidate subsequence pairs in list arrays
#define EXPERIMENT 1    // 0: Performance evaluation, 1: Pattern discovery


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
int			listnum;			/* The number of candidate subsequence pairs in list arrays*/ 


int main(
  int	argc,
  char	*argv[]
) {

  char			dbfile[BUFSIZE], dbfile2[BUFSIZE];
  struct timeval	start, end;
  double		totaltime = 0.0;


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
#if EXPERIMENT==0
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
    cur_column[j].dist = 0;
    cur_column[j].score = 0;
    cur_column[j].x_sp = 0;
    cur_column[j].y_sp = 0;
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
  int		l;
  long int	band_sp, band_ep;
  int		lx, ly;	
  double	L;
  double	dist;
  double	s_best, s_diag, s_vert, s_hori;
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
      s_diag = Bd*epsilon + pre_column[j-1].score;
      s_vert = Bv*epsilon + cur_column[j-1].score;
      s_hori = Bh*epsilon + pre_column[j].score;
      s_best = max (&s_diag, &s_vert, &s_hori);

      cur_column[j].score = s_best - dist;

      // Reset score
      if (cur_column[j].score <= 0) {
        cur_column[j].score = 0;
    	cur_column[j].x_sp = i;
	cur_column[j].y_sp = j;
      // Keep score of neighboring cells
      } else {
        if (s_best == s_diag && pre_column[j-1].score != 0) {  // Select diagonal cell
    	  cur_column[j].x_sp = pre_column[j-1].x_sp;
	  cur_column[j].y_sp = pre_column[j-1].y_sp;
        } else if (s_best == s_vert && cur_column[j-1].score != 0) {  // Select vertical cell	
	  cur_column[j].x_sp = cur_column[j-1].x_sp;
	  cur_column[j].y_sp = cur_column[j-1].y_sp;
        } else if (s_best == s_hori && pre_column[j].score != 0) {  // Select horizontal cell
	  cur_column[j].x_sp = pre_column[j].x_sp;
	  cur_column[j].y_sp = pre_column[j].y_sp;
        } else {  // Set current cell as starting position
    	  cur_column[j].x_sp = i;
	  cur_column[j].y_sp = j;
	}
      }

      // Compute DTW distance
      lx = i - cur_column[j].x_sp + 1;
      ly = j - cur_column[j].y_sp + 1;
      L = (double)(lx + ly) / 2;
      cur_column[j].dist = epsilon*L - cur_column[j].score;


#if EXPERIMENT
      /*
       * Update candidate pairs in list arrays
       */

      if (cur_column[j].score >= epsilon*Lmin) {

	// Add subsequence pair as a new candidate
	if (listnum == 0) { 
	  list[listnum].flag = 0;
	  list[listnum].x_sp = cur_column[j].x_sp;
	  list[listnum].y_sp = cur_column[j].y_sp;
	  list[listnum].x_ep = i;
	  list[listnum].y_ep = j;
	  lx = list[listnum].x_ep - list[listnum].x_sp + 1;
	  ly = list[listnum].y_ep - list[listnum].y_sp + 1;
	  list[listnum].length = (double)(lx + ly)/2;
  	  list[listnum].score = cur_column[j].score;
	  list[listnum].dist = cur_column[j].dist;
	  listnum++;
	} else {
  	  for (l = 0; l < listnum; l++) {
	    if (list[l].x_sp == cur_column[j].x_sp
		&& list[l].y_sp == cur_column[j].y_sp) { 
 	      list[l].flag = 0;
	      // Overwrite maximum score
     	      if (cur_column[j].score-epsilon*Lmin >= list[l].score-epsilon*Lmin) {
		list[l].score = cur_column[j].score;
		list[l].dist = cur_column[j].dist;
	        list[l].x_ep = i;
	        list[l].y_ep = j;
		lx = list[l].x_ep - list[l].x_sp + 1;
		ly = list[l].y_ep - list[l].y_sp + 1;
		list[l].length = (double)(lx + ly)/2;
	      } 
	      break;
	    // Add subsequence pair as a new candidate
	    } else if (l == listnum-1) {
	      list[listnum].flag = 0;
	      list[listnum].x_sp = cur_column[j].x_sp;
	      list[listnum].y_sp = cur_column[j].y_sp;
	      list[listnum].x_ep = i;
	      list[listnum].y_ep = j;
	      lx = list[listnum].x_ep - list[listnum].x_sp + 1;
      	      ly = list[listnum].y_ep - list[listnum].y_sp + 1;
	      list[listnum].length = (double)(lx + ly)/2;
	      list[listnum].score = cur_column[j].score;
	      list[listnum].dist = cur_column[j].dist;
	      listnum++;
	      if (listnum >= LIST) {
	        printf ("The number of list is over %d\n", LIST);
	        exit(1);
	      }
	      break;
	    }
  	  }
        }

      // Processing for report
      } else if (listnum > 0) {
        // Check whether starting positions in current column correspond to the starting positions of candidate subsequences
    	for (l = 0; l < listnum; l++) {
	  if (list[l].flag == 1) {
	    if (list[l].x_sp == cur_column[j].x_sp
		&& list[l].y_sp == cur_column[j].y_sp) {
	      list[l].flag = 0;
	    }
	  }
  	}
      }
#endif

      /* Initialization for next computation */
      L = 0;
      lx = ly = 0;
      dist = 0;
      s_best = 0;
      s_diag = s_vert = s_hori = 0;

    }


#if EXPERIMENT
    /*
     * Report optimal subsequence pair 
     */
    if (listnum > 0) {
      for (l = 0; l < listnum; l++) {
	if (list[l].flag == 1) {
  	  printf ("X[%5d:%5d] \t Y[%5d:%5d] \t DTWdist: %.1f \t score: %.2f\n",
		  list[l].x_sp, list[l].x_ep, list[l].y_sp, list[l].y_ep, list[l].dist, list[l].score);
	  // Delete the reported pair in list arrays 
	  list[l].flag = list[listnum-1].flag;
	  list[l].x_sp = list[listnum-1].x_sp;
	  list[l].y_sp = list[listnum-1].y_sp;
	  list[l].x_ep = list[listnum-1].x_ep;
	  list[l].y_ep = list[listnum-1].y_ep;
	  list[l].length = list[listnum-1].length;
	  list[l].score = list[listnum-1].score;
	  list[l].dist = list[listnum-1].dist;
	  list[listnum-1].flag = 1;
	  list[listnum-1].x_sp = list[listnum-1].y_sp = 0;
	  list[listnum-1].x_ep = list[listnum-1].y_ep = 0;
	  list[listnum-1].length = 0;
	  list[listnum-1].score = 0;
	  list[listnum-1].dist = 0;
	  listnum--;
	  l--;
	} else {
	  list[l].flag = 1;
	}
      }
    }
#endif


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


#if EXPERIMENT
  /*
   * Report optimal subsequence pair in list arrays
   */
  for (l = 0; l < listnum; l++) {
    printf ("X[%5d:%5d] \t Y[%5d:%5d] \t DTWdist: %.1f \t score: %.2f\n",
 	    list[l].x_sp, list[l].x_ep, list[l].y_sp, list[l].y_ep, list[l].dist, list[l].score);
  }
#endif

}


void error (char *toolname) {

  printf ("Usage: %s  'file of Sequence X'  'file of Sequence Y'  'Lmin'  'epsilon'  'scband width'\n", toolname);

}



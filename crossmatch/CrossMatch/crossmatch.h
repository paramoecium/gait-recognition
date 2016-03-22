/*
 *  File name :   crossmatch.h
 */


#define	INITIALIZE_VALUE 1e20

extern void		malloc_double();
extern void		read_file();
extern double		element_dist();
extern long long int 	get_time_usec();
extern double		max();
extern double		min();


/* For data sequence */

struct DataElement {
  int		length;		/* Sequence length */
  double	*value;		/* Sequence value */
};


/* For score and position matrices */

struct MatrixElement {
  double	dist;		/* DTW distance */
  double	score;		/* Score */
  int		x_sp;		/* Starting position of X */
  int		y_sp;		/* Starting position of Y */
};


/* For list arrays, which store candidate subsequence pairs */

struct SeqList {
  int		flag;		/* Flag for report */
  int		x_sp;		/* Starting position of X */
  int		y_sp;		/* Starting position of Y */
  int		x_ep;		/* End position of X */
  int		y_ep;		/* End position of Y */
  double	dist;		/* DTW distance */
  double	score;		/* Score */
  double	length;		/* Length of subsquence pair */
};



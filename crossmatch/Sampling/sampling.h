/*
 *  File name :   sampling.h
 */


#define	INITIALIZE_VALUE 1e20

extern void		malloc_double();
extern void		read_file();
extern double		element_dist();
extern long long int 	get_time_usec();
extern double		gettimeofday_sec();
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
  int		x_sp;		/* Starting position of sampled sequence X */
  int		y_sp;		/* Starting position of sampled sequence Y */
  int		x_sp_org;	/* Starting position of original sequence X */
  int		y_sp_org;	/* Starting position of original sequence Y */
};


/* For list arrays, which store candidate subsequence pairs */

struct SeqList {
  int		flag;		/* Flag for report */
  int		x_sp;		/* Starting position of sampled sequence X */
  int		y_sp;		/* Starting position of sampled sequence Y */
  int		x_ep;		/* End position of sampled sequence X */
  int		y_ep;		/* End position of sampled sequence Y */
  int		x_sp_org;	/* Starting position of original sequence X */
  int		y_sp_org;	/* Starting position of original sequence Y */
  int		x_ep_org;	/* End position of original sequence of X */
  int		y_ep_org;	/* End position of original sequence of Y */
  double	dist;		/* DTW distance */
  double	score;		/* Score */
  double	length;		/* Length of subsquence pair */
};

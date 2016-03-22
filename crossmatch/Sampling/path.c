/*
 * File name : path.c
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define BUFSIZE 256
#define BYTE 1048576
#define Bd 1
#define Bv 0.5
#define Bh 0.5

extern int	count_lines();
extern void	error();

struct Matrix {
  int		x_start;
  int		y_start;
  int		x_end;
  int		y_end;
  double	score;
};


int main (int argc, char * argv[])
{

  int		i, j;
  FILE		*fpin, *fpin2, *fpread;	
  FILE		*fpout, *fpwrite;
  char		filename[BUFSIZE];
  int		filenum;
  char		**readfile;
  char		writefile[BUFSIZE];
  int		*x_start, *y_start;
  int		*x_end, *y_end;	

  char		data[BYTE];
  size_t	size;
  int		count;
  off_t		loc, loc2;
  int		set;

  struct Matrix	**matrix;
  struct Matrix	*swap;
  double	epsilon;
  int		str_len;
  int		x_ep, y_ep;
  int		nxt_x_ep, nxt_y_ep;
  int		init;	
  int		loop;
  off_t		loc3;
  int		judge;
  int		match;
  int		mem_rd;
  int		x_sp_tmp, y_sp_tmp, pre_x_ep_tmp, x_ep_tmp, y_ep_tmp, max_x, max_y;
  double	s_tmp, max_s;


  for (i = 0; i < BYTE; i++) {
    data[BYTE] = '\0';
  }


  if (argc != 4) {
    error (argv[0]);
  }
  strcpy (filename, argv[1]);
  epsilon = atof(argv[2]);
  str_len = atoi(argv[3]);
  filenum = count_lines(filename);


  readfile = (char **)malloc(filenum * sizeof(char*));
  for (i = 0; i < filenum; i++) {
    readfile[i] = (char *)malloc(BUFSIZE * sizeof(char));
  }
  x_start = (int *)malloc(filenum * sizeof(int));
  y_start = (int *)malloc(filenum * sizeof(int));
  x_end = (int *)malloc(filenum * sizeof(int));
  y_end = (int *)malloc(filenum * sizeof(int));

  matrix = (struct Matrix **)malloc(2 * sizeof(struct Matrix*));
  matrix[0] = (struct Matrix *)malloc(str_len * sizeof(struct Matrix));
  matrix[1] = (struct Matrix *)malloc(str_len * sizeof(struct Matrix));
  for (i = 0; i < 2; i++) {
    for (j = 0; j < str_len; j++) {
      matrix[i][j].x_start = 0;
      matrix[i][j].y_start = 0;
      matrix[i][j].x_end = 0;
      matrix[i][j].y_end = 0;
      matrix[i][j].score = 0;
    }
  }


  if ((fpread = fopen64 (filename, "r")) == NULL) {
    printf ("%s cannot open listfile\n", filename);
    exit(1);
  }
  i = 0;
  for (j = 0; j < filenum; j++) {
    fscanf (fpread, "%s %s %d %d %d %d\n", 
	readfile[i], writefile, &x_end[i], &y_end[i], &x_start[i], &y_start[i]);
    i++;
  }
  fclose (fpread);


  if ((fpwrite = fopen64 (writefile, "w")) == NULL) {
    printf ("%s (write) cannot open\n", writefile);
    exit(1);
  } 


  for (i = 0; i < filenum; i++) {

    /*
     * Reverse contents of the file
     */
    printf ("#%d reverse phase\n", i+1);

    set = 0;

    if ((fpin = fopen64 (readfile[i], "r")) == NULL) {
      printf ("%s cannot open\n", readfile[i]);
      exit(1);
    }
    if ((fpout = fopen64 ("temp.txt", "w")) == NULL) {
      printf ("temp.txt cannot open\n");
      exit(1);
    }

    /* Move pointer to the end of original file */
    fseeko (fpin, 0L, SEEK_END);
    loc = ftello(fpin) - BYTE;
    loc2 = ftello(fpin);

    if (loc2 == 0) {
      printf ("%s is empty\n", readfile[i]);
      goto last;
    }

    while (loc2 != 0) {

      if (loc < 0) {
        loc = 0;
        set = 1;
      }

      fseeko (fpin, loc, SEEK_SET);
      if (loc == 0) {
	size = fread (data, 1, loc2, fpin);
      } else {
        size = fread (data, 1, BYTE, fpin);
      }

      count = 0;

      for (j = size-1; j >= 0; j--) {
        if (set == 1) {
  	  if (j == 0) {
	    fwrite (&data[j], count, 1, fpout);
	    count = 0;
	    break;
  	  } else if (data[j] == '\n' && j != size-1) {
	    fwrite (&data[j+1], count, 1, fpout);
	    count = 0;
  	  }
          count++;
        } else {
          if (data[j] == '\n' && j != size-1) {
  	    fwrite (&data[j+1], count, 1, fpout);
	    count = 0;
	  }
          count++;
        }
      }

      loc2 = loc + count;
      loc = loc2 - BYTE;

      for (j = 0; j < BYTE; j++) {
        data[j] = '\0';
      }

    }

    fclose (fpin);
    fclose (fpout);

    /*
     * Trace back warping path 
     */
    printf ("#%d path detection phase\n", i+1);

    x_ep = x_end[i];
    y_ep = y_end[i];
    nxt_x_ep = -1;
    nxt_y_ep = -1;
    init = 0;
    match = 0;
    mem_rd = 1;

    fprintf (fpwrite, "%5d %5d\n", x_ep, y_ep);

    if ((fpin2 = fopen64 ("temp.txt", "r")) == NULL) {
      printf ("temp.txt[%d] cannot open\n", i);
      exit(1);
    }

    while (x_ep != x_start[i] || y_ep != y_start[i]) {

      /* 
       * Find end point for path search from matrix elements
       */
      if (mem_rd == 1) {

        if (init == 0) {
  	  loop = 0;
          init = 1;
        } else {
  	  loop = 1;
        }
        pre_x_ep_tmp = -1;

        for (j = loop; j < 2; j++) {

  	  loc3 = ftello(fpin2);
	  judge = fscanf (fpin2, "%d %d %d %d %lf\n", 
			&x_ep_tmp, &y_ep_tmp, &x_sp_tmp, &y_sp_tmp, &s_tmp);

	  if (match == 0) {
	    if (x_ep_tmp == x_ep) {
	      match = 1;
	    } else {
	      break;
	    }
  	  }

	  if (match == 1 && loop != 0) {
	    match = 2;
  	  }

	  while (judge != EOF) {

	    if (x_ep_tmp == pre_x_ep_tmp || pre_x_ep_tmp == -1) {
	      matrix[j][y_ep_tmp].x_start = x_sp_tmp;
	      matrix[j][y_ep_tmp].y_start = y_sp_tmp;
	      matrix[j][y_ep_tmp].x_end = x_ep_tmp;
	      matrix[j][y_ep_tmp].y_end = y_ep_tmp;
	      matrix[j][y_ep_tmp].score = s_tmp;
	      pre_x_ep_tmp = x_ep_tmp;

	    } else { 
	      fseeko (fpin2, loc3, SEEK_SET);
	      break;
	    }
	    loc3 = ftello(fpin2);
	    judge = fscanf (fpin2, "%d %d %d %d %lf\n", 
			&x_ep_tmp, &y_ep_tmp, &x_sp_tmp, &y_sp_tmp, &s_tmp);
  	  }

        }

      }	

      /* 
       * Search next point on warping path 
       */
      if (match == 2) {

        if (Bd*epsilon+matrix[1][y_ep-1].score == Bd*epsilon 
	     && Bv*epsilon+matrix[0][y_ep-1].score == Bv*epsilon) {
	  if (Bh*epsilon+matrix[1][y_ep].score == Bh*epsilon) { 
	    printf ("error1: #%d fail to detect warping path.\n", i+1);
	    mem_rd = 0;
	    break;
	  } else {
	    max_x = matrix[1][y_ep].x_end;
	    max_y = matrix[1][y_ep].y_end;
	    max_s = matrix[1][y_ep].score;
	    mem_rd = 1;
	  }
        } else if (Bd*epsilon+matrix[1][y_ep-1].score != Bd*epsilon 
	     && (Bd*epsilon+matrix[1][y_ep-1].score >= Bv*epsilon+matrix[0][y_ep-1].score
	         || Bv*epsilon+matrix[0][y_ep-1].score == Bv*epsilon)) { 
  	  max_x = matrix[1][y_ep-1].x_end;
	  max_y = matrix[1][y_ep-1].y_end;
	  max_s = matrix[1][y_ep-1].score;
	  mem_rd = 1;
          if (Bh*epsilon+matrix[1][y_ep].score != Bh*epsilon 
	       && Bh*epsilon+matrix[1][y_ep].score > Bd*epsilon+max_s) {
	    max_x = matrix[1][y_ep].x_end;
	    max_y = matrix[1][y_ep].y_end;
	    max_s = matrix[1][y_ep].score;
	    mem_rd = 1;
          }
        } else if (Bv*epsilon+matrix[0][y_ep-1].score != Bv*epsilon
		    &&(Bd*epsilon+matrix[1][y_ep-1].score < Bv*epsilon+matrix[0][y_ep-1].score
		       || Bd*epsilon+matrix[1][y_ep-1].score == Bd*epsilon)) {
  	  max_x = matrix[0][y_ep-1].x_end;
	  max_y = matrix[0][y_ep-1].y_end;
	  max_s = matrix[0][y_ep-1].score;
	  mem_rd = 0;
          if (Bh*epsilon+matrix[1][y_ep].score != Bh*epsilon 
	       && Bh*epsilon+matrix[1][y_ep].score > Bv*epsilon+max_s) {
	    max_x = matrix[1][y_ep].x_end;
	    max_y = matrix[1][y_ep].y_end;
	    max_s = matrix[1][y_ep].score;
	    mem_rd = 1;
          }
        }

        nxt_x_ep = max_x;
        nxt_y_ep = max_y;
        max_x = max_y = 0;
        max_s = 0;

        if ((nxt_x_ep == 0 && nxt_y_ep == 0 && x_start[i] != 0 && y_start[i] != 0)
	     || (judge == EOF && mem_rd == 1 && nxt_x_ep != x_start[i] && nxt_y_ep != y_start[i])) {
	  printf ("error2: #%d fail to detect warping path.\n", i+1);
	  mem_rd = 1;
	  break;
        }

        x_ep = nxt_x_ep;
        y_ep = nxt_y_ep;
        fprintf (fpwrite, "%5d %5d\n", x_ep, y_ep);

      }

      /* Update pointers for next computation */
      if (mem_rd == 1) {
        swap = matrix[0];
        matrix[0] = matrix[1];
        matrix[1] = swap;
        for (j = 0; j < str_len; j++) {
          matrix[1][j].x_start = 0;
          matrix[1][j].y_start = 0;
          matrix[1][j].x_end = 0;
          matrix[1][j].y_end = 0;
          matrix[1][j].score = 0;
        }
      }

    }

    fclose(fpin2);
    last:;

  }


  free(readfile);
  free(x_start);
  free(y_start);
  free(x_end);
  free(y_end);
  free(matrix[0]);
  free(matrix[1]);
  free(matrix);

  return 0;

}


/*************************************************************
 * Function name	: count_lines
 * Input		: char *filename -  File name
 * note			: Count the number of lines in file
 *************************************************************/
int count_lines (char *filename) {

  int	i = 0;
  char	read_buf[BUFSIZE];
  FILE	*fp;

  if ((fp = fopen (filename, "r")) == NULL) {
    printf ("%s cannot open (count_lines function)\n", filename);
    exit(1);
  }
  while (fgets (read_buf, BUFSIZE, fp) != NULL) {
    if (strcmp (read_buf, "\n") == 0) {
      break;
    }
    i++;
  }
  fclose(fp);
  return i;

}


void error (char *toolname) {

  printf ("Usage: %s  'file name'  'distance threshold'  'length of seq. Y'\n", toolname);

}

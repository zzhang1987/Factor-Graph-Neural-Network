typedef struct {
  int N ;
  int K ; 
  int **C ;
} ulist_matrix ;

typedef struct {
  int N ;
  int u ; /* whether this is an upper triangular matrix, and its inverse, 
	     or vice versa */
  int i ; /* whether the inverse has been computed */
  unsigned char **m ;
  unsigned char **mi ; /* for the inverse */
} triangular_matrix ;

typedef struct {
  int style ;
  double density ; 
  double *f ;
} te_params ; 

typedef struct {
  int N , M ;
  int **mlist;
  int **nlist;
  int *num_mlist;
  int *num_nlist;
  int *l_up_to ;
  int *norder ;
  int biggest_num_m ;       /* actual biggest sizes */
  int biggest_num_n ; 
  int biggest_num_m_alloc ; /* sizes used for memory allocation */
  int biggest_num_n_alloc ; 
  int tot ; 
  int same_length ;  /* whether all vectors in mlist and nlist have same length */
} alist_matrix ;

void report_alist ( alist_matrix * , int  ) ;
void alist_times_cvector_sparse_mod2  ( alist_matrix * , unsigned char * , unsigned char * ) ;
void alist_transpose_cvector_sparse_mod2  ( alist_matrix * , unsigned char * , unsigned char * ) ;
void alist_times_cvector_mod2 ( alist_matrix * , unsigned char * , unsigned char * ) ;
void alist_times_ivector_sparse
( alist_matrix *a , int *x , int *y ) ;
void alist_transpose_ivector_sparse
( alist_matrix *a ,  int *x , int *y ) ;
int finish_off_alist ( alist_matrix * ) ;
void add_to_alist ( alist_matrix * , int n , int  ) ;
void subtract_from_alist ( alist_matrix * , int  , int  ) ;
void list_slide_back ( int * , int , int ) ;
void subtract_col_from_alist ( alist_matrix * , int  ) ;
void kill_blank_col_from_alist ( alist_matrix * , int ) ;
int code0_matrix_alist ( alist_matrix * , int  , int  , int  ) ;
int cmatrix_2_alist ( unsigned char ** , alist_matrix * , int ) ;
int gen_cmatrix_2_alist ( unsigned char ** , alist_matrix * , int , int , int ) ;
int cmatrices_2_alist ( unsigned char ** , unsigned char ** ,
			       alist_matrix * , int , int ) ;
int cmatrix_andI_2_alist ( unsigned char ** , 
			       alist_matrix * , int , int ) ;
void initialize_alist ( alist_matrix * , int , int , int , int ) ;
void free_alist ( alist_matrix *  ) ;
unsigned char cdotprod_mod2_index (  int * , unsigned char * , int  , int  ) ;
int cdotprod_index (  int * , unsigned char * , int , int ) ;
void write_alist ( FILE * , alist_matrix * ) ;
void write_alist_transpose ( FILE * , alist_matrix * ) ;
int read_allocate_alist ( alist_matrix * , char * ) ;
void sparse_random_cmatrix ( unsigned char ** , int , int ) ;
void sparse_random_cmatrix2 ( unsigned char ** , int , int ) ;
void sparse_invertible_cmatrix (  unsigned char **, unsigned char **
				, int , int  ) ;
void sparse_rectangular_cmatrix ( unsigned char ** , int , int , int , int ) ;
int sparse_rectangular_cmatrix2 ( unsigned char ** , int , int , int , int , int ) ;
int 	  read_int_from_dcum ( double * , int , int  , double  ) ; 
int 	  read_int_from_cum ( int * , int , int  , double  ) ; 
int old_sparse_rectangular_alist ( alist_matrix * , int , int , int , int , int , int ) ;
int sparse_rectangular_alist ( alist_matrix * , int , int , int , int , int , int ) ;
void invert_left_hand_end_cmatrix ( unsigned char ** , unsigned char **, int) ;
int fixed_wt_cvector ( unsigned char *v , int n , int lo , int hi ) ;
int random_cvector ( unsigned char *v , double f , int lo , int hi ) ;
int uneven_sparse_rectangular_alist ( alist_matrix * , int  , int  , 
				     double  , int  , int , int ) ;
int staircase_alist ( alist_matrix * , int  , int  , 
				     int , int  , int  ) ;
int slope_alist ( alist_matrix * , int  , int  , 
				     int , int  , int  ) ;

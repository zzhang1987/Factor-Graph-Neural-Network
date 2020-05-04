extern int ran3argument=1 ;
#define ran_seed(s) ran3(-(int)(s))
#define rann()      ran3(&ran3argument)
#define ranu()      ran3(&ran3argument)
#define ranf()      ran3(&ran3argument)

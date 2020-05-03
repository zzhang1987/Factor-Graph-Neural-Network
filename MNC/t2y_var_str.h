  int verbose ;     /* verbosity            <0>    */
  char rfile[100] ; /* file: y vector       <->    */
  int N ;           /* block length         <->    */
  int tsuffix ;     /* whether tfile should have suffices .0001,.0002... appended <0>    */
  int ysuffix ;     /* whether yfile should have suffices .0001,.0002... appended <0>    */
  char tfile[100] ; /* t vector             <->    */
  double gcx ;      /* signal_to_noise ratio <1.0>  */
  double sigma;     /* std of the bursty noise */
  double rho;       /* scale of the bursty noise */
  long int seed ;   /* -                    <1234> */

  int verbose ;     /* verbosity            <0>    */
  char rfile[100] ; /* file: r vector       <->    */
  int N ;           /* block length         <->    */
  int tsuffix ;     /* whether tfile should have suffices .0001,.0002... appended <0>    */
  int rsuffix ;     /* whether rfile should have suffices .0001,.0002... appended <0>    */
  char tfile[100] ; /* t vector             <->    */
  double fn ;       /* noise density        <0.1>  */
  long int seed ;   /* -                    <1234> */

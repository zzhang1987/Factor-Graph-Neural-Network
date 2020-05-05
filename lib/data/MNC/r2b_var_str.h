  int verbose ;     /* verbosity            <0>    */
  char rfile[100] ; /* file: r vector       <->    */
  int N ;           /* number of bits in block <->    */
  int rsuffix ;     /* rfile: add suffices .0001,.0002... <0>    */
  int bsuffix ;     /* whether bfile should have suffices .0001,.0002... appended <0>    */
  char bfile[100] ; /* b vector             <->    */
  double fn ;       /* bias of genuine noise bits <0.1>  */

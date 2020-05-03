  int verbose ;     /* verbosity            <0>    */
  char bfile[100] ; /* file: bias vector    <->    */
  int bnfromfile ;  /* -                    <1>    */
  int bsfromfile ;  /* -                    <1>    */
  int bfromfile ;   /* whether some of b comes from file <1>    */
  double bfn ;      /* constant value for bias of noise bits <0.1>  */
  double bfs ;      /* constant value for bias of signal bits <0.5>  */
  int K ;           /* number of source bits (NMN) <->    */
  int N ;           /* number of noise bits (NMN) <->    */
  char zfile[100] ; /* file: z vector       <->    */
  int zfromfile ;   /* -                    <0>    */
  int zfixed ;      /* constant value for z <0>    */
  char Afile[100] ; /* file: A matrix       <->    */
  int bsuffix ;     /* whether bfile should have suffices .0001,.0002... appended <0>    */
  int zsuffix ;     /* whether zfile should have suffices .0001,.0002... appended <0>    */
  int xsuffix ;     /* whether xfile should have suffices .0001,.0002... appended <0>    */
  char xfile[100] ; /* x vector             <->    */
  int xtofile ;     /* -                    <1>    */
  int xsourceonly ; /* write decoded source bits only <0>    */

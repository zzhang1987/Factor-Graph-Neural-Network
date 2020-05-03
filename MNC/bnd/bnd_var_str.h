  int loops ;       /* Max number of loops  <20>   */
  int loop ;        /* -                    <0>    */
  int writelog ;    /* number of states to log <0>    */
  char logfile[100] ;/* for state info       <->    */
  int doclip ;      /* whether to clip extreme probablities <1>    */
  double clip ;     /* clip probs here      <0.9999999999> */
  double tinydiv ;  /* value below which division should not be done <1e-40> */
  int dofudge ;     /* whether to fudge q values up/down or not <0>    */
  double fudge ;    /* fudge scale          <1.0>  */

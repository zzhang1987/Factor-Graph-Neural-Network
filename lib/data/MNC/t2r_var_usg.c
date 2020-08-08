  DNT; fprintf( fp, "-v verbose <%d>                (verbosity           )", c->verbose);
  DNT; fprintf( fp, "-rfile rfile                  (file: r vector      )");
  DNT; fprintf( fp, "-n n                          (block length        )");
  DNT; fprintf( fp, "-tsuffix ts <%d>               (whether tfile should have suffices .0001,.0002... appended)", c->tsuffix);
  DNT; fprintf( fp, "-rsuffix rs <%d>               (whether rfile should have suffices .0001,.0002... appended)", c->rsuffix);
  DNT; fprintf( fp, "-tfile tfile                  (t vector            )");
  DNT; fprintf( fp, "-fn fn <%9.3g>         (noise density       )", c->fn);
  DNT; fprintf( fp, "-seed seed <%ld>            (-                   )", c->seed);

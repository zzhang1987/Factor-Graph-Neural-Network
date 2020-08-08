  DNT; fprintf( fp, "-v verbose <%d>                (verbosity           )", c->verbose);
  DNT; fprintf( fp, "-yfile yfile                  (file: y vector      )");
  DNT; fprintf( fp, "-n n                          (block length        )");
  DNT; fprintf( fp, "-tsuffix ts <%d>               (whether tfile should have suffices .0001,.0002... appended)", c->tsuffix);
  DNT; fprintf( fp, "-ysuffix ys <%d>               (whether yfile should have suffices .0001,.0002... appended)", c->ysuffix);
  DNT; fprintf( fp, "-tfile tfile                  (t vector            )");
  DNT; fprintf( fp, "-gcx gcx <%9.3g>       (signal_to_noise ratio)", c->gcx);
  DNT; fprintf( fp, "-seed seed <%ld>            (-                   )", c->seed);

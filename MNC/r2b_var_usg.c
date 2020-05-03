  DNT; fprintf( fp, "-v verbose <%d>                (verbosity           )", c->verbose);
  DNT; fprintf( fp, "-rfile rfile                  (file: r vector      )");
  DNT; fprintf( fp, "-n n                          (number of bits in block)");
  DNT; fprintf( fp, "-rsuffix rs <%d>               (rfile: add suffices .0001,.0002...)", c->rsuffix);
  DNT; fprintf( fp, "-bsuffix bs <%d>               (whether bfile should have suffices .0001,.0002... appended)", c->bsuffix);
  DNT; fprintf( fp, "-bfile bfile                  (b vector            )");
  DNT; fprintf( fp, "-fn fn <%9.3g>         (bias of genuine noise bits)", c->fn);

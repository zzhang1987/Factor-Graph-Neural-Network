    else if ( strcmp (argv[i], "-v") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->verbose));
       }
    }
    else if ( strcmp (argv[i], "-yfile") == 0 ) {
       if ( i + 1 == argc ) { ERROR1; }
       else {
         strcpy(c->yfile, argv[++i]);
       }
    }
    else if ( strcmp (argv[i], "-n") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->N));
       }
    }
    else if ( strcmp (argv[i], "-bsuffix") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->bsuffix));
       }
    }
    else if ( strcmp (argv[i], "-ysuffix") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->ysuffix));
       }
    }
    else if ( strcmp (argv[i], "-bfile") == 0 ) {
       if ( i + 1 == argc ) { ERROR1; }
       else {
         strcpy(c->bfile, argv[++i]);
       }
    }
    else if ( strcmp (argv[i], "-gcx") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%lf", &(c->gcx));
       }
    }

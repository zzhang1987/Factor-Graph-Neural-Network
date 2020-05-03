    else if ( strcmp (argv[i], "-v") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->verbose));
       }
    }
    else if ( strcmp (argv[i], "-yfile") == 0 ) {
       if ( i + 1 == argc ) { ERROR1; }
       else {
         strcpy(c->rfile, argv[++i]);
       }
    }
    else if ( strcmp (argv[i], "-n") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->N));
       }
    }
    else if ( strcmp (argv[i], "-tsuffix") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->tsuffix));
       }
    }
    else if ( strcmp (argv[i], "-ysuffix") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->ysuffix));
       }
    }
    else if ( strcmp (argv[i], "-tfile") == 0 ) {
       if ( i + 1 == argc ) { ERROR1; }
       else {
         strcpy(c->tfile, argv[++i]);
       }
    }
    else if ( strcmp (argv[i], "-gcx") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%lf", &(c->gcx));
       }
    }
    else if ( strcmp (argv[i], "-seed") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%ld", &(c->seed));
       }
    }
    else if ( strcmp (argv[i], "-rho") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%ld", &(c->rho));
       }
    }
    else if ( strcmp (argv[i], "-sigma") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%lf", &(c->sigma));
        printf("%.6f\n", c->sigma);
       }
    }

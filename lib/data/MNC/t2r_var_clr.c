    else if ( strcmp (argv[i], "-v") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->verbose));
       }
    }
    else if ( strcmp (argv[i], "-rfile") == 0 ) {
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
    else if ( strcmp (argv[i], "-rsuffix") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->rsuffix));
       }
    }
    else if ( strcmp (argv[i], "-tfile") == 0 ) {
       if ( i + 1 == argc ) { ERROR1; }
       else {
         strcpy(c->tfile, argv[++i]);
       }
    }
    else if ( strcmp (argv[i], "-fn") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%lf", &(c->fn));
       }
    }
    else if ( strcmp (argv[i], "-seed") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%ld", &(c->seed));
       }
    }

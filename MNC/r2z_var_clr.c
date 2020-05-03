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
    else if ( strcmp (argv[i], "-Cnfile") == 0 ) {
       if ( i + 1 == argc ) { ERROR1; }
       else {
         strcpy(c->Cnfile, argv[++i]);
       }
    }
    else if ( strcmp (argv[i], "-rsuffix") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->rsuffix));
       }
    }
    else if ( strcmp (argv[i], "-zsuffix") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->zsuffix));
       }
    }
    else if ( strcmp (argv[i], "-zfile") == 0 ) {
       if ( i + 1 == argc ) { ERROR1; }
       else {
         strcpy(c->zfile, argv[++i]);
       }
    }

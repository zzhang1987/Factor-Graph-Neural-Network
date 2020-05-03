    else if ( strcmp (argv[i], "-v") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->verbose));
       }
    }
    else if ( strcmp (argv[i], "-sfile") == 0 ) {
       if ( i + 1 == argc ) { ERROR1; }
       else {
         strcpy(c->sfile, argv[++i]);
       }
    }
    else if ( strcmp (argv[i], "-k") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->K));
       }
    }
    else if ( strcmp (argv[i], "-n") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->N));
       }
    }
    else if ( strcmp (argv[i], "-Gfile") == 0 ) {
       if ( i + 1 == argc ) { ERROR1; }
       else {
         strcpy(c->Gfile, argv[++i]);
       }
    }
    else if ( strcmp (argv[i], "-tsuffix") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->tsuffix));
       }
    }
    else if ( strcmp (argv[i], "-tfile") == 0 ) {
       if ( i + 1 == argc ) { ERROR1; }
       else {
         strcpy(c->tfile, argv[++i]);
       }
    }
    else if ( strcmp (argv[i], "-smn") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->smn));
       }
    }

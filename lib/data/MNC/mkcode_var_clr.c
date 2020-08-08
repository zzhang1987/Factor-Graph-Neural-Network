    else if ( strcmp (argv[i], "-v") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->verbose));
       }
    }
    else if ( strcmp (argv[i], "-Ain") == 0 ) {
       if ( i + 1 == argc ) { ERROR1; }
       else {
         strcpy(c->Ain, argv[++i]);
       }
    }
    else if ( strcmp (argv[i], "-Aout") == 0 ) {
       if ( i + 1 == argc ) { ERROR1; }
       else {
         strcpy(c->Aout, argv[++i]);
       }
    }
    else if ( strcmp (argv[i], "-G") == 0 ) {
       if ( i + 1 == argc ) { ERROR1; }
       else {
         strcpy(c->G, argv[++i]);
       }
    }
    else if ( strcmp (argv[i], "-Cn") == 0 ) {
       if ( i + 1 == argc ) { ERROR1; }
       else {
         strcpy(c->Cn, argv[++i]);
       }
    }
    else if ( strcmp (argv[i], "-force") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->force));
       }
    }

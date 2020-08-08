    else if ( strcmp (argv[i], "-v") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->verbose));
       }
    }
    else if ( strcmp (argv[i], "-bfile") == 0 ) {
       if ( i + 1 == argc ) { ERROR1; }
       else {
         strcpy(c->bfile, argv[++i]);
         c->bnfromfile=1;
       }
    }
    else if ( strcmp (argv[i], "-bfn") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%lf", &(c->bfn));
         c->bnfromfile=0;
       }
    }
    else if ( strcmp (argv[i], "-bfs") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%lf", &(c->bfs));
         c->bsfromfile=0;
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
    else if ( strcmp (argv[i], "-zfile") == 0 ) {
       if ( i + 1 == argc ) { ERROR1; }
       else {
         strcpy(c->zfile, argv[++i]);
         c->zfromfile=1;
       }
    }
    else if ( strcmp (argv[i], "-zfixed") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->zfixed));
         c->zfromfile=0;
       }
    }
    else if ( strcmp (argv[i], "-Afile") == 0 ) {
       if ( i + 1 == argc ) { ERROR1; }
       else {
         strcpy(c->Afile, argv[++i]);
       }
    }
    else if ( strcmp (argv[i], "-bsuffix") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->bsuffix));
       }
    }
    else if ( strcmp (argv[i], "-zsuffix") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->zsuffix));
       }
    }
    else if ( strcmp (argv[i], "-xsuffix") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->xsuffix));
       }
    }
    else if ( strcmp (argv[i], "-xfile") == 0 ) {
       if ( i + 1 == argc ) { ERROR1; }
       else {
         strcpy(c->xfile, argv[++i]);
         c->xtofile=1;
       }
    }
    else if ( strcmp (argv[i], "-xso") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(c->xsourceonly));
       }
    }

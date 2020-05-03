    else if ( strcmp (argv[i], "-bndloops") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(bndc->loops));
       }
    }
    else if ( strcmp (argv[i], "-bndlogn") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(bndc->writelog));
       }
    }
    else if ( strcmp (argv[i], "-bndlog") == 0 ) {
       if ( i + 1 == argc ) { ERROR1; }
       else {
         strcpy(bndc->logfile, argv[++i]);
         if(!bndc->writelog)bndc->writelog=10;
       }
    }
    else if ( strcmp (argv[i], "-doclip") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(bndc->doclip));
       }
    }
    else if ( strcmp (argv[i], "-clip") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%lf", &(bndc->clip));
         bndc->doclip=1;
       }
    }
    else if ( strcmp (argv[i], "-tinydiv") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%lf", &(bndc->tinydiv));
       }
    }
    else if ( strcmp (argv[i], "-dofudge") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%d", &(bndc->dofudge));
       }
    }
    else if ( strcmp (argv[i], "-fudge") == 0 ) {
      if ( i + 1 == argc ) { ERROR1; }
      else {
        cs *= sscanf(argv[++i], "%lf", &(bndc->fudge));
         bndc->dofudge=1;
       }
    }

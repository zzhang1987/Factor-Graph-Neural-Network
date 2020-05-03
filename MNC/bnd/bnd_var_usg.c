  DNT; fprintf( fp, "-bndloops loops <%d>          (Max number of loops )", bndc->loops);
  DNT; fprintf( fp, "-bndlogn n <%d>                (number of states to log)", bndc->writelog);
  DNT; fprintf( fp, "-bndlog logfile               (for state info      )");
  DNT; fprintf( fp, "-doclip doclip <%d>            (whether to clip extreme probablities)", bndc->doclip);
  DNT; fprintf( fp, "-clip clip <%9.3g>     (clip probs here     )", bndc->clip);
  DNT; fprintf( fp, "-tinydiv tiny <%9.3g>  (value below which division should not be done)", bndc->tinydiv);
  DNT; fprintf( fp, "-dofudge dofudge <%d>          (whether to fudge q values up/down or not)", bndc->dofudge);
  DNT; fprintf( fp, "-fudge fudge <%9.3g>   (fudge scale         )", bndc->fudge);

PRO convert_csv, date, sc, bflag

  nevents=5

  read_input, date, nevents, sc, bflag=bflag
  track_statistics_maike, date, sc, /make_dir, /save_file, /no_cor, bflag=bflag

END

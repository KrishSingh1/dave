
run 10834:

- initial implementation of gzip output copied from kutabiri svn branch
- gzip output enabled
- several daily patch-level outputs enabled during spinup
- 71% of simulation spent in outdaily(), 54% of simulation spent in
  GZFileOutputChannel::finish_row(), 38% of simulation spent in xtring::printf()
  as called from this finish_row() function

run 14775:

- improved implementation of gzip output involving repeated calls to gzprintf(),
  to eliminate the time wasted in xtring::printf().

run 16186:

- same as 14775 with plaintext output

run 19135:

- hoping to eliminate the overhead associated with std::unique_ptr<OutputFile>

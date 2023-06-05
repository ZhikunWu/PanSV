
```
$ pairtools dedup --help
Usage: pairtools dedup [OPTIONS] [PAIRS_PATH]

  Find and remove PCR/optical duplicates.

  Find PCR duplicates in an upper-triangular flipped sorted pairs/pairsam
  file. Allow for a +/-N bp mismatch at each side of duplicated molecules.

  PAIRS_PATH : input triu-flipped sorted .pairs or .pairsam file.  If the path
  ends with .gz/.lz4, the input is decompressed by pbgzip/lz4c.  By default,
  the input is read from stdin.

Options:
  -o, --output TEXT               output file for pairs after duplicate
                                  removal. If the path ends with .gz or .lz4,
                                  the output is pbgzip-/lz4c-compressed. By
                                  default, the output is printed into stdout.
  --output-dups TEXT              output file for duplicated pairs.  If the
                                  path ends with .gz or .lz4, the output is
                                  pbgzip-/lz4c-compressed. If the path is the
                                  same as in --output or -, output duplicates
                                  together  with deduped pairs. By default,
                                  duplicates are dropped.
  --output-unmapped TEXT          output file for unmapped pairs. If the path
                                  ends with .gz or .lz4, the output is
                                  pbgzip-/lz4c-compressed. If the path is the
                                  same as in --output or -, output unmapped
                                  pairs together with deduped pairs. If the
                                  path is the same as --output-dups, output
                                  unmapped reads together with dups. By
                                  default, unmapped pairs are dropped.
  --output-stats TEXT             output file for duplicate statistics.  If
                                  file exists, it will be open in the append
                                  mode. If the path ends with .gz or .lz4, the
                                  output is pbgzip-/lz4c-compressed. By
                                  default, statistics are not printed.
  --max-mismatch INTEGER          Pairs with both sides mapped within this
                                  distance (bp) from each other are considered
                                  duplicates.  [default: 3]
  --method [max|sum]              define the mismatch as either the max or the
                                  sum of the mismatches ofthe genomic
                                  locations of the both sides of the two
                                  compared molecules  [default: max]
  --sep TEXT                      Separator (\t, \v, etc. characters are
                                  supported, pass them in quotes)
  --comment-char TEXT             The first character of comment lines
  --send-header-to [dups|dedup|both|none]
                                  Which of the outputs should receive header
                                  and comment lines
  --c1 INTEGER                    Chrom 1 column; default 1
  --c2 INTEGER                    Chrom 2 column; default 3
  --p1 INTEGER                    Position 1 column; default 2
  --p2 INTEGER                    Position 2 column; default 4
  --s1 INTEGER                    Strand 1 column; default 5
  --s2 INTEGER                    Strand 2 column; default 6
  --unmapped-chrom TEXT           Placeholder for a chromosome on an unmapped
                                  side; default !
  --mark-dups                     If specified, duplicate pairs are marked as
                                  DD in "pair_type" and as a duplicate in the
                                  sam entries.
  --extra-col-pair TEXT...        Extra columns that also must match for two
                                  pairs to be marked as duplicates. Can be
                                  either provided as 0-based column indices or
                                  as column names (requires the "#columns"
                                  header field). The option can be provided
                                  multiple times if multiple column pairs must
                                  match. Example: --extra-col-pair "phase1"
                                  "phase2"
  --nproc-in INTEGER              Number of processes used by the auto-guessed
                                  input decompressing command.  [default: 3]
  --nproc-out INTEGER             Number of processes used by the auto-guessed
                                  output compressing command.  [default: 8]
  --cmd-in TEXT                   A command to decompress the input file. If
                                  provided, fully overrides the auto-guessed
                                  command. Does not work with stdin. Must read
                                  input from stdin and print output into
                                  stdout. EXAMPLE: pbgzip -dc -n 3
  --cmd-out TEXT                  A command to compress the output file. If
                                  provided, fully overrides the auto-guessed
                                  command. Does not work with stdout. Must
                                  read input from stdin and print output into
                                  stdout. EXAMPLE: pbgzip -c -n 8
  -h, --help                      Show this message and exit.

```




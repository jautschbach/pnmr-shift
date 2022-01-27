cp TAPE21.bak TAPE21
cp TAPE10.bak TAPE10

$ADFBIN/cpl << eor
gga
hybrid
hyperfine-scalar
  scf converge 1e-8 iterations 25
  atoms 1 2
end
end input
eor

rm TAPE21 TAPE10 logfile

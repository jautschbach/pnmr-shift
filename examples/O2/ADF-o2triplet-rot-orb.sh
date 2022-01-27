cp TAPE21.bak TAPE21
cp TAPE10.bak TAPE10

$ADFBIN/nmr << eor
allinone
nmr
 u1k best
 calc all
 out iso tens
 atoms 1 2
end
end input
eor

rm TAPE21 TAPE10 TAPE15 logfile

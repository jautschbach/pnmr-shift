cp TAPE21.bak TAPE21
cp TAPE10.bak TAPE10

$ADFBIN/nmr << eor
allinone
nmr
 gfactors
 u1k best
 calc all
 out iso tens
end
end input
eor

rm TAPE21 TAPE10 TAPE15 logfile

@echo off

echo "Creating pnmr.nod.out"
call ..\..\windows\PNMRShift.exe -f pnmrshift.nod.inp --splitinp --detail -c geo.xyz -t 298 -s 1.0 > pnmr.nod.out

echo "Creating pnmr.out"
call ..\..\windows\PNMRShift.exe -f pnmrshift.inp --splitinp --detail -c geo.xyz -t 298 -s 1.0 > pnmr.out

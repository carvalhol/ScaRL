#!/bin/bash

instruction='tail -n2 out_run.txt'
echo; echo ----------------------------
echo; echo ----------------------------
echo nPROCS     =9
echo nSELECT    =3
echo nPROCSTOTAL=27
echo nMEMORY    =22500
echo; echo ----------------------------
echo ------------ size = 500 ------------
echo \(cd TESTS/BIG-P027-SZ500\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ500; $instruction )
echo ------------ size = 490 ------------
echo \(cd TESTS/BIG-P027-SZ490\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ490; $instruction )
echo ------------ size = 480 ------------
echo \(cd TESTS/BIG-P027-SZ480\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ480; $instruction )
echo ------------ size = 470 ------------
echo \(cd TESTS/BIG-P027-SZ470\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ470; $instruction )
echo ------------ size = 460 ------------
echo \(cd TESTS/BIG-P027-SZ460\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ460; $instruction )
echo ------------ size = 450 ------------
echo \(cd TESTS/BIG-P027-SZ450\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ450; $instruction )
echo ------------ size = 440 ------------
echo \(cd TESTS/BIG-P027-SZ440\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ440; $instruction )
echo ------------ size = 430 ------------
echo \(cd TESTS/BIG-P027-SZ430\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ430; $instruction )
echo ------------ size = 420 ------------
echo \(cd TESTS/BIG-P027-SZ420\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ420; $instruction )
echo ------------ size = 410 ------------
echo \(cd TESTS/BIG-P027-SZ410\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ410; $instruction )
echo ------------ size = 400 ------------
echo \(cd TESTS/BIG-P027-SZ400\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ400; $instruction )
echo ------------ size = 390 ------------
echo \(cd TESTS/BIG-P027-SZ390\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ390; $instruction )
echo ------------ size = 380 ------------
echo \(cd TESTS/BIG-P027-SZ380\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ380; $instruction )
echo ------------ size = 370 ------------
echo \(cd TESTS/BIG-P027-SZ370\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ370; $instruction )
echo ------------ size = 360 ------------
echo \(cd TESTS/BIG-P027-SZ360\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ360; $instruction )
echo ------------ size = 350 ------------
echo \(cd TESTS/BIG-P027-SZ350\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ350; $instruction )
echo ------------ size = 340 ------------
echo \(cd TESTS/BIG-P027-SZ340\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ340; $instruction )
echo ------------ size = 330 ------------
echo \(cd TESTS/BIG-P027-SZ330\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ330; $instruction )
echo ------------ size = 320 ------------
echo \(cd TESTS/BIG-P027-SZ320\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ320; $instruction )
echo ------------ size = 310 ------------
echo \(cd TESTS/BIG-P027-SZ310\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ310; $instruction )
echo ------------ size = 300 ------------
echo \(cd TESTS/BIG-P027-SZ300\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ300; $instruction )
echo ------------ size = 290 ------------
echo \(cd TESTS/BIG-P027-SZ290\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ290; $instruction )
echo ------------ size = 280 ------------
echo \(cd TESTS/BIG-P027-SZ280\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ280; $instruction )
echo ------------ size = 270 ------------
echo \(cd TESTS/BIG-P027-SZ270\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ270; $instruction )
echo ------------ size = 260 ------------
echo \(cd TESTS/BIG-P027-SZ260\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ260; $instruction )
echo ------------ size = 250 ------------
echo \(cd TESTS/BIG-P027-SZ250\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ250; $instruction )
echo ------------ size = 240 ------------
echo \(cd TESTS/BIG-P027-SZ240\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ240; $instruction )
echo ------------ size = 230 ------------
echo \(cd TESTS/BIG-P027-SZ230\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ230; $instruction )
echo ------------ size = 220 ------------
echo \(cd TESTS/BIG-P027-SZ220\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ220; $instruction )
echo ------------ size = 210 ------------
echo \(cd TESTS/BIG-P027-SZ210\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ210; $instruction )
echo ------------ size = 200 ------------
echo \(cd TESTS/BIG-P027-SZ200\; qsub run.pbs\)
(cd TESTS/BIG-P027-SZ200; $instruction )
echo; echo ----------------------------
echo; echo ----------------------------
echo nPROCS     =8
echo nSELECT    =8
echo nPROCSTOTAL=64
echo nMEMORY    =20000
echo; echo ----------------------------
echo ------------ size = 500 ------------
echo \(cd TESTS/BIG-P064-SZ500\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ500; $instruction )
echo ------------ size = 490 ------------
echo \(cd TESTS/BIG-P064-SZ490\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ490; $instruction )
echo ------------ size = 480 ------------
echo \(cd TESTS/BIG-P064-SZ480\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ480; $instruction )
echo ------------ size = 470 ------------
echo \(cd TESTS/BIG-P064-SZ470\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ470; $instruction )
echo ------------ size = 460 ------------
echo \(cd TESTS/BIG-P064-SZ460\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ460; $instruction )
echo ------------ size = 450 ------------
echo \(cd TESTS/BIG-P064-SZ450\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ450; $instruction )
echo ------------ size = 440 ------------
echo \(cd TESTS/BIG-P064-SZ440\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ440; $instruction )
echo ------------ size = 430 ------------
echo \(cd TESTS/BIG-P064-SZ430\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ430; $instruction )
echo ------------ size = 420 ------------
echo \(cd TESTS/BIG-P064-SZ420\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ420; $instruction )
echo ------------ size = 410 ------------
echo \(cd TESTS/BIG-P064-SZ410\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ410; $instruction )
echo ------------ size = 400 ------------
echo \(cd TESTS/BIG-P064-SZ400\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ400; $instruction )
echo ------------ size = 390 ------------
echo \(cd TESTS/BIG-P064-SZ390\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ390; $instruction )
echo ------------ size = 380 ------------
echo \(cd TESTS/BIG-P064-SZ380\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ380; $instruction )
echo ------------ size = 370 ------------
echo \(cd TESTS/BIG-P064-SZ370\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ370; $instruction )
echo ------------ size = 360 ------------
echo \(cd TESTS/BIG-P064-SZ360\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ360; $instruction )
echo ------------ size = 350 ------------
echo \(cd TESTS/BIG-P064-SZ350\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ350; $instruction )
echo ------------ size = 340 ------------
echo \(cd TESTS/BIG-P064-SZ340\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ340; $instruction )
echo ------------ size = 330 ------------
echo \(cd TESTS/BIG-P064-SZ330\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ330; $instruction )
echo ------------ size = 320 ------------
echo \(cd TESTS/BIG-P064-SZ320\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ320; $instruction )
echo ------------ size = 310 ------------
echo \(cd TESTS/BIG-P064-SZ310\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ310; $instruction )
echo ------------ size = 300 ------------
echo \(cd TESTS/BIG-P064-SZ300\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ300; $instruction )
echo ------------ size = 290 ------------
echo \(cd TESTS/BIG-P064-SZ290\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ290; $instruction )
echo ------------ size = 280 ------------
echo \(cd TESTS/BIG-P064-SZ280\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ280; $instruction )
echo ------------ size = 270 ------------
echo \(cd TESTS/BIG-P064-SZ270\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ270; $instruction )
echo ------------ size = 260 ------------
echo \(cd TESTS/BIG-P064-SZ260\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ260; $instruction )
echo ------------ size = 250 ------------
echo \(cd TESTS/BIG-P064-SZ250\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ250; $instruction )
echo ------------ size = 240 ------------
echo \(cd TESTS/BIG-P064-SZ240\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ240; $instruction )
echo ------------ size = 230 ------------
echo \(cd TESTS/BIG-P064-SZ230\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ230; $instruction )
echo ------------ size = 220 ------------
echo \(cd TESTS/BIG-P064-SZ220\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ220; $instruction )
echo ------------ size = 210 ------------
echo \(cd TESTS/BIG-P064-SZ210\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ210; $instruction )
echo ------------ size = 200 ------------
echo \(cd TESTS/BIG-P064-SZ200\; qsub run.pbs\)
(cd TESTS/BIG-P064-SZ200; $instruction )

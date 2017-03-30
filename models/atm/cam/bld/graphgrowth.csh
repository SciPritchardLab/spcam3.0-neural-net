#!/bin/csh -f
#
#	graph_growth.csh			Jim Rosinski
#
#	Graphs the error growth. Input is the log
#	output from "cprnc" the CAM tool that compares
#	history files. Uses "xgraph", must be run on
#	a machine that has "xgraph" available.
#
#	$Id: graphgrowth.csh,v 1.4 2001/11/06 18:42:03 erik Exp $
#

#
# Ask user about a few options
#
set var = "T"
echo "enter variable for RMS [$var]"
set ans = $<
if ( $ans != "" ) set var = $ans
echo "Will plot RMS $var"

echo "Plot max as well as rms [n]?"
set ans = $<
if ( $ans != "y" ) then
  echo "Plotting rms only"
  unset maxtoo
else
  set maxtoo
  echo "Plotting rms and max"
endif

#
# Setup the basic commands for plotting
#
echo BoundBox: true      >! xg.dat
echo LogY: true          >> xg.dat
echo Ticks: true         >> xg.dat
echo XUnitText: NSTEP    >> xg.dat
echo YUnitText: $var >> xg.dat

#
# Enter each case that you want to plot in, loop until stop entering cases
# (Keep list of cprout files to enter for defaults)
#
set caselist = `ls *.cprout`;
set i = 1;
nextcase:

if ( $i < $#caselist )then
  set default =  ${caselist[$i]:r}
else
  set default =  "p"
endif
echo "Enter Case name: default $default (p to plot)"
set case = $<
@ i = $i + 1
if $case == "" set case = $default
if $case == "p" goto plotit

#
# Process each case stripping out data to plot the output logs of "cprnc"
#
if ( ! -f $case.cprout ) then
  echo "$case.cprout not found"
  goto nextcase
endif
grep NSTEPH: $case.cprout | cut -c 19- >! $case.nstep.grepout
grep "RMS $var " $case.cprout | cut -c 26- >! $case.$var.grepout
paste -d"\t\n" $case.nstep.grepout $case.$var.grepout >! $case.pasteout.rms
if ( $?maxtoo ) then
  grep "^ $var " $case.cprout | cut -c 71-77 >! $case.$var.grepout
  paste -d"\t\n" $case.nstep.grepout $case.$var.grepout >! $case.pasteout.max
endif
set nlinesnstep = `wc -l $case.nstep.grepout`
set nlinesvar   = `wc -l $case.$var.grepout`
\rm $case.nstep.grepout $case.$var.grepout
if ( $nlinesvar[1] < 1 ) then
  echo "No lines containing RMS $var were found in $case.cprout"
  goto nextcase
endif
if ( $nlinesnstep[1] != $nlinesvar[1] ) then
  echo "$nlinesnstep[1] NSTEP lines but $nlinesvar[1] RMS lines were found"
  echo "Quitting."
  goto cleanup
endif

echo \"$case RMS\"         >> xg.dat
cat $case.pasteout.rms     >> xg.dat
\rm $case.pasteout.rms
echo ""                    >> xg.dat
if ( $?maxtoo ) then
  echo \"$case MAX\"       >> xg.dat
  cat $case.pasteout.max   >> xg.dat
  \rm $case.pasteout.max
  echo ""                  >> xg.dat
endif
#
# If find a zero replace it with a small number, so that log plots will still work.
#
sed s/0.0000E+00/1.0000E-15/ xg.dat >! xg.dat.tmp
/bin/mv -f xg.dat.tmp xg.dat
goto nextcase

#
# Plotting, loop as long as ask to change limits
#
plotit:

unset userlim
echo -n "let xgraph pick ordinate limit? [y]"
set ans = $<
if ( $ans == n ) set userlim

if ( $?userlim ) then
  set ymin = -14
  echo -n "enter exponent on y min [1.e$ymin]"
  set ans = $<
  if ( $ans != "" ) set ymin = $ans

  set ymax = -6
  echo -n "enter exponent on y max [1.e$ymax]"
  set ans = $<
  if ( $ans != "" ) set ymax = $ans
  xgraph -ly $ymin,$ymax xg.dat
else
  xgraph xg.dat
endif

echo "Do you want to change ordinate limits? [n]"
set ans = $<
if ( $ans == y ) goto plotit

#
# End, cleanup files
#
cleanup:

\rm $case.pasteout.rms >& /dev/null
\rm $case.pasteout.max >& /dev/null

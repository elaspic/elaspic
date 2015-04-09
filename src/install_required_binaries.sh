#!/usr/bin/env bash

# This shell scripts downloads and compiles, and installs (by moving the executables 
# to the`elaspic/bin/` folder) most programs that are required by ELASPIC. 
 
# Note that Tcoffee, Modeller, Blast+, and FoldX, have to be installed manually.

ELASPIC_ROOT=`pwd`
if [[ "$ELASPIC_ROOT" == */src ]] ; then
  ELASPIC_ROOT="${ELASPIC_ROOT::-4}"
fi
PATH_TO_SRC="$ELASPIC_ROOT/src"
PATH_TO_BIN="$ELASPIC_ROOT/bin"


# DSSP 
cd "$PATH_TO_SRC"
printf '\n\n'
printf '*%.0s' {1..100}
printf "\nPreparing DSSP...\n"
if [ ! -e ./dssp*.tgz ] ; then 
  wget ftp://ftp.cmbi.ru.nl/pub/software/dssp/dssp-2.2.1.tgz
fi
tar xzvf ./dssp*.tgz
cd dssp*
make
cp -f "./mkdssp" "$PATH_TO_BIN/"


# MSMS
# http://mgltools.scripps.edu/
cd "$PATH_TO_SRC"
printf '\n\n'
printf '*%.0s' {1..100}
printf "\nPreparing MSMS...\n"
mkdir -p msms
if [ ! -e ./msms*.tar.gz ] ; then
  wget http://mgltools.scripps.edu/downloads/tars/releases/MSMSRELEASE/REL2.6.1/msms_i86_64Linux2_2.6.1.tar.gz
fi
mkdir -p ./msms
tar xzvf msms*.tar.gz -C ./msms
cd ./msms
# Change `nawk` to `gawk`, since `nawk` is not installed on some platforms
sed -i 's|nawk|gawk|g' pdb_to_xyzrn
sed -i 's|resnum=substr($0,23,4);|resnum=substr($0,23,4);\n\tchain=substr($0,21,2);|' pdb_to_xyzrn
sed -i 's|gsub(" ", "", aname);|gsub(" ", "", aname);\n\tgsub(" ", "", resnum);\n\tgsub(" ", "", chain);|' pdb_to_xyzrn
sed -i 's|"%f %f %f %f %d %s_%s_%d\\n", x, y, z, |"%f %f %f %f %d %s_%s_%d_%s\\n", x, y, z, |g' pdb_to_xyzrn
sed -i 's|1, aname, resname, resnum);|1, aname, resname, resnum, chain);|g' pdb_to_xyzrn
# Export files to ``elaspic/bin``
cp -f ./pdb_to_xyzr $PATH_TO_BIN
cp -f ./pdb_to_xyzrn $PATH_TO_BIN
cp -f ./atmtypenumbers $PATH_TO_BIN
cp -f ./msms.x86_64Linux2.2.6.1 $PATH_TO_BIN/msms


# CD-HIT
# https://code.google.com/p/cdhit/
cd $PATH_TO_SRC
printf '\n\n'
printf '*%.0s' {1..100}
printf "\nPreparing CD-HIT...\n"
if [ ! -e ./cdhit ] ; then
    hg clone https://code.google.com/p/cdhit/
fi
cd cdhit
make clean
make
cp -f cd-hit $PATH_TO_BIN


# PROVEAN
# http://provean.jcvi.org/downloads.php
cd $PATH_TO_SRC
printf '\n\n'
printf '*%.0s' {1..100}
printf "\nPreparing Provean...\n"
if [ ! -e ./provean.tar.gz ] ; then
    wget http://ufpr.dl.sourceforge.net/project/provean/provean-1.1.5.tar.gz
fi
tar xzvf provean*.tar.gz
cd provean*
./configure
make clean
make
cp -f ./src/provean $PATH_TO_BIN


# POPS
# http://mathbio.nimr.mrc.ac.uk/wiki/Software#POPS.2A
cd $PATH_TO_SRC
printf '\n\n'
printf '*%.0s' {1..100}
printf "\nPreparing POPS...\n"
if [ ! -e ./pops*.tar.gz ] ; then
    wget http://mathbio.nimr.mrc.ac.uk/download/POPS/pops-1.6.2.tar.gz
fi
tar xzvf pops*.tar.gz
cd pops*
./configure
make clean
make
cp -f ./src/pops $PATH_TO_BIN


# KNOT
# http://mathbio.nimr.mrc.ac.uk/wiki/Software#KNOT
cd $PATH_TO_SRC
printf '\n\n'
printf '*%.0s' {1..100}
printf "\nPreparing KNOT...\n"
if [ ! -e ./KNOT*.tgz ] ; then
    wget http://mathbio.nimr.mrc.ac.uk/download/KNOT/KNOT-1.0.0.tgz
fi
tar xzvf KNOT*.tgz
cd KNOT*
make
cp -f topol $PATH_TO_BIN


# Stride
# http://webclu.bio.wzw.tum.de/stride/
cd $PATH_TO_SRC
printf '\n\n'
printf '*%.0s' {1..100}
printf "\nPreparing Stride...\n"
if [ ! -e ./stride*.tar.gz ] ; then
    wget http://webclu.bio.wzw.tum.de/stride/stride.tar.gz
fi
mkdir -p ./stride
tar xzvf stride.tar.gz -C ./stride
cd stride
make
cp -f stride $PATH_TO_BIN


# Libfaketime
cd $PATH_TO_SRC
printf '\n\n'
printf '*%.0s' {1..100}
printf "\nPreparing Libfaketime...\n"
if [ ! -e ./libfaketime ] ; then
    git clone git@github.com:wolfcw/libfaketime.git
fi
cd libfaketime
make
make test
cp ./src/libfaketime.so.1 $PATH_TO_BIN

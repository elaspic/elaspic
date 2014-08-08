# Tcoffee and modeller have to be installed separately
# To install tcoffee download the *.bin file from http://www.tcoffee.org/Packages/Stable/Latest
# and follow instructions in http://www.tcoffee.org/Projects/tcoffee/#DOWNLOAD on how to install


PATH_TO_SRC=/home/kimlab1/strokach/working/elaspic/src
PATH_TO_BIN=/home/kimlab1/strokach/working/elaspic/bin_ubuntu

# FoldX, DSSP and MSMS come as pre-compiled binaries that should work on most
# x86_64 linux systems. Those binaries are included in the elaspic/src/bin
# folder, and we copy them to the working bin folder.
cp -f $PATH_TO_SRC/bin/* $PATH_TO_BIN


# MSMS
# http://mgltools.scripps.edu/downloads/tars/releases/MSMSRELEASE/REL2.6.1/msms_i86_64Linux2_2.6.1.tar.gz
# While MSMS binaries are provided in the elaspic/src/bin folder, the 
# pdb_to_xyzr file nay need to be edited depending on whether your system
# uses awk or gawk.
cd $PATH_TO_SRC
WHICH_GAWK=`which gawk`
WHICH_AWK=`which awk`
if [ ${WHICH_GAWK:${#WHICH_GAWK}-4} != "gawk" ]
then
    if [ ${WHICH_AWK:${#WHICH_AWK}-3} == "awk" ]
    then
        sed 's/gawk/awk/g' $PATH_TO_BIN/pdb_to_xyzrn > $PATH_TO_BIN/pdb_to_xyzrn
        cp -f ./bin/pdb_to_xyzrn.awk $PATH_TO_BIN/pdb_to_xyzrn
    else
        echo "No gawk or awk found! MSMS will not work!"
        exit 1
    fi
fi


# CD-HIT
# https://code.google.com/p/cdhit/
cd $PATH_TO_SRC
if ! [ -e ./cd-hit*.tgz ]
    then wget https://cdhit.googlecode.com/files/cd-hit-v4.6.1-2012-08-27.tgz
fi
tar xzvf cd-hit*.tgz
cd cd-hit*
make clean
make
cp -f cd-hit $PATH_TO_BIN


# PROVEAN
# http://provean.jcvi.org/downloads.php
cd $PATH_TO_SRC
if ! [ -e ./provean.tar.gz ]
    then wget -O provean.tar.gz http://sourceforge.net/projects/provean/files/latest/download?source=dlp
fi
tar xzvf provean.tar.gz
cd provean*
./configure
make clean
make
cp -f ./src/provean $PATH_TO_BIN


# POPS
# http://mathbio.nimr.mrc.ac.uk/wiki/Software#POPS.2A
cd $PATH_TO_SRC
if ! [ -e ./pops*.tar.gz ]
    then wget http://mathbio.nimr.mrc.ac.uk/download/POPS/pops-1.6.2.tar.gz
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
if ! [ -e ./KNOT*.tgz ]
    then wget http://mathbio.nimr.mrc.ac.uk/download/KNOT/KNOT-1.0.0.tgz
fi
tar xzvf KNOT*.tgz
cd KNOT*
make
cp -f topol $PATH_TO_BIN



# Stride
# http://webclu.bio.wzw.tum.de/stride/
if ! [ -e ./stride*.tar.gz ]
    then wget http://webclu.bio.wzw.tum.de/stride/stride.tar.gz
fi
tar xzvf stride.tar.gz
cd stride
make
cp -f stride $PATH_TO_BIN



### Turns out Tcoffee installs these automatically in its plugins folder
#~ # SAP
#~ # http://mathbio.nimr.mrc.ac.uk/wiki/Software#SAP
#~ cd $PATH_TO_SRC
#~ if ! [ -e ./sap*.tgz ]
    #~ then wget http://mathbio.nimr.mrc.ac.uk/download/SAP/sap-1.1.3.tgz
#~ fi
#~ tar xzvf sap*.tgz
#~ cd sap
#~ ./configure
#~ make clean
#~ make
#~ make check
#~ cp -f ./src/sap $PATH_TO_BIN
#~ 
#~ 
#~ # Mustang
#~ # http://www.csse.monash.edu.au/~karun/Site/mustang.html
#~ cd $PATH_TO_SRC
#~ if ! [ -e ./mustang*.tgz ]
    #~ then wget http://www.csse.monash.edu.au/~karun/mustang/mustang_v3.2.1.tgz
#~ fi
#~ tar xzvf mustang*.tgz
#~ cd MUSTANG*
#~ make clean
#~ make
#~ cp -f ./bin/mustang-3.2.1 $PATH_TO_BIN
#~ 
#~ 
#~ # TMalign
#~ # http://zhanglab.ccmb.med.umich.edu/TM-align/
#~ cd $PATH_TO_SRC
#~ mkdir -p TMalign
#~ cd TMalign
#~ if ! [ -e ./TMalign.f ]
    #~ then wget http://zhanglab.ccmb.med.umich.edu/TM-align/TMalign.f
#~ fi
#~ gfortran -static -O3 -ffast-math -lm -o TMalign TMalign.f
#~ # using gfortran is recommended, but you can also use g77:
#~ # g77 -static -O3 -lm -o TMalign TMalign.f 
#~ cp -f TMalign $PATH_TO_BIN

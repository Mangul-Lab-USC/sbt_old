#!/bin/bash

source $(dirname $0)/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-r', '--forse', action='store_true', default=False,
help='Forse [default %(default)s]')

EOF



DIR=`dirname $(readlink -f "$0")`

# Set default options.
CLEAN_ONLY=false
FORCE=false
NATIVE=false
REINSTALL=false


if [[ $FORSE ]]
then
REINSTALL=true
fi



# ------------------------------------------------------------------------------
# DOWNLOAD TOOLS
# ------------------------------------------------------------------------------

cd "$DIR/tools"

pwd
ls
echo $REINSTALL

# Skip this section if neither -c nor -r are selected and there is a previous
# installation (as indicated by the presence of the imrep directory).
echo '----- Checking for existing installations --------------------------------------'
if [ $CLEAN_ONLY = false ] && [ $REINSTALL = false ] && [ -d 'imrep' ]; then
    echo 'Existing installation found. Skipping tools download. To reinstall,' \
        'please use the -r option.'
else
    echo '----- Removing previous versions -----------------------------------------------'
    rm -fr imrep metaphlan2 MiniConda
    if [ $CLEAN_ONLY = true ]; then
        echo 'Done: Cleaning complete.'
        exit 0
    fi

    # Download ImReP.
    echo '----- Downloading ImRep --------------------------------------------------------'
    git clone https://github.com/mandricigor/imrep.git
    cd imrep
    ./install.sh
    cd ..

    #Download megahit
    echo '----- Downloading Megahit --------------------------------------------------'
    git clone https://github.com/voutcn/megahit.git
    cd megahit
    make
    cd ..

    # Download MetaPhlAn 2.
    echo '----- Downloading MetaPhlAn 2 --------------------------------------------------'
    hg clone https://bitbucket.org/biobakery/metaphlan2
    cd metaphlan2
    ln -s ../../db_human/databases
    cd ..

    # Download MiniConda and add shebangs.
    echo '----- Setting up Python environment --------------------------------------------'
    if [ $NATIVE = false ]; then
        ./install-MiniConda.sh
        cd MiniConda/lib
        ln -s libncursesw.so.5 libtinfow.so.5
        cd ../..
        MiniConda="$PWD/MiniConda/bin/python"
    #    sed -i "1c #!$MiniConda" metaphlan2/metaphlan2.py
    #    sed -i "1c #!$MiniConda" metaphlan2/strainphlan.py
    #    sed -i "1c #!$MiniConda" metaphlan2/utils/read_fastx.py
    #else
    #    sed -i '1c #!/usr/bin/env python2.7' metaphlan2/metaphlan2.py
    #    sed -i '1c #!/usr/bin/env python2.7' metaphlan2/strainphlan.py
    #    sed -i '1c #!/usr/bin/env python2.7' metaphlan2/utils/read_fastx.py
    fi
fi

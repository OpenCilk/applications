#!/bin/bash
#overwrite the next two variables for your setup
export top_dir="/home/angelee/cilk_dev/cilk-M/src_for_IntelCC"
export gcc_lib_dir="/home/angelee/cilk_dev/gcc_cilk/libcilkrts/.libs"
export icc_lib_dir="/home/angelee/intel/lib/intel64"
#export icc_lib_dir="/opt/intel/lib/intel64"

#overwrite this if the name of cilk-m runtime library changes
export cilkmrts_lib="libcilkm.so.0.0.0"
export cilkmrts_lib_in_rts=$top_dir"/runtime/.libs/"$cilkmrts_lib
#overwrite this if the name of cilk plus runtime library changes
export cilkplusrts_lib="libcilkrts.so.5.orig"
export gcccilkrts_lib="libcilkrts.so.5.0.2068"

export cilkmrts_sym_link="libcilkm.so.0"
export cilkplusrts_sym_link="libcilkrts.so.5"
export lib_input_file="libcilkrts.so"


if [ "$1" != "cilkplus" -a "$1" != "cilkm" -a "$1" != "gcc" ]; then 
    echo "Error: Unknown input value $1.  Accepted values: cilkplus, cilkm, gcc"
    exit 1;
fi

# setting up cilkplus
if [ "$1" == "cilkplus" -o "$1" == "gcc" ]; then
    if [ -d $icc_lib_dir ]; then 
        echo "cd into $icc_lib_dir"
        cd $icc_lib_dir
    else
        echo "Intel compiler lib directory does not exist, QUIT."
        exit 1;
    fi

	# make backup for the real cilkplusrts lib
	if [ ! -e $cilkplusrts_lib ]; then
		cp $cilkplusrts_sym_link $cilkplusrts_lib
		cp $cilkplusrts_sym_link $cilkplusrts_sym_link.bak
	fi

    echo "INPUT ($cilkplusrts_lib)" >& $lib_input_file
	if [ "$1" == "gcc" ]; then
	    echo "Setting up $cilkplusrts_sym_link with gcc cilk rts lib"
	    echo "cp -f $gcc_lib_dir/$gcccilkrts_lib $icc_lib_dir"
	    cp -f $gcc_lib_dir/$gcccilkrts_lib $icc_lib_dir
	    echo "ln -sf $gcccilkrts_lib $cilkplusrts_sym_link"
	    ln -sf $gcccilkrts_lib $cilkplusrts_sym_link
	else
	    echo "Setting up $cilkplusrts_sym_link with cilk plus lib"
	    echo "ln -sf $cilkplusrts_lib $cilkplusrts_sym_link"
	    ln -sf $cilkplusrts_lib $cilkplusrts_sym_link
	fi
    exit 0;
fi

# else setting up cilkm
if [ ! -e $cilkmrts_lib_in_rts ]; then
    make -j
fi

# check if the intel compiler lib dir actually exists
if [ -d $icc_lib_dir ]; then 
    echo "Copying Cilk-M runtime library to Intel compiler lib directory: $icc_lib_dir"
    cp -f $cilkmrts_lib_in_rts $icc_lib_dir
else
    echo "Intel compiler lib directory does not exist, QUIT."
    exit 1;
fi

echo "cd into $icc_lib_dir"
cd $icc_lib_dir

# check if a symbolic link for the cilk-m rts library exists
# if so, move it to backup file
if [ -e $cilkmrts_sym_link ]; then
    echo "sym link $cilkmrts_sym_link exists - remove."
    rm -f $cilkmrts_sym_link
fi

# create the symbolic link for the cilk-M rts library
echo "create a sym link $cilkmrts_sym_link to point to $cilkmrts_lib"
ln -sf $cilkmrts_lib $cilkmrts_sym_link

# now create the libcilkrts.so to point to the right sym link 
echo "Setting up $lib_input_file with cilkm lib"
echo "INPUT ($cilkmrts_sym_link)" >& $lib_input_file

#!/usr/bin/env tcsh

if ($#argv != 2) then
    echo "Usage: proc_hsqc.com file_root output.ft2"
    echo "Reads file_root_pos.fid and file_root_neg.fid."
    exit 1
endif

set file_root=$1
set out_file=$2
set work_dir="${file_root}_nmrpipe_work"

rm -rf -- "$work_dir"
mkdir -- "$work_dir"
if ($status != 0) exit 1

if ($?NMRBIN) then
    set nmrpipe="$NMRBIN/nmrPipe"
    set addnmr="$NMRBIN/addNMR"
else
    set nmrpipe=nmrPipe
    set addnmr=addNMR
endif

$nmrpipe -in "${file_root}_pos.fid" -fn FT                   \
          -out "$work_dir/pos_f2.ft" -ov
if ($status != 0) exit 1

$nmrpipe -in "${file_root}_neg.fid" -fn FT                   \
          -out "$work_dir/neg_f2.ft" -ov
if ($status != 0) exit 1

$nmrpipe -in "$work_dir/neg_f2.ft" -fn SIGN -i               \
          -out "$work_dir/neg_f2_conj.ft" -ov
if ($status != 0) exit 1

$addnmr -in1 "$work_dir/pos_f2.ft"                           \
        -in2 "$work_dir/neg_f2_conj.ft"                      \
        -out "$work_dir/states_f2.ft" -add
if ($status != 0) exit 1

$nmrpipe -in "$work_dir/states_f2.ft"                        \
| $nmrpipe -fn TP                                            \
| $nmrpipe -fn FT                                            \
| $nmrpipe -fn TP                                            \
          -out "$out_file" -ov
if ($status != 0) exit 1

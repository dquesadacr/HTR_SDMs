#!/bin/sh

dirname=$1

# for i in 1 2 3 4; do # different predictors in the sets, set to 1 (all variables) after testing
for i in 1; do
        base=`pwd`
        repo="HTR_SDMs"
        mkdir -p "$dirname"/C"$i" && cd "$_"

        git clone https://github.com/dquesadacr/"$repo".git
        mv "$repo"/{*,.*} ./
        mv ./0_Code/{*,.*} ./
        rmdir "$repo"/ 0_Code/

        cd $base

        echo sbatch -c 6 --time=6:00:00 --mem=11G -J colinvar_C"$i"_"$dirname" -o logs/colinvar_C"$i"_"$dirname"_%j.out -e logs/colinvar_C"$i"_"$dirname"_%j.err --mail-user dannell.quesada@pik-potsdam.de --mail-type END ./run_colinvar.sh $i $dirname "$2"

        sbatch -c 6 --time=6:00:00 --mem=11G -J colinvar_C"$i"_"$dirname" -o logs/colinvar_C"$i"_"$dirname"_%j.out -e logs/colinvar_C"$i"_"$dirname"_%j.err --mail-user dannell.quesada@pik-potsdam.de --mail-type END ./run_colinvar.sh $i $dirname "$2"
done

# case $1 in
#     1) # full -> all variables
#     filt=" ";;
#     2) # no radiation
#     filt="Rn_";;
#     3) # no indexes based on days
#     filt="days|growing|FFP";;
#     4) # no Rn and no days
#     filt="days|growing|FFP|Rn_";;
#     *)
#     echo Undefined case, aborting...
#     exit;;
# esac

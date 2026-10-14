#!/bin/bash
#------------------------------------------------------------------------------
# License
#     This file is part of moorFV, licensed under the GNU General Public
#     License <http://www.gnu.org/licenses/>.
#
# Script
#     bash_tip.sh
#
# Description
#     Write the x-component of W on beam patch left at each time to
#     left_W_x.dat.
#
#------------------------------------------------------------------------------

outFile="left_W_x.dat"
tmpFile=$(mktemp)

for Wfile in [0-9]*/beam_0/W; do
    [ -f "$Wfile" ] || continue

    timeDir=${Wfile%%/*}

    [[ "$timeDir" == "0" || "$timeDir" == "0.orig" ]] && continue

    Wx=$(awk '
        $1 == "left" {inLeft=1; next}

        inLeft && /value/ {
            line=$0
            sub(/.*\(/, "", line)
            sub(/\).*/, "", line)
            split(line, a, /[[:space:]]+/)
            print a[1]
            exit
        }

        inLeft && /^[[:space:]]*}/ {inLeft=0}
    ' "$Wfile")

    [ -n "$Wx" ] && echo "$timeDir $Wx" >> "$tmpFile"
done

{
    echo "Time W"
    sort -g -k1,1 "$tmpFile"
} > "$outFile"

rm "$tmpFile"

echo "Written to $outFile"

#!/bin/bash

outFile="right_W_x.dat"

echo "Time W" > "$outFile"

for Wfile in [0-9]*/beam_0/W; do
    timeDir=$(echo "$Wfile" | cut -d'/' -f1)

    Wx=$(awk '
        /right/ {inRight=1}
        inRight && /value/ {
            match($0, /\(([eE0-9+\.-]+)[[:space:]]+/, a)
            if (a[1] != "") {
                print a[1]
                exit
            }
        }
        inRight && /}/ {inRight=0}
    ' "$Wfile")

    if [ -n "$Wx" ]; then
        echo "$timeDir $Wx" >> "$outFile"
    fi
done

sort -n "$outFile" -o "$outFile"

echo "Written to $outFile"
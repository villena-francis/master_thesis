#!/bin/bash

extension=("aux" "log" "out" "toc" "bbl" "bcf" "blg" "synctex.gz" "fls" "fdb_latexmk" "xdv" "xml")

for ext in "${extension[@]}"
do
  find . -name "*.$ext" -type f -delete
done

echo "Temp LaTeX files deleted."

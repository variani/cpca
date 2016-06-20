#!/bin/bash

dirs=("demo")

dir_out="inst/doc"
name="output"
file_Rmd="$name.Rmd"

author="https://github.com/variani/cpca"
R_package="cpca"

opt_pandoc="-s --toc -mathjax -c normalize.css -c screen.css --self-contained" 

date_curr=$(date +"%d/%m/%Y")


header_str_1="% Output from R code
% $author
% $date_curr"

header_str_2='
```{r setup, include = FALSE}
opts_chunk$set(fig.path = "figure/", tidy = FALSE, dev = "svg")
# upload images automatically? (identity/imgur_upload)
opts_knit$set(upload.fun = identity)
```'
for dir in ${dirs[@]} ; do
  echo " dir: $dir"
  files=$(ls $dir | grep R) 

  echo " * file processing: R > Rmd"
  echo "$header_str_1" > $file_Rmd
  echo "$header_str_2" >> $file_Rmd

  pat='```'
  
  # Section Include
  echo "# Include" >> $file_Rmd
  echo "$pat{r}" >> $file_Rmd
  echo "library($R_package)" >> $file_Rmd
  echo "$pat" >> $file_Rmd

  echo "$pat{r}" >> $file_Rmd
  echo "packageVersion(\"$R_package\")" >> $file_Rmd
  echo "$pat" >> $file_Rmd

  # Section $dir
  echo "# $dir" >> $file_Rmd
  
  for file in $files ; do
    echo $file
  
    echo "## $file" >> $file_Rmd
    echo "$pat{r $file}" >> $file_Rmd
    cat $dir/$file >> $file_Rmd
    echo "$pat" >> $file_Rmd
  done

  echo "# R session info" >> $file_Rmd
  echo "$pat{r sessionInfo}" >> $file_Rmd
  echo "sessionInfo()" >> $file_Rmd
  echo "$pat" >> $file_Rmd

  echo " * knitr: Rmd -> md"
  R -q -e 'library(knitr);knit("output.Rmd")'

  echo " * pandoc: md -> html"
  cmd="pandoc $opt_pandoc $name.md -o $dir_out/$dir.html"
  #echo "cmd: $cmd"
  $cmd
  
  echo " * clean"
  rm -rf figure/ cache/ $name.md 
  rm -rf $name.Rmd
  
  echo " * ls"
  ls $dir_out/ | grep html
done


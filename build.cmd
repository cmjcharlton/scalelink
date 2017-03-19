R CMD build Scalelink
IF NOT EXIST ".\package\" mkdir package
R CMD INSTALL -l package --build Scalelink
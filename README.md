# jsonExforUtils - R package

This package is intended to contain functions 
useful to process and convert data in EXFOR fields.
At the moment it only contains functions to
parse the REACTION field.

## Installation

The package can be installed via
```
git clone https://github.com/gschnabel/jsonExforUtils.git
R CMD INSTALL jsonExforUtils
```

## Basic usage

The function `parseReacStr` extracts the information 
from a REACTION string and puts the result into a
list:

```
library(jsonExforUtils)

reacStr <- "(100-FM-257(N,ABS),,RI,,LIM,RECOM) ignored free text"
reacStruc <- parseReacStr(reacStr)
print(reacStruc)
```

Some reaction strings are arithmetic expressions of 
elementary reaction processes.
They can be processed with the function `parseReacExpr`,
e.g.,

```
reacStr <- "((64-GD-158(N,G)64-GD-159,,SIG,,MXW,RECOM)/ (64-GD-158(N,G)64-GD-159,,SIG,,,RECOM))"
reacStruc <- parseReacExpr(reacStr)
print(reacStruc)
```


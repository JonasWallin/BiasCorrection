# BiasCorrection

Data and code from article by

1. [David Bolin](http://www.math.chalmers.se/~bodavid/), 
2. [Arnoldo Frigessi](http://www.med.uio.no/imb/english/people/aca/frigessi/)
3. [Peter Guttorp](https://www.stat.washington.edu/peter/)
4. [Ola Haug](https://www.nr.no/~ola/)
5. Elisabeth Orksaug
6. [Ida Scheel](http://www.mn.uio.no/math/personer/vit/idasch/)
7. [Jonas Wallin](http://jonaswallin.github.io/)


### Code (Read before running code!)
The main code to run model 0,1,2 of the crossvalidation code is __Rcode/crossval.R__.
This code will download all the needed data, but not store it. If one would like to store the
data one needs to run the script __Rcode/DownloadData.R__, and then specify the location of the files in  __Rcode/crossval.R__ by changing the variable `data_location`.

The reference model can be runed using the script __Rcode/comparisonModel.R__

##### installing INLA
The R code requires use of INLA which is not aviable on CRAN, instead follows the instruction from
[rINLA](http://www.r-inla.org/download). Or in R, just run
`
install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable")
`


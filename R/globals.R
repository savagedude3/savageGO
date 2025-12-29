# Fix for R CMD Check "no visible binding for global variable"
if(getRversion() >= "2.15.1")  {
  utils::globalVariables(c(".", ".data", "desc"))
}

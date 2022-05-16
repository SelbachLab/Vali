

##########################################################################
#	 "Copyright (C) <2018>  <Henrik Zauber; MDC-Berlin in der Helmholtzgesellschaft>
#
#   This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>"
##########################################################################

#### Please run this script!!!!
# mac: cmd+e
# windows: drag and drop script in the console and press enter!find

args = commandArgs(trailingOnly=TRUE)

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
other.name <- file.path(script.basename, "other.R")
print(paste("Sourcing",other.name,"from",script.name))

#####
# checking for packages:
for(pck in c("shiny","RSQLite","pracma","parallel","gtools","h2o","Peptides","enviPat","shinyWidgets","compiler","tcltk","rstudioapi",
             "ggplot2","dplyr","hexbin","protViz","scales","tidyr","tidyverse","yaml"
             )){
  
if(!require(pck,character.only = T)){
  try({
    install.packages(pck)
  })
  missing <- !require(pck,character.only = T)
  if(missing){
    print("WARNING: Could not install",pck,". Please install manually.")
  }
}else{
  print(paste("found",pck))
}

}
# Installing RawDiag:
try({
  rdcheck <- require(rawDiag)
  if(rdcheck){
    print("Succesfully loaded rdcheck")
  }else{
    install.packages('http://fgcz-ms.uzh.ch/~cpanse/rawDiag_0.0.41.tar.gz', repo=NULL) 
    rdcheck <- require(rawDiag)
    if(rdcheck){
      print("Succesfully loaded rdcheck")
    }else{
      "Failed in installing RawDiag, please install manually."
    }
  }
})

#####
if(length(script.name)>0){
  print(args)
  SystemPath <- dirname(script.name)
}else{
  SystemPath <- dirname(rstudioapi::documentPath())
}
# Evaluating Libraries:



print(SystemPath)
print(getwd())
print(list.files())
ShinyStarter <- paste(SystemPath,"R/app_vali.R",sep = "/")
if(file.exists(ShinyStarter)){
  print("Starting Vali")
  print(ShinyStarter)
  runApp(ShinyStarter)
}





 

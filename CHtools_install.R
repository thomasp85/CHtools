pack <- c('xlsx', 'gWidgets', 'gWidgetsRGtk2', 'RGtk2', 'gridExtra', 'reshape2', 'plyr', 'RColorBrewer', 'scales', 'gtable')
install.packages(pack)
install.packages('ggplot2', repos='file:///R:/ASSAYS/@Pers_dir/TDP/R/')
install.packages('CHtools', repos='file:///R:/ASSAYS/@Pers_dir/TDP/R/')

library(RGtk2)

homedir <- path.expand('~')
rProfile <- file.path(homedir, '.Rprofile')
if(file.exists(rProfile)){
	cat('\n\nPlease include the line \"source(\'R:/ASSAYS/@Pers_dir/TDP/R_scripts/first.R\')\" in the .First function in your .Rprofile file.\n')
} else {
	file.rename('R:/ASSAYS/@Pers_dir/TDP/R_scripts/.Rprofile', homedir)
}
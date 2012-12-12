# TODO: Add comment
# 
# Author: Thomas
###############################################################################

### Import function for GC-MET
### Allows for merging of several xlsx files
importProcessReport <- function(files, startCol=6, GUI=FALSE){
	if(missing(files)){
		files <- choose.files()
	}
	data <- list()
	meta <- list()
	
	# Parse xlsx files
	for(i in 1:length(files)){
		tmp <- read.xlsx2(file=files[i], header=T, startRow=12, sheetIndex=1, check.names=F)
		tmp <- tmp[tmp[,1] != '',]
		unit <- read.xlsx2(file=files[i], startRow=11, endRow=11, sheetIndex=1, check.names=F, colClasses='character', header=F)
		unit <- unit[,!grepl("^c.+\\d*$", names(tmp), perl=T)]
		unit <- apply(unit, 2, as.character)
		tmp <- tmp[, !grepl("^c.+\\d*$", names(tmp), perl=T)]
		tmp[,startCol:ncol(tmp)] <- apply(tmp[,startCol:ncol(tmp)], 2, function(x) as.numeric(as.character(x)))
		names(tmp)[startCol:ncol(tmp)] <- paste(names(tmp)[startCol:ncol(tmp)], unit[startCol:ncol(tmp)])
		meta[[i]] <- tmp[,1:(startCol-1)]
		data[[i]] <- tmp[, startCol:ncol(tmp)]
	}
	
	# Test consistency between files
	testnames <- sapply(meta, names)
	if(!all(table(testnames) == length(meta))){
		if(GUI){
			answer <- gconfirm('Column names for sample info\ndoes not match between files', title='Warning', icon='warning', parent=window)
			if(!answer){
				return(NULL)
			} else {
				meta <- rbind.fill(meta)
				meta <- data.frame(lapply(meta, function(x) if(class(x) == 'factor') as.character(x) else x), stringsAsFactors=F, check.names=F)
				data <- rbind.fill(data)
				list(meta, data)
			}
		} else {
			cat('Column names for sample info does not match between files\n\n')
			answer <- readline('Continue? (y/n): ')
			answer <- tolower(answer)
			if(answer == 'n'){
				return(NULL)
			} else if(answer != 'y'){
				while(!answer %in% c('y', 'n')){
					answer <- readline('Wrong input. type \'y\' or \'n\': ')
				}
				if(answer == 'n'){
					return(NULL)
				} else {
					meta <- rbind.fill(meta)
					meta <- data.frame(lapply(meta, function(x) if(class(x) == 'factor') as.character(x) else x), stringsAsFactors=F, check.names=F)
					data <- rbind.fill(data)
					list(meta, data)
				}
			} else {
				meta <- rbind.fill(meta)
				meta <- data.frame(lapply(meta, function(x) if(class(x) == 'factor') as.character(x) else x), stringsAsFactors=F, check.names=F)
				data <- rbind.fill(data)
				list(meta, data)
			}
		}
	} else {
		meta <- rbind.fill(meta)
		meta <- data.frame(lapply(meta, function(x) if(class(x) == 'factor') as.character(x) else x), stringsAsFactors=F, check.names=F)
		data <- rbind.fill(data)
		list(meta, data)
	}
}

### Format a readline input
### Creates a numeric vector based on a mix of sequence and single number calls
stringToSeq <- function(x, sep= ' '){
	x <- unlist(strsplit(x, sep))
	x <- unlist(sapply(x, function(call) eval(parse(text=call))))
	names(x) <- NULL
	x
}

### Facetted plot of MCA data
### Asks for all the required information
plotProcessReport <- function(data, numeric=TRUE, errorbar=TRUE, pdf, plot='bar', ...){
	
	# Runs data import if missing
	if(missing(data)){
		data <- importProcessReport(...)
	}
	
	# Get information on important metadata columns
	metanames <- paste(1:ncol(data[[1]]), names(data[[1]]), sep=': ')
	for (i in 1:length(metanames)){
		cat(metanames[i], "\n")
	}
	cat('\n')
	category <- as.numeric(readline('Which column should the data be grouped by?: '))
	order <- as.numeric(readline('Which column should the data be ordered by?: '))
	
	# Removes selected samples
	samples <- unique(data[[1]][,category])
	cat('\nSamples in set:\n\n')
	for (i in 1:length(samples)){
		cat(i, ': ', samples[i], "\n")
	}
	remove <- readline('Samples to remove (space separeted - leave empty to keep all): ')
	if(remove != ''){
		remove <- stringToSeq(remove)
		remove <- which(data[[1]][,category] %in% samples[remove])
		data[[1]] <- data[[1]][-remove, ]
		data[[2]] <- data[[2]][-remove, ]
	}
	
	# Reformats the ordering column to numeric
	if(all(grepl('\\d+', data[[1]][, order])) && numeric){
		m <- regexpr('\\d+', data[[1]][, order], perl=T)
		data[[1]][, order] <- as.numeric(substr(data[[1]][, order], m, m+attr(m, 'match.length')-1))
	}
	
	# Selects selected compounds
	datanames <- paste(1:ncol(data[[2]]), names(data[[2]]), sep=': ')
	cat('\n')
	for (i in 1:length(datanames)){
		cat(datanames[i], "\n")
	}
	cat('\n')
	compounds <- readline('Compunds to plot (space separeted - leave empty to plot all): ')
	if(compounds != ''){
		compounds <- stringToSeq(compounds)
		data[[2]] <- data[[2]][, compounds]
		
		## Finds columns with only 0 and NA
		zeroCol <- which(apply(data[[2]], 2, max, na.rm=TRUE) < 0)
		if(length(zeroCol != 0)){
			cat('Also removed ', paste(names(data[[2]])[zeroCol], collapse=', '), 'with zero value data...\n')
			data[[2]] <- data[[2]][, -zeroCol]
		}
	} else {
		
		## Finds columns with only 0 and NA
		zeroCol <- which(apply(data[[2]], 2, max, na.rm=TRUE) < 0)
		if(length(zeroCol != 0)){
			cat('Removed ', paste(names(data[[2]])[zeroCol], collapse=', '), 'with zero value data...\n')
			data[[2]] <- data[[2]][, -zeroCol]
		}
	}
	
	# Formats plotting data
	plotdata <- cbind(data[[1]][, c(category, order)], data[[2]])
	plotdata <- melt(plotdata, id=1:2)
	names(plotdata) <- c('Cat', 'Ord', 'variable', 'value')
	plotdataM <- dcast(plotdata, Cat + Ord ~ variable, mean)
	plotdataM <- melt(plotdataM, id=1:2)
	plotdataS <- dcast(plotdata, Cat + Ord ~ variable, sd)
	plotdataS <- melt(plotdataS, id=1:2)
	plotdataM$ymin <- plotdataM$value - plotdataS$value
	plotdataM$ymax <- plotdataM$value + plotdataS$value
	allComb <- expand.grid(Cat=unique(plotdataM$Cat), Ord=unique(plotdataM$Ord), variable=unique(plotdataM$variable))
	plotdataM <- merge(allComb, plotdataM, all=TRUE)
	plotdataM$Ord <- factor(plotdataM$Ord)
	
	# Splits data if plot can't fit on one page
	subs <- split(unique(plotdataM$variable), rep(1:ceiling(length(unique(plotdataM$variable))/20), each=20)[1:length(unique(plotdataM$variable))])
	plots <- list()
	for(i in 1:length(subs)){
		
		# Plotting
		if(plot == 'bar'){
			plotdataTEMP <- subset(plotdataM, plotdataM$variable %in% subs[[i]])
			p <- ggplot(data=plotdataTEMP, aes(x=Cat, y=value, fill=Cat, colour=Ord))
			p <- p + geom_bar(position=position_dodge(0.9), stat='identity', size=0)
			if(errorbar){
				p <- p + geom_errorbar(aes(ymin=ymin, ymax=ymax), position=position_dodge(0.9), width=0)
			}
			if(length(subs) == 1 || missing(pdf)){
				p <- p + facet_wrap(~variable, scales='free') + theme_bw()
			} else {
				p <- p + facet_wrap(~variable, scales='free', ncol=4, nrow=5) + theme_bw()
			}
			p <- p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())
			p <- p + theme(strip.text=element_text(size=6)) + labs(x='', y='')
			
			# Selects colouring scheme based on number of samples
			if(length(unique(plotdataM$Cat)) < 10){
				p <- p + scale_fill_brewer(names(data[[1]])[category], type='qual', palette='Set1')
			} else if(length(unique(plotdataM$Cat)) < 13){
				p <- p + scale_fill_brewer(names(data[[1]])[category], type='qual', palette='Paired')
			} else {
				p <- p + scale_fill_hue(names(data[[1]])[category])
			}
			p <- p + scale_colour_manual(breaks=unique(plotdataM$Ord), values=rep('black', length(unique(plotdataM$Ord))), guide=guide_legend(title='Bar order', keywidth=0, keyheight=0, direction='horizontal', label.position='bottom', title.position='top', label.theme=element_text(angle=90, size=8), label.hjust=0.5, label.vjust=0.5, override.aes=list(alpha=0)))
			plots[[i]] <- p
		} else if(plot == 'line'){
			if(!numeric){
				stop('X-axis must be numeric for lineplots...')
			}
			plotdataTEMP <- subset(plotdataM, plotdataM$variable %in% subs[[i]])
			plotdataTEMP$Ord <- as.numeric(as.character(plotdataTEMP$Ord))
			p <- ggplot(data=plotdataTEMP, aes(x=Ord, y=value, colour=Cat, group=Cat))
			p <- p + geom_line()
			if(errorbar){
				p <- p + geom_errorbar(aes(ymin=ymin, ymax=ymax), width=0.1)
			}
			if(length(subs) == 1 || missing(pdf)){
				p <- p + facet_wrap(~variable, scales='free_y') + theme_bw()
			} else {
				p <- p + facet_wrap(~variable, scales='free_y', ncol=4, nrow=5) + theme_bw()
			}
			p <- p + theme(strip.text=element_text(size=6)) + labs(x=names(data[[1]])[order], y='')
			p <- p + theme(axis.text.x=element_text(size=6, angle=45, vjust=1, hjust=1))
			
			# Selects colouring scheme based on number of samples
			if(length(unique(plotdataM$Cat)) < 10){
				p <- p + scale_colour_brewer(names(data[[1]])[category], type='qual', palette='Set1')
			} else if(length(unique(plotdataM$Cat)) < 13){
				p <- p + scale_colour_brewer(names(data[[1]])[category], type='qual', palette='Paired')
			} else {
				p <- p + scale_colour_hue(names(data[[1]])[category])
			}
			plots[[i]] <- p
		} else {
			stop('Unknown plot type...')
		}
		
	}
	
	# Creating output
	if(missing(pdf)){
		for(i in 1:length(plots)){
			print(plots[[i]])
			if(i != length(plots)){
				readline('Press return for next page...')
			} else {
				cat('Done...\n')
			}
		}
	} else {
		pdf(file=pdf, width=9, height=12)
		for(i in 1:length(plots)){
			print(plots[[i]])
		}
		dev.off()
		cat('PDF file written to ', file.path(getwd(), pdf), '\n', sep='')
	}
	invisible(plotdataM)
}


formatProcessReport <- function(data, numeric, category, order, sample, compound){
	data[[1]] <- data[[1]][sample, , drop=FALSE]
	data[[2]] <- data[[2]][sample, , drop=FALSE]
	data[[2]] <- data[[2]][, compound, drop=FALSE]
	if(all(grepl('\\d+', data[[1]][, order])) && numeric){
		m <- regexpr('\\d+', data[[1]][, order], perl=T)
		data[[1]][, order] <- as.numeric(substr(data[[1]][, order], m, m+attr(m, 'match.length')-1))
	}
	zeroCol <- which(apply(data[[2]], 2, function(x) if(all(is.na(x)) || !all(is.numeric(x))) TRUE else max(x, na.rm=TRUE) < 0))
	if(length(zeroCol != 0)){
		data[[2]] <- data[[2]][, -zeroCol, drop=FALSE]
	}
	
	plotdata <- cbind(data[[1]][, c(category, order)], data[[2]])
	plotdata <- melt(plotdata, id=1:2)
	names(plotdata) <- c('Cat', 'Ord', 'variable', 'value')
	plotdataM <- dcast(plotdata, Cat + Ord ~ variable, mean)
	plotdataM <- melt(plotdataM, id=1:2)
	plotdataS <- dcast(plotdata, Cat + Ord ~ variable, sd)
	plotdataS <- melt(plotdataS, id=1:2)
	plotdataM$ymin <- plotdataM$value - plotdataS$value
	plotdataM$ymax <- plotdataM$value + plotdataS$value
	allComb <- expand.grid(Cat=unique(plotdataM$Cat), Ord=unique(plotdataM$Ord), variable=unique(plotdataM$variable))
	plotdataM <- merge(allComb, plotdataM, all=TRUE)
	plotdataM$Ord <- factor(plotdataM$Ord)
	
	plotdataM
}


#### GUI

MCAplot <- function(){
	
	options(guiToolkit='RGtk2')
	
	fplotProcessReport <- function(data, category, order, plottype, errorbar=TRUE){
		if(plottype == 'Bar'){
			p <- ggplot(data=data, aes(x=Cat, y=value, fill=Cat, colour=Ord))
			p <- p + geom_bar(position=position_dodge(0.9), stat='identity', size=0)
			if(errorbar){
				p <- p + geom_errorbar(aes(ymin=ymin, ymax=ymax), position=position_dodge(0.9), width=0)
			}
			if(length(unique(data$variable)) < 17){
				p <- p + facet_wrap(~variable, scales='free') + theme_bw()
			} else {
				p <- p + facet_wrap(~variable, scales='free', ncol=4, nrow=5) + theme_bw()
			}
			p <- p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())
			p <- p + theme(strip.text=element_text(size=6)) + labs(x='', y='')
			
			# Selects colouring scheme based on number of samples
			if(length(levels(plotdata$Cat)) < 10){
				p <- p + scale_fill_brewer(category, drop=FALSE, type='qual', palette='Set1')
			} else if(length(levels(plotdata$Cat)) < 13){
				p <- p + scale_fill_brewer(category, drop=FALSE, type='qual', palette='Paired')
			} else {
				p <- p + scale_fill_hue(category, drop=FALSE)
			}
			p <- p + scale_colour_manual(breaks=levels(plotdata$Ord), values=rep('black', length(unique(plotdata$Ord))), guide=guide_legend(title='Bar order', keywidth=0, keyheight=0, direction='horizontal', label.position='bottom', title.position='top', label.theme=element_text(angle=90, size=8), label.hjust=0.5, label.vjust=0.5, override.aes=list(alpha=0)))
			p
		} else if(plottype == 'Line'){
			data$Ord <- as.numeric(as.character(data$Ord))
			p <- ggplot(data=data, aes(x=Ord, y=value, colour=Cat, group=Cat))
			p <- p + geom_line()
			if(errorbar){
				p <- p + geom_errorbar(aes(ymin=ymin, ymax=ymax), width=0.1)
			}
			if(length(unique(data$variable)) < 17){
				p <- p + facet_wrap(~variable, scales='free_y') + theme_bw()
			} else {
				p <- p + facet_wrap(~variable, scales='free_y', ncol=4, nrow=5) + theme_bw()
			}
			p <- p + theme(strip.text=element_text(size=6)) + labs(x=order, y='')
			p <- p + theme(axis.text.x=element_text(size=6, angle=45, vjust=1, hjust=1))
			
			# Selects colouring scheme based on number of samples
			if(length(levels(plotdata$Cat)) < 10){
				p <- p + scale_colour_brewer(category, drop=FALSE, type='qual', palette='Set1')
			} else if(length(levels(plotdata$Cat)) < 13){
				p <- p + scale_colour_brewer(category, drop=FALSE, type='qual', palette='Paired')
			} else {
				p <- p + scale_colour_hue(category, drop=FALSE)
			}
			p
		}
	}
	importProcessReport <- function(files, startCol=6){
		data <- list()
		meta <- list()
		
		# Parse xlsx files
		for(i in 1:length(files)){
			tmp <- read.xlsx2(file=files[i], header=T, startRow=12, sheetIndex=1, check.names=F)
			tmp <- tmp[tmp[,1] != '',]
			unit <- read.xlsx2(file=files[i], startRow=11, endRow=11, sheetIndex=1, check.names=F, colClasses='character', header=F)
			unit <- unit[,!grepl("^c.+\\d*$", names(tmp), perl=T)]
			unit <- apply(unit, 2, as.character)
			tmp <- tmp[, !grepl("^c.+\\d*$", names(tmp), perl=T)]
			tmp[,startCol:ncol(tmp)] <- apply(tmp[,startCol:ncol(tmp)], 2, function(x) as.numeric(as.character(x)))
			names(tmp)[startCol:ncol(tmp)] <- paste(names(tmp)[startCol:ncol(tmp)], unit[startCol:ncol(tmp)])
			meta[[i]] <- tmp[,1:(startCol-1)]
			data[[i]] <- tmp[, startCol:ncol(tmp)]
		}
		
		# Test consistency between files
		testnames <- sapply(meta, names)
		if(!all(table(testnames) == length(meta))){
			answer <- gconfirm('Column names for sample info\ndoes not match between files', title='Warning', icon='warning', parent=window)
			if(!answer){
				return(NULL)
			} else {
				meta <- rbind.fill(meta)
				meta <- data.frame(lapply(meta, function(x) if(class(x) == 'factor') as.character(x) else x), stringsAsFactors=F, check.names=F)
				data <- rbind.fill(data)
				list(meta, data)
			}
		} else {
			meta <- rbind.fill(meta)
			meta <- data.frame(lapply(meta, function(x) if(class(x) == 'factor') as.character(x) else x), stringsAsFactors=F, check.names=F)
			data <- rbind.fill(data)
			list(meta, data)
		}
	}
	
	filelist <- character()
	category <- ''
	order <- ''
	samples <- character()
	compounds <- character()
	data <- list()
	plotdata <- data.frame()
	datasplit <- list()
	previewplot <- 1
	
	
	window <- gwindow('MCA plot', visible=FALSE)
	mainGroup <- ggroup(horizontal=TRUE, cont=window)
	controlGroup <- ggroup(horizontal=FALSE, cont=mainGroup)
	
	importBox <- gframe('Report files to use:', horizontal=FALSE, cont=controlGroup)
	addSpace(importBox, 3)
	fileSelect <- ggroup(horizontal=TRUE, cont=importBox)
	addSpace(fileSelect, 3)
	files <- gtable(filelist, cont=fileSelect)
	names(files) <- 'Files'
	size(files) <- c(300, 100)
	fbuttons <- ggroup(horizontal=FALSE, cont=fileSelect)
	addFile <- gbutton('Add', cont=fbuttons)
	addSpace(fbuttons, 6, horizontal=FALSE)
	removeFile <- gbutton('Remove', cont=fbuttons)
	enabled(removeFile) <- FALSE
	addSpace(fileSelect, 3)
	importBox2 <- ggroup(cont=importBox)
	addSpace(importBox2, 3)
	glabel('Data starts at column: ', cont=importBox2)
	startCol <- gspinbutton(from=1, to=50, by=1, digits=0, value=6, cont=importBox2)
	size(startCol) <- c(50,20)
	import <- gbutton('Import file(s)', cont=importBox2, expand=F)
	enabled(import) <- FALSE
	addSpace(importBox, 6)
	
	addSpace(controlGroup, 12)
	
	optionsBox <- ggroup(cont=controlGroup)
	addSpring(optionsBox)
	options <- glayout(cont=optionsBox, spacing=5, fill='y', expand=TRUE)
	optList <- list()
	options[1, 1, anchor=c(1,0)] <- 'Category:'
	options[1, 2, anchor=c(-1,0), expand=TRUE] <- optList$category <- gcombobox(category, cont=options)
	options[2, 1, anchor=c(1,0)] <- 'Order:'
	options[2, 2, anchor=c(-1,0), expand=TRUE] <- optList$order <- gcombobox(order, cont=options)
	options[2, 3, anchor=c(-1,0)] <- optList$numeric <- gcheckbox('numeric', checked=TRUE)
	options[3, 1, anchor=c(1,0)] <- 'Plottype:'
	options[3, 2, anchor=c(0,0)] <- optList$plottype <- gradio(c('Bar', 'Line'), selected=1, horizontal=TRUE, cont=options)
	options[3, 3, anchor=c(-1,0)] <- optList$errorbar <- gcheckbox('errorbars', checked=TRUE)
	addSpring(optionsBox)
	
	addSpace(controlGroup, 12)
	
	filterFrame <- gframe('Filtering', pos=0.5, cont=controlGroup, expand=TRUE)
	addSpace(filterFrame, 3)
	filterSamplesBox <- ggroup(cont=filterFrame, horizontal=FALSE, expand=T)
	addSpace(filterSamplesBox, 6)
	labelgroup1 <- ggroup(cont=filterSamplesBox)
	glabel('Samples:', cont=labelgroup1)
	addSpring(labelgroup1)
	all1 <- gcheckbox('All', cont=labelgroup1)
	filterSamples <- gcheckboxgroup(samples, checked=T, expand=T, use.table=T, cont=filterSamplesBox)
	addSpace(filterSamplesBox, 3)
	addSpace(filterFrame, 3)
	gseparator(horizontal=F, cont=filterFrame)
	addSpace(filterFrame, 3)
	filterCompoundBox <- ggroup(cont=filterFrame, horizontal=FALSE, expand=T)
	addSpace(filterCompoundBox, 6)
	labelgroup2 <- ggroup(cont=filterCompoundBox)
	glabel('Compounds:', cont=labelgroup2)
	addSpring(labelgroup2)
	all2 <- gcheckbox('All', cont=labelgroup2)
	filterCompound <- gcheckboxgroup(compounds, checked=T, expand=T, use.table=T, cont=filterCompoundBox)
	addSpace(filterCompoundBox, 3)
	addSpace(filterFrame, 3)
	
	
	plotGroup <- ggroup(horizontal=FALSE, cont=mainGroup)
	
	plotDevice <- ggraphics(cont=plotGroup)
	size(plotDevice) <- c(500*sqrt(2), 500)
	visible(plotDevice) <- TRUE
	
	plotButtons <- ggroup(cont=plotGroup, spacing=0)
	preview <- gbutton('Preview', cont=plotButtons)
	enabled(preview) <- FALSE
	addSpace(plotButtons, 6)
	back <- gbutton('<', cont=plotButtons)
	enabled(back) <- FALSE
	forward <- gbutton('>', cont=plotButtons)
	enabled(forward) <- FALSE
	addSpace(plotButtons, 6)
	count <- glabel('', cont=plotButtons)
	addSpring(plotButtons)
	save <- gbutton('Save', cont=plotButtons)
	enabled(save) <- FALSE
	
	
	addHandlerClicked(addFile, handler=function(h, ...){
				file <- gfile(type='open', multi=TRUE, filter=list('Excel (.xlsx)'=list(patterns=c('*.xlsx')), 'All files'=list(patterns=c('*'))))
				if(!is.na(file)){
					files[] <- as.character(c(files[], file))
					setwd(dirname(file))
				}
				if(length(files[]) != 0){
					enabled(removeFile) <- TRUE
					enabled(import) <- TRUE
				} else {
					enabled(removeFile) <- FALSE
					enabled(import) <- FALSE
				}
			}
	)
	addHandlerClicked(removeFile, handler=function(h, ...){
				remove <- svalue(files, index=TRUE)
				if(length(remove) != 0){
					files[] <- files[-remove,]
				}
				if(length(files[]) != 0){
					enabled(removeFile) <- TRUE
					enabled(import) <- TRUE
				} else {
					enabled(removeFile) <- FALSE
					enabled(import) <- FALSE
				}
			}
	)
	addHandlerClicked(import, handler=function(h, ...){
				data <<- list()
				filterSamples[] <- character()
				data <<- importProcessReport(files[], startCol=svalue(startCol))
				optList$category[] <- c('', names(data[[1]]))
				svalue(optList$category) <-  ''
				optList$order[] <- c('', names(data[[1]]))
				svalue(optList$order) <-  ''
				filterCompound[] <- names(data[[2]])
				svalue(filterCompound) <- TRUE
				svalue(all2) <- TRUE
				visible(plotDevice) <- TRUE
				frame()
				enabled(preview) <- FALSE
				enabled(back) <- FALSE
				enabled(forward) <- FALSE
				enabled(save) <- FALSE
				svalue(count) <- ''
			}
	)
	addHandlerChanged(optList$category, handler=function(h, ...){
				if(length(svalue(optList$category)) != 0){
					if(svalue(optList$category) == ''){
						filterSamples[] <- character()
					} else {
						filterSamples[] <- unique(data[[1]][,svalue(h$obj)])
						svalue(filterSamples) <- TRUE
						svalue(all1) <- TRUE
					}
					if(svalue(optList$category) == '' || svalue(optList$order) == ''){
						enabled(preview) <- FALSE
						enabled(save) <- FALSE
					} else {
						enabled(preview) <- TRUE
						enabled(save) <- TRUE
					}
				}
			}
	)
	addHandlerChanged(optList$order, handler=function(h, ...){
				if(length(svalue(optList$order)) != 0){
					if(svalue(optList$order) != ''){
						if(!all(grepl('\\d+', data[[1]][, svalue(h$obj)])) && svalue(optList$numeric)){
							gmessage(paste('Missing numeric information in ', svalue(h$obj)), 'Warning', parent=window)
						}
					}
					if(svalue(optList$category) == '' || svalue(optList$order) == ''){
						enabled(preview) <- FALSE
						enabled(save) <- FALSE
					} else {
						enabled(preview) <- TRUE
						enabled(save) <- TRUE
					}
				}
			}
	)
	addHandlerChanged(optList$numeric, handler=function(h, ...){
				if(svalue(h$obj) && !all(grepl('\\d+', data[[1]][, svalue(optList$order)]))){
					gmessage(paste('Missing numeric information in ', svalue(optList$order)), 'Warning', parent=window)
				}
				if(!svalue(h$obj) && !svalue(optList$plottype) == 'Bar'){
					enabled(preview) <- FALSE
					enabled(save) <- FALSE
				} else {
					enabled(preview) <- TRUE
					enabled(save) <- TRUE
				}
			}
	)
	addHandlerChanged(optList$plottype, handler=function(h, ...){
				if(!svalue(h$obj) == 'Bar' && !svalue(optList$numeric)){
					enabled(preview) <- FALSE
					enabled(save) <- FALSE
				} else {
					enabled(preview) <- TRUE
					enabled(save) <- TRUE
				}
			}
	)
	addHandlerChanged(all1, handler=function(h, ...){
				svalue(filterSamples) <- svalue(h$obj)
			}
	)
	addHandlerChanged(all2, handler=function(h, ...){
				svalue(filterCompound) <- svalue(h$obj)
			}
	)
	addHandlerClicked(preview, handler=function(h, ...){
				previewplot <<- 1
				sample <- which(data[[1]][, svalue(optList$category)] %in% svalue(filterSamples))
				plotdata <<- formatProcessReport(data, numeric=svalue(optList$numeric), category=svalue(optList$category), order=svalue(optList$order), sample=sample, compound=svalue(filterCompound, index=TRUE))
				datasplit <<- split(unique(plotdata$variable), rep(1:ceiling(length(unique(plotdata$variable))/20), each=20)[1:length(unique(plotdata$variable))])
				p <- fplotProcessReport(subset(plotdata, plotdata$variable %in% datasplit[[previewplot]]), svalue(optList$category), svalue(optList$order), svalue(optList$plottype), errorbar=svalue(optList$errorbar))
				visible(plotDevice) <- TRUE
				print(p)
				visible(plotDevice) <- FALSE
				svalue(count) <- paste('(', previewplot, '/', length(datasplit), ')', sep='')
				if(previewplot != 1){
					enabled(back) <- TRUE
				}
				if(previewplot != length(datasplit)){
					enabled(forward) <- TRUE
				}
			}
	)
	addHandlerClicked(back, handler=function(h, ...){
				previewplot <<- previewplot - 1
				p <- fplotProcessReport(subset(plotdata, plotdata$variable %in% datasplit[[previewplot]]), svalue(optList$category), svalue(optList$order), svalue(optList$plottype), errorbar=svalue(optList$errorbar))
				visible(plotDevice) <- TRUE
				print(p)
				visible(plotDevice) <- FALSE
				svalue(count) <- paste('(', previewplot, '/', length(datasplit), ')', sep='')
				if(previewplot != 1){
					enabled(back) <- TRUE
				} else {
					enabled(back) <- FALSE
				}
				if(previewplot != length(datasplit)){
					enabled(forward) <- TRUE
				} else {
					enabled(forward) <- FALSE
				}
			}
	)
	addHandlerClicked(forward, handler=function(h, ...){
				previewplot <<- previewplot + 1
				p <- fplotProcessReport(subset(plotdata, plotdata$variable %in% datasplit[[previewplot]]), svalue(optList$category), svalue(optList$order), svalue(optList$plottype), errorbar=svalue(optList$errorbar))
				visible(plotDevice) <- TRUE
				print(p)
				visible(plotDevice) <- FALSE
				svalue(count) <- paste('(', previewplot, '/', length(datasplit), ')', sep='')
				if(previewplot != 1){
					enabled(back) <- TRUE
				} else {
					enabled(back) <- FALSE
				}
				if(previewplot != length(datasplit)){
					enabled(forward) <- TRUE
				} else {
					enabled(forward) <- FALSE
				}
			}
	)
	addHandlerClicked(save, handler=function(h, ...){
				file <- gfile(type='save', filter=list('PDF'=list(patterns=c('*.pdf')), 'All files'=list(patterns=c('*'))))
				sample <- which(data[[1]][, svalue(optList$category)] %in% svalue(filterSamples))
				pdata <- formatProcessReport(data, numeric=svalue(optList$numeric), category=svalue(optList$category), order=svalue(optList$order), sample=sample, compound=svalue(filterCompound, index=TRUE))
				dsplit <- split(unique(pdata$variable), rep(1:ceiling(length(unique(pdata$variable))/20), each=20)[1:length(unique(pdata$variable))])
				if(!grepl('*.pdf', file, ignore.case=TRUE)){
					file <- paste(file, '.pdf', sep='')
				}
				pdf(file, height=9, width=12)
				for(i in 1:length(dsplit)){
					p <- fplotProcessReport(subset(pdata, pdata$variable %in% dsplit[[i]]), svalue(optList$category), svalue(optList$order), svalue(optList$plottype), errorbar=svalue(optList$errorbar))
					print(p)
				}
				dev.off()
			}
	)
	visible(window) <- TRUE
	visible(plotDevice) <- TRUE
	frame()
}
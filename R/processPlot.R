# TODO: Add comment
# 
# Author: Thomas
###############################################################################

### Import function for GC-MET
### Allows for merging of several xlsx files
importProcessReport <- function(files, startCol=6){
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
		cat('Column names for sample info does not match between files\n\n')
		answer <- readline('Continue? (y/n): ')
		answer <- tolower(answer)
		if(answer == 'n'){
			stop('Import cancelled')
		} else if(answer != 'y'){
			while(!answer %in% c('y', 'n')){
				answer <- readline('Wrong input. type \'y\' or \'n\': ')
			}
			if(answer == 'n'){
				stop('Import cancelled')
			}
		}
	}
	
	# Format output
	meta <- rbind.fill(meta)
	meta <- data.frame(lapply(meta, function(x) if(class(x) == 'factor') as.character(x) else x), stringsAsFactors=F, check.names=F)
	data <- rbind.fill(data)
	list(meta, data)
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
plotProcessReport <- function(data, numeric=TRUE, errorbar=TRUE, pdf, ...){
	
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
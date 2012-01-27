
### CIS: keywords of titles with
### "break point" OR "breakpoint" OR "change point" OR "changepoint" OR 
### "structural change" OR "structural break"
files <- list.files(pattern = "\\.txt")

### convert to DCF format
citations <- unlist(sapply(files, function(file) readLines(file)))
citations <- gsub("^(%[A-Z][a-z]*)", "\\1:", citations, perl = TRUE)
citations <- gsub("^%", "", citations)
citations <- gsub(" = ", "", citations)
tmpfile <- tempfile()
writeLines(citations, con = tmpfile)
citations <- read.dcf(tmpfile)

### remove duplicates
dup <- duplicated(paste(citations[,"Author"], citations[,"Title"], citations[,"Year"], sep = ""))
citations <- citations[!dup,]

### select Journals with more than 15 citations
ncitations <- table(citations[,"Journal"])
refdata <- as.data.frame(citations[citations[,"Journal"] %in% names(ncitations[ncitations > 15]),])

### tabulate data and aggregate over journals
yearnref <- xtabs(~ Journal + Year, data = refdata)
nrefs <- colSums(yearnref)
year <- ordered(colnames(yearnref))

### test
library("coin")
maxstat_test(nrefs ~ year)


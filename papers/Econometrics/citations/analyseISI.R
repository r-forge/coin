
### CIS: keywords of titles with
### "break point" OR "breakpoint" OR "change point" OR "changepoint" OR 
### "structural change" OR "structural break"
files <- list.files(pattern = "save")

### convert to DCF format
citations <- unlist(sapply(files, function(file) readLines(file)))
citations <- gsub("^([A-Z][A-Z,0-9])", "\\1:", citations, perl = TRUE)
tmpfile <- tempfile()
writeLines(citations, con = tmpfile)
citations <- read.dcf(tmpfile)

refdata <- as.data.frame(citations)
save(refdata, file = "ISI_strucchange.rda")

### tabulate data and aggregate over journals
yearnref <- xtabs(~ SO + PY, data = refdata)
nrefs <- colSums(yearnref)
year <- ordered(colnames(yearnref))

### test
library("coin")
maxstat_test(nrefs ~ year)


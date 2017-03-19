data(FOI)
data(LDFCOMP)

idcols <- c("Day", "Month", "Year", "Sex")
result <- calcScores(FOI[, idcols], LDFCOMP[, idcols])

print(result$scores)
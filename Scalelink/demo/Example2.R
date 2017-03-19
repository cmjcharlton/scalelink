data(FOI)
data(LDF)

idcols <- c("Day", "Month", "Year", "Sex")
result <- calcScores(FOI[, idcols], LDF[, idcols])

print(result$scores)
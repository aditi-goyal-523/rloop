#get 5 number summary of the processed information

df <- read.delim("C:/Users/aditi/Desktop/korf_lab/rloop/processed_narrowpeak/LS61A.neg.processed.txt")

attach(df)
summary(df$length)
sd(df$length)
boxplot(df$length)

Q1 <- quantile(length, .25)
Q3 <- quantile(length, .75)
IQR <- IQR(length)

detach(df)

no_length_outliers <- subset(df, length > (Q1 - 1.5*IQR) & length < (Q3 + 1.5*IQR))

attach(no_length_outliers)
summary(length)
sd(length)
boxplot(length)

#-------------------------------------

summary(count)
sd(count)
boxplot(count)

Q1 <- quantile(count, .25)
Q3 <- quantile(count, .75)
IQR <- IQR(count)

no_count_outliers <- subset(no_length_outliers, count > (Q1 - 1.5*IQR) & count < (Q3 + 1.5*IQR))

detach(no_length_outliers)
attach(no_count_outliers)

summary(count)
sd(count)
boxplot(count)

View(no_count_outliers)

#--------------------------------------

hist(count)
hist(length)

summary(length)
sd(length)

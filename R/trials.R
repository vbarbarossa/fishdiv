df <- data.frame(
  y = sort(round(rbeta(100,1,10)*100,0)+1),
  x = sort(rbeta(100,1,1))+rnorm(100,0.3,0.15)
)
plot(df$x,df$y)
plot(log10(df$x),log10(df$y))

summary(lm(df$y ~ df$x))
summary(lm(log10(df$y) ~ log10(df$x)))

m <- glm(y~x,data = df)
summary(m)
plot(df$y,predict(m))
abline(0,1)

summary(glm(log10(y)~x,data = df))

mp <- glm(y~x,data = df,family = 'poisson')
summary(mp)
plot(df$y,predict(mp))
abline(0,1)

ml <- glm(log(y)~x,data = df)
summary(ml)
plot(df$y,predict(ml))
abline(0,1)



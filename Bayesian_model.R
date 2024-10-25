
library(nimble)
library(readr)
par_reference <- read_csv("~/par_reference.csv")

numleft = 1
numright = 300
par_reference <- par_reference[numleft:numright,]


genelist <- par_reference$...1
for (gene in genelist){
    #gene = 'LINC01128'
    print(gene)
    LINC00128 <- read_csv(paste0("~/",gene,".csv",sep=''))
    par_df = par_reference[par_reference$...1 == gene,]

    pumpCode <- nimbleCode({ 
      for (i in 1:N){
        gamma[i] <- rgamma(1, 1 + 1 / beta)
        y[i] <- 1 / (2 * scale * gamma[i]) * exp(-(abs(x[i] - loc) / scale) ^ beta)
      }

      loc ~ dunif(locmin,locmax) 
      scale ~ dunif(0,scalemax) 
      beta ~ dunif(betamin,betamax)
    })


    pumpConsts <- list(N = length(LINC00128$x),
                       x = LINC00128$x,
                       locmin = as.double(par_df[1:1,'loc'])-0.1,
                       locmax = as.double(par_df[1:1,'loc'])+0.1,
                       scalemax = as.double(par_df[1:1,'scale'])+0.01,
                       betamin = as.double(par_df[1:1,'beta'])-0.1,
                       betamax = as.double(par_df[1:1,'beta'])+0.1)

    pumpData <- list(y = LINC00128$y)

    pumpInits <- list(loc = 1,
                      scale = 0.1,
                      beta = 0.3)


    pumpModel <- nimbleModel(code = pumpCode, 
                             name = 'pump', 
                             constants = pumpConsts,
                             data = pumpData, 
                             inits = pumpInits)
    pumpConf <- configureMCMC(pumpModel,
                              monitors = c('loc','scale','beta'),
                              thin = 100)
    pumpMCMC <- buildMCMC(pumpConf)
    CpumpModel <- compileNimble(pumpModel)
    CpumpMCMC <- compileNimble(pumpMCMC, project = pumpModel)

    fit <- runMCMC(mcmc = CpumpMCMC, niter = 100000, nburnin =1000)
    print(summary(fit))
    plot(fit[ , 'loc'], type = 'l', xlab = 'iteration',  ylab = expression(loc))
    plot(fit[ , 'scale'], type = 'l', xlab = 'iteration',  ylab = expression(scale))
    plot(fit[ , 'beta'], type = 'l', xlab = 'iteration',  ylab = expression(beta))
    fit <- as.data.frame(fit)
    par_reference[par_reference$...1==gene,'fit.beta'] <- mean(fit$beta)
    par_reference[par_reference$...1==gene,'fit.loc'] <- mean(fit$loc)
    par_reference[par_reference$...1==gene,'fit.scale'] <- mean(fit$scale)

}

write.table(par_reference,paste0("~/par_reference_MCMC_",numleft,"_",numright,".csv",sep=''),
          row.names = FALSE, col.names = TRUE, sep = ',')

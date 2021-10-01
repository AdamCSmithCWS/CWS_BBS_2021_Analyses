model

{

	#### counts and overdispersion effects  ###### Hierarchical GAM model with additional random year-effects

	####  builds on the GAM model used in 2016 work

	### not yet applied or tested





	for( k in 1 : ncounts )

	{

		log(lambda[k]) <- obs[strat[k],obser[k]] + eta*firstyr[k] + strata[strat[k]] + yeareffect[strat[k],year[k]] + yy[strat[k],year[k]] + noise[k]

	 	#noise[k] ~ dnorm(0, taunoise)

		noise[k] ~ dt(0, taunoise, 3) #alternative t-distributed noise = heavy-tailed overdispersion
		
		count[k] ~ dpois(lambda[k])




	}






	### fixed effect priors



   sdnoise ~ dt(0, 1, 3) T(0,)

   taunoise <- 1 / pow(sdnoise, 2)


	eta ~ dnorm( 0.0,1)

	STRATA ~ dnorm( 0.0,1)



	sdstrata ~ dgamma(2,2)  #<- 1/pow(sdbeta,2)#

	taustrata <- 1/pow(sdstrata,2)#~ dunif(0.001,10)



  sdobs ~ dnorm(0,4) T(0,)#<- 

  tauobs <- 1/pow(sdobs, 2)#dgamma(0.001,0.0001)



	#### stratum-level effects  ######

	for( i in 1 : nstrata )

	{

		#### observer effects  ######



		for( o in 1 : nobservers[i] )

		{

			obs[i,o] ~ dnorm(0.0, tauobs)

		}

		#### end observer effects  ######





		strata[i] ~ dnorm(STRATA,taustrata) #<-  + strata.p[i]


	}# end s strata loop and stratum-level effects









	###########COMPUTING GAMs for yeareffects##############

  sdBETA ~ dt(0, 1, 3) T(0,) #dgamma(2,0.5) #alternate prior, original from Cainiceanu et al. second gamma parameter == 0.0001 << (abs(mean(BETA[]))^2)/2, mean(BETA[]) ~ 0.2
  
  tauBETA <- 1/pow(sdBETA,2) # prior on precision of gam hyperparameters
  
  
  

	


	for(k in 1:nknots)

	{
	  BETA[k] ~ dnorm(0,tauBETA)
	  
	  sdbeta[k] ~ dt(0, 1, 3) T(0,) #dgamma(2,0.5) # prior on precision of gam coefficients(
	  
	  taubeta[k] <- 1/(pow(sdbeta[k],2))
	  
		# Computation of GAM components



	  
	  for(i in 1:nstrata){
	    beta[i,k] ~ dnorm(BETA[k],taubeta[k])# T(-20,20) #avoiding log density calculation errors
	  }
	  
		for(i in 1:nstrata)

		{


			for ( t in ymin : ymax )

			{

				X.part[i,k,t] <- beta[i,k]*(X.basis[t,k])

			}#t

		}#i

	}#k



		for(i in 1:nstrata)

		{

		for (t in ymin : ymax )

	{

			yeareffect[i,t] <- sum(X.part[i,1:nknots,t])

		}#t

	}#i





	#-------------------------------------------------#





	#### additional random year effects  ######


		for( i in 1 : nstrata )

    {

		  	for( t in (ymin) : ymax )

		{

			yy[i,t] ~ dnorm(0,tauyy[i])
      
		  	  	}
		  
		  
		  sdyy[i] ~ dgamma(2,2)
		  tauyy[i] <- 1/pow(sdyy[i],2)
		  




	}











	#### summary statistics  ######

	sdn_ret <- 0.5*sdnoise*sdnoise

	sdobs_ret <- 0.5*sdobs*sdobs



	for( i in 1 : nstrata )

	{
	  
	  sdyy_ret[i] <- 0.5*sdyy[i]*sdyy[i]
	  
		for( t in ymin : ymax )

		{

			for(o in 1 : nobservers[i])

			{

				no[i,t,o] <- exp(strata[i]+yeareffect[i,t] + yy[i,t] + obs[i,o] + sdn_ret)

				nosmooth[i,t,o] <- exp(strata[i]+yeareffect[i,t] + obs[i,o] + sdn_ret + sdyy_ret[i])

			}



			mn[i,t] <- mean(no[i,t,(1 : nobservers[i])])

			mnsmooth[i,t] <- mean(nosmooth[i,t,(1 : nobservers[i])])

			n[i,t] <- nonzeroweight[i]*(mn[i,t])

			nsmooth[i,t] <- nonzeroweight[i]*(mnsmooth[i,t])



			n2[i,t] <- nonzeroweight[i]*exp(strata[i]+yeareffect[i,t] + yy[i,t] + sdn_ret + sdobs_ret) #n2 is an alternative way of calculating n

			n4[i,t] <- nonzeroweight[i]*exp(strata[i]+yeareffect[i,t] + sdn_ret + sdobs_ret) #n4 is an alternative way of calculating n3

		}

	}





	#-------------------------------------------------#

}



#==========================================================================      
# Given genetical paramters, calculate power for three model selection strategies
#==========================================================================      
require(mvtnorm)
require(adapt)
require(corpcor) # for func is.positive.definite().
                         
##*** Main function ***##

markerSearchPower = function(mainEff1, mainEff2, epistasisEff, noiseSD, alleleFreq1, alleleFreq2, alleleFreq3=0.5, 
                                DetectionN, obsN, markerN, 
                                strategy = c("marginal", "exhaustive", "forward"),
                                powerDef = c("both", "either"), samplePointN=20000,
                                decompMethod = c("svd", "chol", "eigen")) {
        b1 = mainEff1; b2 = mainEff2; b3 = epistasisEff; sigm = noiseSD; 
        p1 = alleleFreq1; p2 = alleleFreq2; p3 = alleleFreq3;
        R = DetectionN; n = obsN; p = markerN; 
        
        ### calculate parameters for relevant distributions
        meanT1 = meanT1Func(b1, b2, b3, p1, p2, sigm, n);
        meanT2 = meanT1Func(b2, b1, b3, p2, p1, sigm, n);
        varT1 = varT1Func(b1, b2, b3, p1, p2, sigm);
        varT2 = varT1Func(b2, b1, b3, p2, p1, sigm);
        covT1T2 = covT1T2Func(b1, b2, b3, p1, p2, sigm);
        covMT1T2 = matrix(c(varT1, covT1T2, covT1T2, varT2), ncol=2);
        
        meanT12 = meanT12Func(b1, b2, b3, p1, p2, sigm, n) ;  
        varT12  = varT12Func(b1, b2, b3, p1, p2, sigm); 
        covT12T1= covT12T1Func(b1, b2, b3, p1, p2, sigm);
        covT12T2= covT12T1Func(b2, b1, b3, p2, p1, sigm); 
        covMT12s= matrix(c(varT12, covT12T1, covT12T2, covT12T1, varT1, 
                            covT1T2, covT12T2, covT1T2, varT2), ncol=3);
        
        meanF3l1 = meanF3l1Func(b1, b2, b3, p1, p2, p3, sigm) ; 
        varF3l1 = varF3l1Func(b1, b2, b3, p1, p2, p3, sigm);
        c3l1=varF3l1/(2*meanF3l1); d3l1=2*meanF3l1^2/varF3l1;

        meanF3l2 = meanF3l1Func(b2, b1, b3, p2, p1, p3, sigm) ; 
        varF3l2 = varF3l1Func(b2, b1, b3, p2, p1, p3, sigm);
        c3l2=varF3l2/(2*meanF3l2); d3l2=2*meanF3l2^2/varF3l2;

        #distn parameters for Tvector=c(T1, T2, T1l3, T2l3)
        meanT1l3 = meanT1l3Func(b1, b2, b3, p1, p2, sigm, n);
        meanT2l3 = meanT1l3Func(b2, b1, b3, p2, p1, sigm, n);
        varT1l3 = varT1l3Func(b1, b2, b3, p1, p2, sigm) ;
        varT2l3 = varT1l3Func(b2, b1, b3, p2, p1, sigm);
        covT1T1l3 = covT1T1l3Func(b1, b2, b3, p1, p2, sigm);
        covT2T2l3 = covT1T1l3Func(b2, b1, b3, p2, p1, sigm);
        covT1T2l3 = covT1T2l3Func(b1, b2, b3, p1, p2, sigm) ;
        covT2T1l3 = covT1T2l3Func(b2, b1, b3, p2, p1, sigm);
        covT1l3T2l3 = covT1l3T2l3Func(b1, b2, b3, p1, p2, sigm);
        covMTvector = matrix( c(varT1, covT1T2, covT1T1l3, covT1T2l3,
                                covT1T2, varT2, covT2T1l3, covT2T2l3,
                                covT1T1l3, covT2T1l3, varT1l3, covT1l3T2l3,
                                covT1T2l3, covT2T2l3, covT1l3T2l3, varT2l3), ncol=4);
        
        ### calculate power for model selection strategies.
        if (strategy=="marginal" && powerDef=="both") {
            R = DetectionN -1;
            return( power.marginal(meanT1, meanT2, covMT1T2, R, p, TRUE, samplePointN) )
        }
        if (strategy=="marginal" && powerDef=="either") 
            return( power.marginal(meanT1, meanT2, covMT1T2, R, p, FALSE, samplePointN) )
    
        if (strategy=="exhaustive" && powerDef=="both") 
            return( power.exhaustive.both(meanT12, meanT1, meanT2, covMT12s, c3l1, d3l1, c3l2, d3l2, R, p, n, samplePointN))
        if (strategy=="exhaustive" && powerDef=="either") 
            return( power.exhaustive.atleast1(meanT12, meanT1, meanT2, covMT12s, c3l1, d3l1, c3l2, d3l2, R, p, n, samplePointN))
    
        if (strategy=="forward" && powerDef=="both") 
            return( power.forward.both(meanT12, meanT1, meanT2, covMT12s, c3l1, d3l1, c3l2, d3l2, R, p, n, samplePointN) )
        if (strategy=="forward" && powerDef=="either") 
            return( power.forward.atleast1(meanT1, meanT2, meanT1l3,meanT2l3, covMTvector, R, p, samplePointN, decompMethod[1]) )
            
        cat(paste("Power calculation for finding ", powerDef, " true markers with ", strategy, " search strategy failed", sep=""));
}




##*** Internal Functions ***##


######################################################################################
###   relevant distribution functions
######################################################################################

    CDF.T3  = function(x) {pnorm(x,0,1)}               #cdf of any Tj, j=3,...,p
    CDF.F3l1= function(x, c, d) {pchisq((1/c)*x, d)}   #cdf of any Fjli, i=1,2, j=3,...,p
    CDF.F34 = function(x) {pchisq(3*x, 3)}              #cdf of any Fjk, j<k, and j,k=3,...,p
    PDF.F34 = function(x) {3*dchisq(x*3, 3)}            #pdf of any Fjk, j<k, and j,k=3,...,p
    CDF.F4l3= function(x) {pchisq(2*x, 2)}              #cdf of any Fklj, j<k, and j,k=3,...,p



######################################################################################
###   power calculation functions
######################################################################################

    ### Return marginal selection power for finding both true SNPs or at least one true SNP
    #Input: mean and covariance variance of T1 and T2, 
    #       number of false positive models R, number of markers p, boolean isBoth indicate power type.
    power.marginal = function(meanT1, meanT2, covMT1T2, R, p, isBoth, samplePointN) {
        condProb.both = function(x) { 
            N=p-2; r=N+1-R; #R[i]; #control signif level by R
            j=0:(N-r);
            summ=sum( choose(r+j-1, j)*(2*CDF.T3(-min(abs(x))))^j );
            return( (1-2*CDF.T3(-min(abs(x))))^r * summ );
        }
        condProb.atleast1 = function(x) { 
            N=p-2; r=N+1-R; #R[i];
            j=0:(N-r);
            summ=sum( choose(r+j-1, j)*(2*CDF.T3(-max(abs(x))))^j );
            return( (1-2*CDF.T3(-max(abs(x))))^r * summ );
        }

        x <- rmvnorm(samplePointN, mean=c(meanT1,meanT2), sigma=covMT1T2)#, method="chol")
        if (isBoth) result=apply(x, 1, condProb.both)
        if (!isBoth) result=apply(x, 1, condProb.atleast1)
        return(mean(result))
    }


    ### Return exhaustive selection power for finding both true SNPs 
    #Input: mean and covariance variance of T12, T1 and T2, 
    #       parameter for distn of F3l1 and F3l2: c3l1, d3l1, c3l2, d3l2
    #       number of false positive models R, number of markers p.
    power.exhaustive.both = function(meanT12, meanT1, meanT2, covMT12s, c3l1, d3l1, c3l2, d3l2, R, p, n, samplePointN) {
        condProb = function(x) { 
            y1 = (3*x[1]^2-x[2]^2)/(2*(1+meanT1^2/n));
            y2 = (3*x[1]^2-x[3]^2)/(2*(1+meanT2^2/n));
            z  = x[1]^2; F34N = choose(p-2,2);
            sumTriple = array(0, R); #R[Ri]); #sum of probs over all triples (r1,r2,r3) given r=r1+r2+r3
            sumTriple[1] = (CDF.F3l1(y1, c3l1, d3l1)*CDF.F3l1(y2, c3l2, d3l2))^(p-2) * CDF.F34(z)^F34N;
            if (R>1) { #[Ri]>1) { #to use apply(rs, 1, function{...}), rs has to have at least two rows.
                for (r in 1:(R-1)) { #[Ri]-1)) {
                    rs = as.array(tripleSet(r, p, F34N));
                    sumTriple[r+1] = sum( apply(rs, 1, function(rsi) {
                                                    (choose(p-2, rsi[1])*(1-CDF.F3l1(y1, c3l1, d3l1))^rsi[1]*CDF.F3l1(y1, c3l1, d3l1)^(p-2-rsi[1])
                                                    *choose(p-2, rsi[2])*(1-CDF.F3l1(y2, c3l2, d3l2))^rsi[2]*CDF.F3l1(y2, c3l2, d3l2)^(p-2-rsi[2])
                                                    *choose(F34N, rsi[3])*(1-CDF.F34(z))^rsi[3]*CDF.F34(z)^(F34N-rsi[3])) }) )
                }
            }
            return( sum(sumTriple) );
        }
        x <- rmvnorm(samplePointN, mean=c(meanT12, meanT1, meanT2), sigma=covMT12s)
        result=apply(x, 1, condProb)
        return(mean(result))
    }

            
    ### Return exhaustive selection power for finding at least one true SNPs 
    #Input: mean and covariance variance of T12, T1 and T2, 
    #       parameter for distn of F3l1 and F3l2: c3l1, d3l1, c3l2, d3l2
    #       number of false positive models R, number of markers p.
    power.exhaustive.atleast1 = function(meanT12, meanT1, meanT2, covMT12s, c3l1, d3l1, c3l2, d3l2, R, p, n, samplePointN) {
        cutoff = qchisq(1-(R-0.5)/choose(p-2,2), 3)/3
        condProb = function(x) { 
            u1 = (3*cutoff-x[2]^2)/(2*(1+meanT1^2/n));
            u2 = (3*cutoff-x[3]^2)/(2*(1+meanT2^2/n));
            return( (CDF.F3l1(u1, c3l1, d3l1)*CDF.F3l1(u2, c3l2, d3l2))^(p-2) * (x[1]^2<cutoff) );      
        }
        x <- rmvnorm(n=samplePointN, mean=c(meanT12, meanT1, meanT2), sigma=covMT12s)
        result=apply(x, 1, condProb)
        1- mean(result)
    }



    ### Forward selection
    ### Return forward selection power for finding both true SNPs 
    #Input: mean and covariance variance of T12, T1 and T2, 
    #       parameter for distn of F3l1 and F3l2: c3l1, d3l1, c3l2, d3l2
    #       number of false positive models R, number of markers p.
    power.forward.both = function(meanT12, meanT1, meanT2, covMT12s, c3l1, d3l1, c3l2, d3l2, R, p, n, samplePointN) {
        condProb = function(x) { #approx using quantiles
            maxT12 = max(abs(x[2:3]))
            if (abs(x[2]) >= abs(x[3])) {
                cutoff = qchisq(1-(R-1+0.5)/(p-2), d3l1)*c3l1 #distn of F3l1 is c3l1*chisq(d3l1).
                u=(3*x[1]^2-maxT12^2)/(2*(1+meanT1^2/n));
                return( (1-2*CDF.T3(-maxT12))^(p-2) * (u>cutoff) )
            }
            if (abs(x[2]) < abs(x[3])) {
                cutoff = qchisq(1-(R-1+0.5)/(p-2), d3l2)*c3l2 #distn of F3l2 is c3l2*chisq(d3l2).
                u=(3*x[1]^2-maxT12^2)/(2*(1+meanT2^2/n));
                return( (1-2*CDF.T3(-maxT12))^(p-2) * (u>cutoff) )
            }
        }
        #condProb = function(x) { #use cdf of order statistics
        #    maxT12 = max(abs(x[2:3]))
        #    
        #    N=p-2; r=N+1-R;
        #    j=0:(N-r);
        #    if (abs(x[2]) >= abs(x[3])) { #distn of F3l1 is c3l1*chisq(d3l1).
        #       u=(3*x[1]^2-maxT12^2)/(2*(1+meanT1^2/n));
        #       summ=sum( choose(r+j-1, j)*(1-CDF.F3l1(u, c3l1, d3l1))^j );
        #       return( CDF.F3l1(u, c3l1, d3l1)^r * summ * (1-2*CDF.T3(-maxT12))^(p-2));
        #    }
        #    if (abs(x[2]) < abs(x[3])) {  #distn of F3l2 is c3l2*chisq(d3l2).
        #       u=(3*x[1]^2-maxT12^2)/(2*(1+meanT2^2/n));
        #       summ=sum( choose(r+j-1, j)*(1-CDF.F3l1(u, c3l2, d3l2))^j );
        #       return( CDF.F3l1(u, c3l2, d3l2)^r * summ * (1-2*CDF.T3(-maxT12))^(p-2));
        #    }
        #}
        x <- rmvnorm(samplePointN, mean=c(meanT12,meanT1, meanT2), sigma=covMT12s) 
        result=apply(x, 1, condProb)
        mean(result) 
    }
        
        
 

    ### Return forward selection power for finding at least one true SNPs 
    #Input: mean and covariance of Tvector=c(T1, T2, T1l3, T2l3), 
    #       number of false positive models R, number of markers p.
    power.forward.atleast1 = function(meanT1, meanT2, meanT1l3,meanT2l3, covMTvector, R, p, samplePointN, decompMethod) {    
        condProb = function(x) { 
            return( (1-2*CDF.T3(-max(abs(x))))^(p-2) );
        }
        x <- rmvnorm(n=samplePointN, mean=c(meanT1,meanT2), sigma=covMTvector[1:2, 1:2])
        result=apply(x, 1, condProb)
        pow1 = mean(result)

        cutoff = qchisq(1-(R-1+0.5)/(p-3), 2)/2 #distn of F4l3 is chisq(2)/2.
        condProb = function(x) { 
            (1-(1-2*CDF.T3(-max(abs(x[1:2]))))^(p-2)) * (max(x[3]^2, x[4]^2)>cutoff)
        }
        #condProb = function(x) { 
        #    u=max(x[3]^2, x[4]^2);
        #    N=p-3; # # of null stats: F4l3,...,Fpl3
        #    r=N+1-R;
        #    j=0:(N-r);
        #    summ=sum( choose(r+j-1, j)*(1-CDF.F4l3(u))^j );
        #    return( CDF.F4l3(u)^r * summ * (1-(1-2*CDF.T3(-max(abs(x[1:2]))))^(p-2)) );
        #}
        if (sum(is.nan(covMTvector))>0) { #condi means F1l3 =(d)= F4l3
            pow2 = 0; # prob of pick up sign is cR/p -> 0 as R samll & p -> Inf
        }
        if (sum(is.nan(covMTvector))==0) {
            x <- rmvnorm(n=samplePointN, mean=c(meanT1, meanT2, meanT1l3,meanT2l3), sigma=covMTvector, method=decompMethod)  
            result=apply(x, 1, condProb)
            pow2 = mean(result)
        }
        return(pow1+ pow2)
    }

 
 
 
 
 
 
 
 ########################################################
 # Functions need in calculation.
 ########################################################
  
    # Function return all possible triples (r1, r2, r3) 
    #        s.t. r1+r2+r3=r; 0<=r1<=p-2, 0<=r2<=p-2, 0<=r3<=F34N.
    tripleSet = function(r, p, F34N) {
        result = array(0,3)
        for (r1 in 0:min(p-2, r)) {
            for (r2 in 0:min(p-2, r-r1)) {
                r3 = r-r1-r2;
                if (r3<=F34N) {
                    result = rbind(result, c(r1, r2, r3))
                } 
                else break;
            }
        }
        return(result[2:dim(result)[1],])
    }  #test: tripleSet (4, 10, 10)
    
    # Function return all possible triples (r1, r2, r3) 
    #        s.t. r1+r2+r3=r; 0<=r1<=p1, 0<=r2<=p2, 0<=r3<=F34N.
    tripleSet2 = function(r, p1, p2, F34N) {
        result = array(0,3)
        for (r1 in 0:min(p1, r)) {
            for (r2 in 0:min(p2, r-r1)) {
                r3 = r-r1-r2;
                if (r3<=F34N) {
                    result = rbind(result, c(r1, r2, r3))
                } 
                else break;
            }
        }
        return(result[2:dim(result)[1],])
    }  #test: 
    #tripleSet2 (4, 8, 8, 10) #same as above
 
 





 
 ########################################################
 # Formulas of popn distribution params.
 ########################################################
 
#Get the mean of T1 test statistics
meanT1Func = function(b1, b2, b3, p1, p2, sigm, n) {  
    q1=1-p1; q2=1-p2; 
     return( 
     (sqrt(2)*sqrt(n)*p1*q1*(b1 + b3*(p2 - q2)))/
     sqrt(p1*q1*(2*p2*((b2 + b3*p1)^2 - 2*b2*b3*q1 + b3^2*q1^2)*q2 + sigm^2))
     )
}

#Get the variance of T1 test statistics
varT1Func = function(b1, b2, b3, p1, p2, sigm) { 
q1=1-p1; q2=1-p2; 
return(
(-6*b3^6*p1^4*p2^4*q1^2 + 8*b3^6*p1^4*p2^6*q1^2 - 2*b3^6*p1^6*p2^6*q1^2 - 
  4*b3^6*p1^5*p2^6*q1^3 - 6*b3^6*p1^2*p2^4*q1^4 + 8*b3^6*p1^2*p2^6*q1^4 - 
  4*b3^6*p1^4*p2^6*q1^4 - 4*b3^6*p1^3*p2^6*q1^5 - 2*b3^6*p1^2*p2^6*q1^6 + 
  12*b3^6*p1^4*p2^3*q1^2*q2 + 4*b3^6*p1^6*p2^5*q1^2*q2 + 
  8*b3^6*p1^5*p2^5*q1^3*q2 + 12*b3^6*p1^2*p2^3*q1^4*q2 + 
  8*b3^6*p1^4*p2^5*q1^4*q2 + 8*b3^6*p1^3*p2^5*q1^5*q2 + 
  4*b3^6*p1^2*p2^5*q1^6*q2 - 8*b3^6*p1^6*p2^4*q2^2 + 8*b3^6*p1^8*p2^4*q2^2 + 
  32*b3^6*p1^5*p2^4*q1*q2^2 - 12*b3^6*p1^7*p2^4*q1*q2^2 - 
  12*b3^6*p1^4*p2^2*q1^2*q2^2 - 32*b3^6*p1^4*p2^4*q1^2*q2^2 + 
  10*b3^6*p1^6*p2^4*q1^2*q2^2 + 64*b3^6*p1^3*p2^4*q1^3*q2^2 - 
  32*b3^6*p1^5*p2^4*q1^3*q2^2 - 12*b3^6*p1^2*p2^2*q1^4*q2^2 - 
  32*b3^6*p1^2*p2^4*q1^4*q2^2 + 4*b3^6*p1^4*p2^4*q1^4*q2^2 + 
  32*b3^6*p1*p2^4*q1^5*q2^2 - 32*b3^6*p1^3*p2^4*q1^5*q2^2 - 
  8*b3^6*p2^4*q1^6*q2^2 + 10*b3^6*p1^2*p2^4*q1^6*q2^2 - 
  12*b3^6*p1*p2^4*q1^7*q2^2 + 8*b3^6*p2^4*q1^8*q2^2 + 
  32*b2^6*p1*p2^3*q1*q2^3 - 8*b3^6*p1^7*p2^3*q1*q2^3 + 
  12*b3^6*p1^4*p2*q1^2*q2^3 + 40*b3^6*p1^6*p2^3*q1^2*q2^3 - 
  40*b3^6*p1^5*p2^3*q1^3*q2^3 + 12*b3^6*p1^2*p2*q1^4*q2^3 + 
  80*b3^6*p1^4*p2^3*q1^4*q2^3 - 40*b3^6*p1^3*p2^3*q1^5*q2^3 + 
  40*b3^6*p1^2*p2^3*q1^6*q2^3 - 8*b3^6*p1*p2^3*q1^7*q2^3 - 
  8*b3^6*p1^6*p2^2*q2^4 + 8*b3^6*p1^8*p2^2*q2^4 + 32*b3^6*p1^5*p2^2*q1*q2^4 - 
  12*b3^6*p1^7*p2^2*q1*q2^4 - 6*b3^6*p1^4*q1^2*q2^4 - 
  32*b3^6*p1^4*p2^2*q1^2*q2^4 + 10*b3^6*p1^6*p2^2*q1^2*q2^4 + 
  64*b3^6*p1^3*p2^2*q1^3*q2^4 - 32*b3^6*p1^5*p2^2*q1^3*q2^4 - 
  6*b3^6*p1^2*q1^4*q2^4 - 32*b3^6*p1^2*p2^2*q1^4*q2^4 + 
  4*b3^6*p1^4*p2^2*q1^4*q2^4 + 32*b3^6*p1*p2^2*q1^5*q2^4 - 
  32*b3^6*p1^3*p2^2*q1^5*q2^4 - 8*b3^6*p2^2*q1^6*q2^4 + 
  10*b3^6*p1^2*p2^2*q1^6*q2^4 - 12*b3^6*p1*p2^2*q1^7*q2^4 + 
  8*b3^6*p2^2*q1^8*q2^4 + 4*b3^6*p1^6*p2*q1^2*q2^5 + 
  8*b3^6*p1^5*p2*q1^3*q2^5 + 8*b3^6*p1^4*p2*q1^4*q2^5 + 
  8*b3^6*p1^3*p2*q1^5*q2^5 + 4*b3^6*p1^2*p2*q1^6*q2^5 + 
  8*b3^6*p1^4*q1^2*q2^6 - 2*b3^6*p1^6*q1^2*q2^6 - 4*b3^6*p1^5*q1^3*q2^6 + 
  8*b3^6*p1^2*q1^4*q2^6 - 4*b3^6*p1^4*q1^4*q2^6 - 4*b3^6*p1^3*q1^5*q2^6 - 
  2*b3^6*p1^2*q1^6*q2^6 - 16*b2^5*b3*p2^2*(p1 - q1)*q2^2*
   (p2^2*(-1 + p1 + q1)*(1 + p1 + q1) - 8*p1*p2*q1*q2 + 
    (-1 + p1 + q1)*(1 + p1 + q1)*q2^2) + 
  4*b3^4*(p1^5*p2*q1*q2*(-3*p2^2 + 2*p2*q2 - 3*q2^2) + 
    2*p1^6*p2*q2*(p2^2 + q2^2) + 2*p2*q1^4*(-1 + q1^2)*q2*(p2^2 + q2^2) + 
    p1*p2*q1^3*q2*(p2^2*(8 - 3*q1^2) + 2*p2*q1^2*q2 + (8 - 3*q1^2)*q2^2) + 
    2*p1^3*q1*(-2*p2^4*q1^2 + p2^3*(4 + 5*q1^2)*q2 - 10*p2^2*q1^2*q2^2 + 
      p2*(4 + 5*q1^2)*q2^3 - 2*q1^2*q2^4) - 
    2*p1^4*(p2^4*q1^2 + p2*(1 - 4*q1^2)*q2^3 + q1^2*q2^4 + 
      p2^3*(q2 - 4*q1^2*q2)) - 2*p1^2*q1^2*(p2^4*(-2 + q1^2) + 
      p2^3*(6 - 4*q1^2)*q2 + q2^2 + (-2 + q1^2)*q2^4 + p2^2*(1 - 4*q2^2) - 
      2*p2*(q2 + (-3 + 2*q1^2)*q2^3)))*sigm^2 + 
  b3^2*(p2^2*(2*p1^4 - 3*p1^3*q1 + 2*q1^2*(-1 + q1^2) - 2*p1^2*(1 + 4*q1^2) + 
      p1*(8*q1 - 3*q1^3)) + 2*p1*p2*q1*(7*p1^2 + 12*p1*q1 + 7*q1^2)*q2 + 
    (2*p1^4 - 3*p1^3*q1 + 2*q1^2*(-1 + q1^2) - 2*p1^2*(1 + 4*q1^2) + 
      p1*(8*q1 - 3*q1^3))*q2^2)*sigm^4 + 4*p1*q1*sigm^6 + 
  4*b2^4*p2*q2*(b3^2*(-14*p1^4*p2*q2*(p2^2 + q2^2) - 
      14*p2*q1^2*(-1 + q1^2)*q2*(p2^2 + q2^2) + 
      p1^2*(p2^4*q1^2 - 78*p2^2*q1^2*q2^2 + 14*p2*(1 + 3*q1^2)*q2^3 + 
        q1^2*q2^4 + 14*p2^3*(q2 + 3*q1^2*q2)) + 
      p1*(4*p2^4*q1^3 + 8*p2*q1^3*q2 - 3*p2^3*q1*(8 + q1^2)*q2 - 
        3*p2*q1*(8 + q1^2)*q2^3 + 4*q1^3*q2^2*(-1 + q2^2) + 
        2*p2^2*q1^3*(-2 + 19*q2^2)) + p1^3*q1*(4*p2^4 - 3*p2^3*q2 + 
        4*q2^2*(-1 + q2^2) + p2^2*(-4 + 38*q2^2) + p2*(8*q2 - 3*q2^3))) + 
    12*p1*p2*q1*q2*sigm^2) - 8*b2^3*b3*p2*(p1 - q1)*q2*
   (b3^2*(p1*p2^2*q1*(3 + p1^2*(2 + p2^2) + 4*p1*p2^2*q1 + 2*q1^2 + 
        p2^2*(-6 + q1^2)) + 2*p2*(4*p1^4*p2^2 + p1^3*(-2 + p2^2)*q1 + 
        4*p2^2*q1^2*(-1 + q1^2) - 2*p1^2*p2^2*(2 + 11*q1^2) + 
        p1*q1*(-3 - 2*q1^2 + p2^2*(6 + q1^2)))*q2 + 
      p1*q1*(3 + 2*p1^2*(1 + 5*p2^2) + 48*p1*p2^2*q1 + 2*q1^2 + 
        2*p2^2*(-6 + 5*q1^2))*q2^2 + 2*p2*(4*p1^4 + p1^3*q1 + 
        4*q1^2*(-1 + q1^2) + p1*q1*(6 + q1^2) - 2*p1^2*(2 + 11*q1^2))*q2^3 + 
      p1*q1*(-6 + p1^2 + 4*p1*q1 + q1^2)*q2^4) + 
    2*(p2^2*(-1 + p1 + q1)*(1 + p1 + q1) - 8*p1*p2*q1*q2 + 
      (-1 + p1 + q1)*(1 + p1 + q1)*q2^2)*sigm^2) - 
  4*b2^2*(2*b3^4*(2*p1^6*p2^2*q2^2*(p2^2 + q2^2) + 2*p2^2*q1^4*(-1 + q1^2)*
       q2^2*(p2^2 + q2^2) + p1^5*p2*q1*q2*(8*p2^4 - 7*p2^3*q2 - 2*q2^2 + 
        8*q2^4 + p2^2*(-2 + 34*q2^2) + p2*(4*q2 - 7*q2^3)) - 
      p1^4*(3*p2^6*q1^2 + p2^5*q1^2*q2 - 2*p2^3*q1^2*q2*(-3 + q2^2) + 
        3*q1^2*q2^4*(-1 + q2^2) + p2*q1^2*q2^3*(6 + q2^2) + 
        p2^4*(-3*q1^2 + (2 + 39*q1^2)*q2^2) + 
        p2^2*(-6*q1^2*q2^2 + (2 + 39*q1^2)*q2^4)) + 
      p1^2*q1^2*(-3*p2^6*q1^2 - p2^5*(-20 + q1^2)*q2 - 
        3*q1^2*q2^4*(-1 + q2^2) + 3*p2^4*(q1^2 - (4 + 13*q1^2)*q2^2) + 
        p2^2*q2^2*(20 - 12*q2^2 + q1^2*(6 - 39*q2^2)) + 
        2*p2^3*q2*(-5 + 20*q2^2 + q1^2*(-3 + q2^2)) - 
        p2*q2^3*(10 - 20*q2^2 + q1^2*(6 + q2^2))) + 
      2*p1^3*p2*q1*q2*(-2*p2^4*(3 + q1^2) + p2^3*(4 + 29*q1^2)*q2 + 
        q2^2*(3 - 2*q1^2 - 2*(3 + q1^2)*q2^2) - 
        p2^2*(-3 + 2*q1^2 + 6*(2 + 3*q1^2)*q2^2) + 
        p2*q2*(-6 + 4*q2^2 + q1^2*(4 + 29*q2^2))) + 
      p1*p2*q1^3*q2*(4*p2^4*(-3 + 2*q1^2) + p2^3*(8 - 7*q1^2)*q2 - 
        2*(-3 + q1^2)*q2^2 + 4*(-3 + 2*q1^2)*q2^4 + 
        p2*q2*(-12 + 8*q2^2 + q1^2*(4 - 7*q2^2)) + 
        p2^2*(6 - 24*q2^2 + q1^2*(-2 + 34*q2^2)))) + 
    b3^2*(6*p1^4*p2*q2*(p2^2 + q2^2) + 6*p2*q1^2*(-1 + q1^2)*q2*
       (p2^2 + q2^2) - 2*p1^2*p2*q2*(p2^2*(3 + 10*q1^2) - 16*p2*q1^2*q2 + 
        (3 + 10*q1^2)*q2^2) + p1*(-2*p2^4*q1^3 - 4*p2*q1^3*q2 + 
        p2^3*q1*(8 + 3*q1^2)*q2 + p2*q1*(8 + 3*q1^2)*q2^3 + 
        2*p2^2*q1^3*(1 - 13*q2^2) - 2*q1^3*q2^2*(-1 + q2^2)) + 
      p1^3*q1*(-2*p2^4 + 3*p2^3*q2 + 2*q2^2 - 2*q2^4 + p2^2*(2 - 26*q2^2) + 
        p2*q2*(-4 + 3*q2^2)))*sigm^2 - 6*p1*p2*q1*q2*sigm^4) + 
  4*b2*b3*(p1 - q1)*(2*b3^4*(2*p1^6*p2^2*q2^2*(p2^2 + q2^2) + 
      2*p2^2*q1^4*(-1 + q1^2)*q2^2*(p2^2 + q2^2) - p1^4*(p2^2 + q2^2)*
       (p2^4*q1^2 - 2*p2^2*(-1 + 4*q1^2)*q2^2 + q1^2*q2^4) + 
      p1^2*q1^2*(-(p2^4*(3 + p2^2*(-4 + q1^2))) - 2*p2^3*(-5 + 4*p2^2)*q2 + 
        p2^2*(-14 + p2^2*(8 + 7*q1^2))*q2^2 - 2*p2*(-5 + 8*p2^2)*q2^3 + 
        (-3 + p2^2*(8 + 7*q1^2))*q2^4 - 8*p2*q2^5 - (-4 + q1^2)*q2^6) + 
      p1^5*p2*q1*q2*(-5*p2^4 + 2*p2^3*q2 + 2*q2^2 - 5*q2^4 + 
        p2^2*(2 - 18*q2^2) + 2*p2*q2*(-2 + q2^2)) + 
      p1*p2*q1^3*q2*(p2^4*(6 - 5*q1^2) + 2*p2^3*(2 + q1^2)*q2 + 
        2*p2*q2*(3 - 2*q1^2 + (2 + q1^2)*q2^2) + 
        p2^2*(-3 + 2*q1^2 - 6*(-2 + 3*q1^2)*q2^2) + 
        q2^2*(-3 + 6*q2^2 + q1^2*(2 - 5*q2^2))) + 
      p1^3*q1*(-2*p2^6*q1^2 + 2*p2^5*(3 + q1^2)*q2 - 2*p2^4*(-2 + 13*q1^2)*
         q2^2 - 2*q1^2*q2^6 + p2*q2^3*(-3 + 4*q1^2 + 2*(3 + q1^2)*q2^2) + 
        p2^3*q2*(-3 + 4*q1^2 + 4*(3 + q1^2)*q2^2) + 
        2*p2^2*q2^2*(3 + 2*q2^2 - q1^2*(4 + 13*q2^2)))) - 
    b3^2*p1*q1*(p2^4*(-6 + 5*p1^2 + 6*p1*q1 + 5*q1^2) + 
      2*p2^3*(-2 + p1^2 - 8*p1*q1 + q1^2)*q2 + 
      2*p2*q2*(-3 + 2*p1^2 + 2*q1^2 + (-2 + p1^2 - 8*p1*q1 + q1^2)*q2^2) + 
      p2^2*(3 - 2*p1^2 - 2*q1^2 + 2*(-6 + p1^2 + 10*p1*q1 + q1^2)*q2^2) + 
      q2^2*(3 - 2*p1^2 - 2*q1^2 + (-6 + 5*p1^2 + 6*p1*q1 + 5*q1^2)*q2^2))*
     sigm^2 - (p2^2*(-1 + p1 + q1)*(1 + p1 + q1) - 8*p1*p2*q1*q2 + 
      (-1 + p1 + q1)*(1 + p1 + q1)*q2^2)*sigm^4) + 
  b1^2*p1*q1*(4*b2^4*p2*q2*(p2*q1 + p1*q2)*(p1*p2 + q1*q2) + 
    2*b3^4*(p1^2 + q1^2)*
     (-(p1*p2^2*q1*(3 + p2^2*(-2 + p1 + q1)*(2 + p1 + q1))) + 
      8*p1*p2*q1*(1 + p2^2*(-1 + p1 + q1)*(1 + p1 + q1))*q2 + 
      (2*p1^4*p2^2 - 10*p1^3*p2^2*q1 - 56*p1^2*p2^2*q1^2 + 2*p2^2*q1^4 + 
        p1*q1*(-3 + 2*p2^2*(4 - 5*q1^2)))*q2^2 + 8*p1*p2*q1*(-1 + p1 + q1)*
       (1 + p1 + q1)*q2^3 - p1*q1*(-2 + p1 + q1)*(2 + p1 + q1)*q2^4) + 
    16*b2^3*b3*p2*(p1 - q1)*q2*(p1^2*p2*q2 + p2*q1^2*q2 + 
      p1*q1*(p2 + q2)^2) + 4*b3^2*p2*(p1^4 + 6*p1^3*q1 - 6*p1^2*q1^2 + 
      6*p1*q1^3 + q1^4)*q2*sigm^2 + (p1^2 + 4*p1*q1 + q1^2)*sigm^4 + 
    4*b2^2*(2*b3^2*(3*p1*p2^2*(-1 + p2^2)*q1*(p1^2 + q1^2) + 
        p1*p2*(8 + p2^2)*q1*(p1^2 + q1^2)*q2 + 
        (3*p1^4*p2^2 - p1^3*(3 + 4*p2^2)*q1 - 10*p1^2*p2^2*q1^2 - 
          p1*(3 + 4*p2^2)*q1^3 + 3*p2^2*q1^4)*q2^2 + p1*p2*q1*(p1^2 + q1^2)*
         q2^3 + 3*p1*q1*(p1^2 + q1^2)*q2^4) + p2*(p1^2 + 4*p1*q1 + q1^2)*q2*
       sigm^2) + 8*b2*b3*(p1 - q1)*
     (b3^2*(-(p1*p2^2*q1*(3 + p2^2*(-2 + p1 + q1)*(2 + p1 + q1))) + 
        8*p1*p2*q1*(1 + p2^2*(-1 + p1 + q1)*(1 + p1 + q1))*q2 + 
        (2*p1^4*p2^2 - 10*p1^3*p2^2*q1 - 40*p1^2*p2^2*q1^2 + 2*p2^2*q1^4 + 
          p1*q1*(-3 + 2*p2^2*(4 - 5*q1^2)))*q2^2 + 8*p1*p2*q1*(-1 + p1 + q1)*
         (1 + p1 + q1)*q2^3 - p1*q1*(-2 + p1 + q1)*(2 + p1 + q1)*q2^4) + 
      p2*(p1^2 + 6*p1*q1 + q1^2)*q2*sigm^2)) + 
  2*b1*b3*p1*q1*(p2 - q2)*(2*b3^4*(p1^2 + q1^2)*
     (-(p1*p2^2*q1*(3 + p2^2*(-2 + p1 + q1)*(2 + p1 + q1))) + 
      4*p1*p2*q1*(1 + p2^2*(p1 + q1)^2)*q2 + 
      (2*p1^4*p2^2 + 6*p1^3*p2^2*q1 - 24*p1^2*p2^2*q1^2 + 2*p2^2*q1^4 + 
        p1*q1*(-3 + p2^2*(8 + 6*q1^2)))*q2^2 + 4*p1*p2*q1*(p1 + q1)^2*q2^3 - 
      p1*q1*(-2 + p1 + q1)*(2 + p1 + q1)*q2^4) - 4*b2^3*b3*p2*(p1 - q1)*q2*
     (3 - 6*p2^2 - 6*q2^2 + 2*p1*q1*(p2^2 - 28*p2*q2 + q2^2) + 
      p1^2*(2 + p2^2 - 14*p2*q2 + q2^2) + q1^2*(2 + p2^2 - 14*p2*q2 + 
        q2^2)) + 4*b2^4*p2*q2*(p1*q1*(p2^2 + 12*p2*q2 + q2^2) + 
      p1^2*(-2 + (2*p2 + q2)*(p2 + 2*q2)) + 
      q1^2*(-2 + (2*p2 + q2)*(p2 + 2*q2))) + 
    4*b3^2*(-(p1*q1*(1 + p2^2*(-2 + (p1 + q1)^2))) + 
      p2*(p1^2 + q1^2)*(p1^2 + 10*p1*q1 + q1^2)*q2 - p1*q1*(-2 + (p1 + q1)^2)*
       q2^2)*sigm^2 + (p1^2 + 4*p1*q1 + q1^2)*sigm^4 - 
    4*b2^2*(2*b3^2*(p1^4*p2*q2*(-1 + 4*p2^2 - 7*p2*q2 + 4*q2^2) + 
        p1^2*p2*q2*(3 - 2*q1^2 - 2*p2^2*(3 + q1^2) + 42*p2*q1^2*q2 - 
          2*(3 + q1^2)*q2^2) + p2*q1^2*q2*(3 + p2^2*(-6 + 4*q1^2) - 
          7*p2*q1^2*q2 - 6*q2^2 + q1^2*(-1 + 4*q2^2)) - 
        p1^3*q1*(3*p2^4 + 4*p2^3*q2 + 3*q2^2*(-1 + q2^2) + 
          p2^2*(-3 + 26*q2^2) + 4*p2*(q2 + q2^3)) + 
        p1*(-3*p2^4*q1^3 + 2*p2^3*q1*(5 - 2*q1^2)*q2 + 
          p2^2*q1^3*(3 - 26*q2^2) - 3*q1^3*q2^2*(-1 + q2^2) - 
          p2*q1*q2*(5 - 10*q2^2 + 4*q1^2*(1 + q2^2)))) - 
      ((-1 + p2^2)*(p1^2 + q1^2) + p2*(3*p1 + q1)*(p1 + 3*q1)*q2 + 
        (p1^2 + q1^2)*q2^2)*sigm^2) - 2*b2*b3*(p1 - q1)*
     (2*b3^2*(2*p1*p2^2*q1*(3 + p2^2*(-2 + p1 + q1)*(2 + p1 + q1)) + 
        p2*(p1^4*(-2 + 5*p2^2) - 6*p1^3*p2^2*q1 - 
          2*p1*q1*(6 + p2^2*(-4 + 3*q1^2)) + q1^2*(3 - 2*q1^2 + 
            p2^2*(-6 + 5*q1^2)) - p1^2*(-3 + 4*q1^2 + 2*p2^2*(3 + 7*q1^2)))*
         q2 - 2*(3*p1^4*p2^2 + 10*p1^3*p2^2*q1 - 22*p1^2*p2^2*q1^2 + 
          3*p2^2*q1^4 + p1*q1*(-3 + 2*p2^2*(4 + 5*q1^2)))*q2^2 + 
        p2*(5*p1^4 - 6*p1^3*q1 + q1^2*(-6 + 5*q1^2) - 2*p1^2*(3 + 7*q1^2) + 
          p1*(8*q1 - 6*q1^3))*q2^3 + 2*p1*q1*(-2 + p1 + q1)*(2 + p1 + q1)*
         q2^4) + (3 - 6*p2^2 - 6*q2^2 + 6*p1*q1*(p2^2 - 8*p2*q2 + q2^2) + 
        p1^2*(-2 + 5*p2^2 - 6*p2*q2 + 5*q2^2) + 
        q1^2*(-2 + 5*p2^2 - 6*p2*q2 + 5*q2^2))*sigm^2)))/
 (4*p1*q1*(2*p2*((b2 + b3*p1)^2 - 2*b2*b3*q1 + b3^2*q1^2)*q2 + sigm^2)^3)
)
}


#Get the covariance b/w T1 and T2 test statistics
covT1T2Func = function(b1, b2, b3, p1, p2, sigm) { 
q1=1-p1; q2=1-p2; 
return(
(p1*p2*q1*q2*(-2*b1^4*b3*p1^2*q1^2*(p2 - q2)*
    (b3*(p1 - q1)*(1 + p2^2*(-2 + p1^2 + 2*p1*q1 + q1^2) - 
       2*p2*(p1^2 + 4*p1*q1 + q1^2)*q2 + (-2 + p1^2 + 2*p1*q1 + q1^2)*q2^2) - 
     2*b2*(4*p1*p2*q1*q2 + p1^2*(-1 + p2^2 + 2*p2*q2 + q2^2) + 
       q1^2*(-1 + p2^2 + 2*p2*q2 + q2^2))) - 
   b1^3*p1*q1*(2*b2*b3^2*(p1^4*p2*q2*(4*p2^2 + 3*p2*q2 + 4*q2^2) + 
       p2*q1^2*q2*(4*p2^2*(-1 + q1^2) + 3*p2*q1^2*q2 + 4*(-1 + q1^2)*q2^2) - 
       2*p1^2*p2*q2*(2*p2^2*(1 + 9*q1^2) - 27*p2*q1^2*q2 + 
         2*(1 + 9*q1^2)*q2^2) - p1^3*q1*(2*p2^4 + p2^3*q2 + p2*q2^3 + 
         2*p2^2*(-1 + q2^2) + 2*q2^2*(-1 + q2^2)) - 
       p1*q1*(2*p2^4*q1^2 + p2^3*(-8 + q1^2)*q2 + p2*(-8 + q1^2)*q2^3 + 
         2*p2^2*q1^2*(-1 + q2^2) + 2*q1^2*q2^2*(-1 + q2^2))) + 
     2*b3^3*(p1 - q1)*(p1^4*p2*q2*(2*p2^2 + p2*q2 + 2*q2^2) + 
       p2*q1^2*q2*(2*p2^2*(-1 + q1^2) + p2*q1^2*q2 + 2*(-1 + q1^2)*q2^2) + 
       p1^3*q1*(3*p2^4 - 5*p2^3*q2 + 12*p2^2*q2^2 - 5*p2*q2^3 + 3*q2^4) + 
       p1^2*(6*p2^4*q1^2 + 46*p2^2*q1^2*q2^2 - 2*p2*(1 + 14*q1^2)*q2^3 + 
         6*q1^2*q2^4 - 2*p2^3*(q2 + 14*q1^2*q2)) + 
       p1*q1*(p2^4*(-4 + 3*q1^2) + p2^3*(8 - 5*q1^2)*q2 + q2^2 + 
         p2*(8 - 5*q1^2)*q2^3 + (-4 + 3*q1^2)*q2^4 + 
         p2^2*(1 + 4*(-2 + 3*q1^2)*q2^2))) + 
     b2*p2*q2*(2*b2^2*(p1^2*p2*q2 + p2*q1^2*q2 + 
         p1*q1*(p2^2 - 8*p2*q2 + q2^2)) + (p1^2 + q1^2)*sigm^2) + 
     b3*(p1 - q1)*(2*b2^2*p2*q2*(p2^2*(-2 + 2*p1^2 + 7*p1*q1 + 2*q1^2) + 
         p2*(3*p1^2 - 8*p1*q1 + 3*q1^2)*q2 + (-2 + 2*p1^2 + 7*p1*q1 + 2*q1^2)*
          q2^2) + (2*p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2) + 
         p2*(p1^2 + 12*p1*q1 + q1^2)*q2 + 2*(-1 + p1^2 + 2*p1*q1 + q1^2)*
          q2^2)*sigm^2)) - b1^2*b3*p1*q1*(p2 - q2)*
    (2*b3^3*(p1 - q1)*(p1^4*p2*q2*(2 + 3*p2^2 + p2*q2 + 3*q2^2) + 
       p1^3*q1*(3*p2^4 + 5*p2^3*q2 + 12*p2^2*q2^2 + 5*p2*q2^3 + 3*q2^4) + 
       p2*q1^2*q2*(-1 + p2^2*(-4 + 3*q1^2) + p2*q1^2*q2 - 4*q2^2 + 
         q1^2*(2 + 3*q2^2)) + p1*q1*(3*p2^4*q1^2 + p2^3*(4 + 5*q1^2)*q2 + 
         3*q2^2*(-1 + q1^2*q2^2) + 3*p2^2*(-1 + 4*q1^2*q2^2) + 
         p2*q2*(10 + (4 + 5*q1^2)*q2^2)) + 
       p1^2*(6*p2^4*q1^2 - 2*p2^3*(2 + 11*q1^2)*q2 + 30*p2^2*q1^2*q2^2 + 
         6*q1^2*q2^4 - p2*(q2 - 4*q1^2*q2 + (4 + 22*q1^2)*q2^3))) + 
     2*b2*b3^2*(p1^4*p2*q2*(6 + p2^2 + 7*p2*q2 + q2^2) + 
       p2*q1^2*q2*(-5 + p2^2*(-2 + q1^2) + 7*p2*q1^2*q2 - 2*q2^2 + 
         q1^2*(6 + q2^2)) + p1^3*q1*(6*p2^4 + 7*p2^3*q2 + 
         6*q2^2*(-1 + q2^2) + 6*p2^2*(-1 + 3*q2^2) + p2*q2*(12 + 7*q2^2)) - 
       p1^2*p2*q2*(5 + p2^2*(2 + 66*q1^2) - 14*p2*q1^2*q2 + 2*q2^2 + 
         6*q1^2*(-2 + 11*q2^2)) + p1*q1*(6*p2^4*q1^2 + p2^3*(8 + 7*q1^2)*q2 + 
         6*q1^2*q2^2*(-1 + q2^2) + 6*p2^2*q1^2*(-1 + 3*q2^2) + 
         p2*q2*(8*(1 + q2^2) + q1^2*(12 + 7*q2^2)))) + 
     b2*p2*q2*(2*b2^2*(p1^2*(-2 + 2*p2^2 + 7*p2*q2 + 2*q2^2) + 
         q1^2*(-2 + 2*p2^2 + 7*p2*q2 + 2*q2^2) + 
         p1*q1*(3*p2^2 - 8*p2*q2 + 3*q2^2)) + (3*p1^2 + 4*p1*q1 + 3*q1^2)*
        sigm^2) + b3*(p1 - q1)*(2*b2^2*p2*q2*(2*p2^2 + 13*p2*q1^2*q2 + 
         p1^2*(2 + 13*p2*q2) + 2*(-2 + q1^2 + q2^2) + 
         p1*q1*(13*p2^2 + 24*p2*q2 + 13*q2^2)) + 
       (-1 + 2*q1^2 + p2^2*(-4 + 3*q1^2) + p2*q1^2*q2 - 4*q2^2 + 
         3*q1^2*q2^2 + 10*p1*q1*(p2^2 + 4*p2*q2 + q2^2) + 
         p1^2*(2 + 3*p2^2 + p2*q2 + 3*q2^2))*sigm^2)) + 
   b1*(-(b3^5*(p1 - q1)*(8*p1^6*p2^2*q2^2*(p2^2 + q2^2) + 
        8*p2^2*q1^4*(-1 + q1^2)*q2^2*(p2^2 + q2^2) + 
        2*p1^4*(p2^6*q1^2 + 6*p2^5*q1^2*q2 - p2^4*(4 + q1^2)*q2^2 + 
          20*p2^3*q1^2*q2^3 - p2^2*(4 + q1^2)*q2^4 + 6*p2*q1^2*q2^5 + 
          q1^2*q2^6) - p1^5*p2*q1*q2*(p2^4 - 4*p2^3*q2 + q2^2*(-8 + q2^2) - 
          4*p2*q2*(-4 + q2^2) + 2*p2^2*(-4 + 5*q2^2)) + 
        2*p1^2*q1^2*(p2^6*(4 + q1^2) + 2*p2^5*(-4 + 3*q1^2)*q2 - 
          p2^4*(5 + (-4 + q1^2)*q2^2) - p2^2*q2^2*(34 + (-4 + q1^2)*q2^2) + 
          q2^4*(-5 + (4 + q1^2)*q2^2) + 2*p2^3*q2*(11 + 2*(-4 + 5*q1^2)*
             q2^2) + p2*(22*q2^3 + (-8 + 6*q1^2)*q2^5)) + 
        p1*p2*q1^3*q2*(-(p2^4*q1^2) + 4*p2^3*(4 + q1^2)*q2 + 
          p2^2*(-7 + q1^2*(8 - 10*q2^2)) - q2^2*(7 + q1^2*(-8 + q2^2)) + 
          4*p2*q2*(q1^2*(-4 + q2^2) + 4*(1 + q2^2))) + 
        p1^3*q1*(4*p2^6*q1^2 - 26*p2^5*q1^2*q2 + 4*p2^4*(4 + 13*q1^2)*q2^2 + 
          4*q1^2*q2^6 + p2^3*q2*(-7 + q1^2*(16 - 68*q2^2)) + 
          p2*q2^3*(-7 + q1^2*(16 - 26*q2^2)) + 4*p2^2*q2^2*
           (4*(1 + q2^2) + q1^2*(-8 + 13*q2^2))))) + 
     b2*b3^4*(-20*p1^6*p2^2*q2^2*(p2^2 + q2^2) - 20*p2^2*q1^4*(-1 + q1^2)*
        q2^2*(p2^2 + q2^2) + p1^5*p2*q1*q2*(7*p2^4 + 4*p2^3*q2 + 
         4*p2*q2*(4 + q2^2) + q2^2*(-12 + 7*q2^2) + 2*p2^2*(-6 + 43*q2^2)) - 
       4*p1^4*(5*p2^6*q1^2 - p2^5*q1^2*q2 - p2*q1^2*q2^3*(-14 + q2^2) + 
         5*q1^2*q2^4*(-1 + q2^2) + 2*p2^3*q1^2*q2*(7 + 13*q2^2) + 
         p2^2*(-5*q2^4 + 2*q1^2*q2^2*(-9 + 4*q2^2)) + 
         p2^4*(-5*q2^2 + q1^2*(-5 + 8*q2^2))) - 4*p1^2*q1^2*
        (5*p2^6*q1^2 - p2^5*(4 + q1^2)*q2 + 5*q1^2*q2^4*(-1 + q2^2) + 
         p2*q2^3*(8 - 4*q2^2 - q1^2*(-14 + q2^2)) + 2*p2^2*q2^2*
          (-8 - 9*q2^2 + q1^2*(-9 + 4*q2^2)) + 
         p2^4*(-18*q2^2 + q1^2*(-5 + 8*q2^2)) + 2*p2^3*q2*
          (4 - 4*q2^2 + q1^2*(7 + 13*q2^2))) + p1*p2*q1^3*q2*
        (p2^4*(-12 + 7*q1^2) + 4*p2^3*(-14 + q1^2)*q2 + 
         q2^2*(17 - 12*q2^2 + q1^2*(-12 + 7*q2^2)) + 
         4*p2*q2*(q1^2*(4 + q2^2) - 2*(4 + 7*q2^2)) + 
         p2^2*(17 - 24*q2^2 + 2*q1^2*(-6 + 43*q2^2))) + 
       p1^3*p2*q1*q2*(2*p2^4*(-6 + 43*q1^2) - 8*p2^3*(7 + 13*q1^2)*q2 - 
         8*p2*q2*(4 + 7*q2^2 + q1^2*(-4 + 13*q2^2)) + 
         p2^2*(17 - 24*q2^2 + 4*q1^2*(-6 + 79*q2^2)) + 
         q2^2*(17 - 12*q2^2 + q1^2*(-24 + 86*q2^2)))) - 
     b2*p1*p2*q1*q2*sigm^2*(b2^2*(p2^2 + q2^2) + 2*sigm^2) + 
     2*b2*b3^2*(b2^2*p2*q2*(2*p1^4*p2*q2*(p2^2 + q2^2) + 
         2*p2*q1^2*(-1 + q1^2)*q2*(p2^2 + q2^2) - 
         p1^2*(3*p2^4*q1^2 - 2*p2^3*(-1 + q1^2)*q2 + 54*p2^2*q1^2*q2^2 - 
           2*p2*(-1 + q1^2)*q2^3 + 3*q1^2*q2^4) + 
         p1^3*q1*(-4*p2^4 + p2^3*q2 + p2*q2*(-8 + q2^2) - 
           4*q2^2*(-1 + q2^2) + 4*p2^2*(1 + 9*q2^2)) + 
         p1*q1^3*(-4*p2^4 + p2^3*q2 + p2*q2*(-8 + q2^2) - 
           4*q2^2*(-1 + q2^2) + 4*p2^2*(1 + 9*q2^2))) + 
       (-3*p1^4*p2*q2*(p2^2 + q2^2) - 3*p2*q1^2*(-1 + q1^2)*q2*
          (p2^2 + q2^2) + p1^2*p2*q2*(p2^2*(3 + 4*q1^2) + 8*p2*q1^2*q2 + 
           (3 + 4*q1^2)*q2^2) + p1^3*q1*(-3*p2^4 + p2^3*q2 + 
           p2*q2*(-8 + q2^2) - 3*q2^2*(-1 + q2^2) + p2^2*(3 + 4*q2^2)) + 
         p1*(-3*p2^4*q1^3 + p2^3*q1*(-8 + q1^2)*q2 - 3*q1^3*q2^2*
            (-1 + q2^2) + p2^2*q1^3*(3 + 4*q2^2) + 
           p2*(-8*q1*q2^3 + q1^3*q2*(-8 + q2^2))))*sigm^2) - 
     b3^3*(p1 - q1)*(2*b2^2*p2*q2*(6*p1^4*p2*q2*(p2^2 + q2^2) + 
         6*p2*q1^2*(-1 + q1^2)*q2*(p2^2 + q2^2) + p1^2*(p2^2 + q2^2)*
          (7*p2^2*q1^2 + 6*p2*(-1 + 3*q1^2)*q2 + 7*q1^2*q2^2) + 
         p1^3*q1*(p2^4 + 7*p2^3*q2 + q2^2*(-2 + q2^2) + p2*q2*(8 + 7*q2^2) - 
           2*p2^2*(1 + 33*q2^2)) + p1*q1*(p2^4*(6 + q1^2) + 
           p2^3*(12 + 7*q1^2)*q2 + q2^2*(-5 + 6*q2^2 + q1^2*(-2 + q2^2)) + 
           p2*q2*(8 + 12*q2^2 + q1^2*(8 + 7*q2^2)) - 
           p2^2*(5 - 12*q2^2 + q1^2*(2 + 66*q2^2)))) + 
       (6*p1^4*p2*q2*(p2^2 + q2^2) + 6*p2*q1^2*(-1 + q1^2)*q2*(p2^2 + q2^2) + 
         p1^2*(6*p2^4*q1^2 + 2*p2^3*(-3 + 20*q1^2)*q2 - 60*p2^2*q1^2*q2^2 + 
           2*p2*(-3 + 20*q1^2)*q2^3 + 6*q1^2*q2^4) - 
         p1^3*q1*(p2^4 - 10*p2^3*q2 + q2^2*(-4 + q2^2) + 2*p2^2*(-2 + q2^2) + 
           p2*(8*q2 - 10*q2^3)) + p1*q1*(-(p2^4*(-2 + q1^2)) + 
           2*p2^3*(2 + 5*q1^2)*q2 + q2^2*(-5 + 2*q2^2 - q1^2*(-4 + q2^2)) + 
           p2^2*(-5 + 4*q2^2 - 2*q1^2*(-2 + q2^2)) + 
           2*p2*q2*(2*(3 + q2^2) + q1^2*(-4 + 5*q2^2))))*sigm^2) + 
     b3*(p1 - q1)*(4*b2^4*p2^2*q2^2*(p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2) + 
         4*p1*p2*q1*q2 + (-1 + p1^2 + 2*p1*q1 + q1^2)*q2^2) - 
       b2^2*p1*p2*q1*q2*(3*p2^2 + 4*p2*q2 + 3*q2^2)*sigm^2 - 
       (p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2) + 8*p1*p2*q1*q2 + 
         (-1 + p1^2 + 2*p1*q1 + q1^2)*q2^2)*sigm^4)) - 
   b3*(p2 - q2)*(b3^5*(p1 - q1)*(4*p1^6*p2^2*q2^2 + 4*p2^2*q1^4*(-1 + q1^2)*
        q2^2 + p1^5*p2*q1*q2*(-3*p2^4 + 4*p2^3*q2 + 4*q2^2 - 3*q2^4 + 
         p2^2*(4 - 6*q2^2) + 4*p2*q2*(-2 + q2^2)) + 
       4*p1^4*p2*q2*(p2^4*q1^2 - 2*p2^3*q1^2*q2 + 2*p2^2*q1^2*q2^2 + 
         q1^2*q2^4 - p2*(q2 - 3*q1^2*q2 + 2*q1^2*q2^3)) + 
       4*p1^2*q1^2*(p2^6 + p2^5*(-2 + q1^2)*q2 + q2^4*(-1 + q2^2) + 
         p2^4*(-1 + (3 - 2*q1^2)*q2^2) + 2*p2^3*q2*(2 + (-2 + q1^2)*q2^2) + 
         p2*q2^3*(4 + (-2 + q1^2)*q2^2) + p2^2*q2^2*(-8 + 3*q2^2 + 
           q1^2*(3 - 2*q2^2))) + p1^3*p2*q1*q2*(p2^4*(4 - 6*q1^2) + 
         8*p2^3*q1^2*q2 + p2^2*(-5 + 8*q2^2 + q1^2*(8 - 12*q2^2)) + 
         q2^2*(-5 + 4*q2^2 + q1^2*(8 - 6*q2^2)) + 
         8*p2*q2*(2 + q1^2*(-2 + q2^2))) + p1*p2*q1^3*q2*
        (p2^4*(4 - 3*q1^2) + 4*p2^3*q1^2*q2 + 
         p2^2*(-5 + 8*q2^2 + q1^2*(4 - 6*q2^2)) + 
         q2^2*(-5 + 4*q2^2 + q1^2*(4 - 3*q2^2)) + 
         4*p2*q2*(4 + q1^2*(-2 + q2^2)))) + 
     b2*b3^4*(2*p1^6*p2^2*q2^2*(4 + p2^2 + 2*p2*q2 + q2^2) + 
       2*p2^2*q1^4*q2^2*(-5 + q1^2*(4 + p2^2 + 2*p2*q2 + q2^2)) - 
       p1^5*p2*q1*q2*(p2^4 - 12*p2^3*q2 + 26*p2^2*q2^2 + q2^4 - 
         4*p2*q2*(-4 + 3*q2^2)) + 2*p1^4*(4*p2^6*q1^2 + 2*p2^5*q1^2*q2 + 
         4*q1^2*q2^4*(-1 + q2^2) - p2^4*q1^2*(4 + q2^2) + 
         2*p2*q1^2*q2^3*(4 + q2^2) + 2*p2^3*q1^2*q2*(4 + 13*q2^2) - 
         p2^2*q2^2*(5 + q1^2*(-4 + q2^2))) - p1*p2*q1^3*q2*
        (p2^4*(-8 + q1^2) - 12*p2^3*q1^2*q2 + q2^2*(7 + (-8 + q1^2)*q2^2) + 
         p2^2*(7 + 2*(-8 + 13*q1^2)*q2^2) - 4*p2*q2*
          (11 + q1^2*(-4 + 3*q2^2))) + p1^3*p2*q1*q2*(p2^4*(8 - 10*q1^2) + 
         40*p2^3*q1^2*q2 - 7*q2^2 + (8 - 10*q1^2)*q2^4 + 
         p2^2*(-7 + (16 - 68*q1^2)*q2^2) + 4*p2*q2*
          (11 + 2*q1^2*(-4 + 5*q2^2))) + 2*p1^2*q1^2*(4*p2^6*q1^2 + 
         2*p2^5*(-4 + q1^2)*q2 + 4*q1^2*q2^4*(-1 + q2^2) - 
         p2^4*q1^2*(4 + q2^2) - p2^2*q2^2*(34 + q1^2*(-4 + q2^2)) + 
         2*p2*q2^3*(4 - 4*q2^2 + q1^2*(4 + q2^2)) + 
         2*p2^3*q2*(4 - 8*q2^2 + q1^2*(4 + 13*q2^2)))) + 
     b2*sigm^2*(b2^2*p2*q2*(2*p1^2*(-1 + p2^2 + 2*p2*q2 + q2^2) + 
         2*q1^2*(-1 + p2^2 + 2*p2*q2 + q2^2) + p1*q1*(p2^2 + 12*p2*q2 + 
           q2^2)) + (8*p1*p2*q1*q2 + p1^2*(-1 + p2^2 + 2*p2*q2 + q2^2) + 
         q1^2*(-1 + p2^2 + 2*p2*q2 + q2^2))*sigm^2) + 
     b2*b3^2*(2*b2^2*p2*q2*(p1^4*p2*q2*(-4 + 3*p2^2 + 6*p2*q2 + 3*q2^2) + 
         p2*q1^2*q2*(1 + q1^2*(-4 + 3*p2^2 + 6*p2*q2 + 3*q2^2)) + 
         p1^3*q1*(2*p2^4 - 5*p2^3*q2 + 2*q2^2*(-1 + q2^2) - 
           2*p2^2*(1 + 14*q2^2) + p2*(8*q2 - 5*q2^3)) + 
         p1*q1^3*(2*p2^4 - 5*p2^3*q2 + 2*q2^2*(-1 + q2^2) - 
           2*p2^2*(1 + 14*q2^2) + p2*(8*q2 - 5*q2^3)) + 
         p1^2*(p2^4*q1^2 + 12*p2^3*q1^2*q2 + 46*p2^2*q1^2*q2^2 + q1^2*q2^4 + 
           p2*(q2 - 8*q1^2*q2 + 12*q1^2*q2^3))) + 
       (-(p1^4*p2*q2*(-2 + p2^2 - 6*p2*q2 + q2^2)) + 
         p1^2*p2*q2*(-5 - 2*p2^2*(-2 + q1^2) - 60*p2*q1^2*q2 + 4*q2^2 - 
           2*q1^2*(-2 + q2^2)) + p2*q1^2*q2*(-5 - p2^2*(-4 + q1^2) + 
           6*p2*q1^2*q2 + 4*q2^2 - q1^2*(-2 + q2^2)) + 
         2*p1^3*q1*(3*p2^4 + 5*p2^3*q2 + 3*q2^2*(-1 + q2^2) + 
           p2*q2*(2 + 5*q2^2) + p2^2*(-3 + 20*q2^2)) + 
         2*p1*q1*(3*p2^4*q1^2 + p2^3*(-4 + 5*q1^2)*q2 + 
           3*q1^2*q2^2*(-1 + q2^2) + p2^2*q1^2*(-3 + 20*q2^2) + 
           p2*q2*(6 - 4*q2^2 + q1^2*(2 + 5*q2^2))))*sigm^2) + 
     2*b3^3*(p1 - q1)*(b2^2*p2*q2*(3*p1^4*p2*q2*(p2 + q2)^2 + 
         3*p2*q1^2*q2*(-1 + p2^2*q1^2 + 2*p2*q1^2*q2 + q1^2*q2^2) + 
         p1^3*q1*(3*p2^4 + 5*p2^3*q2 + q2^2*(-4 + 3*q2^2) + 
           p2*q2*(4 + 5*q2^2) - 2*p2^2*(2 + 11*q2^2)) + 
         p1^2*(p2^4*q1^2 + 12*p2^3*q1^2*q2 + 30*p2^2*q1^2*q2^2 + q1^2*q2^4 + 
           3*p2*q2*(-1 + 4*q1^2*q2^2)) + p1*q1*(p2^4*(2 + 3*q1^2) + 
           5*p2^3*q1^2*q2 + q2^2*(-1 + 2*q2^2 + q1^2*(-4 + 3*q2^2)) + 
           p2*q2*(10 + q1^2*(4 + 5*q2^2)) - p2^2*(1 - 4*q2^2 + 
             q1^2*(4 + 22*q2^2)))) + (-(p1^4*p2*q2*(-2 + p2^2 + q2^2)) - 
         p2*q1^2*(-1 + q1^2)*q2*(-2 + p2^2 + q2^2) + 
         p1^3*q1*(-p2^4 + 3*p2^3*q2 + q2^2 - q2^4 + p2*q2*(-2 + 3*q2^2) + 
           p2^2*(1 + 6*q2^2)) + p1^2*p2*q2*(-2 + p2^2*(1 + 6*q1^2) - 
           16*p2*q1^2*q2 + q2^2 + q1^2*(4 + 6*q2^2)) + 
         p1*q1*(-(p2^4*(-2 + q1^2)) + p2^3*(-2 + 3*q1^2)*q2 - 
           (-2 + q1^2)*q2^2*(-1 + q2^2) + p2*q2*(6 - 2*q2^2 + 
             q1^2*(-2 + 3*q2^2)) + p2^2*(-2 + 4*q2^2 + q1^2*(1 + 6*q2^2))))*
        sigm^2) + b3*(p1 - q1)*(2*b2^4*p2^2*q2^2*
        (1 + p1^2*(-2 + p2^2 + 2*p2*q2 + q2^2) + 
         q1^2*(-2 + p2^2 + 2*p2*q2 + q2^2) - 2*p1*q1*(p2^2 + 4*p2*q2 + 
           q2^2)) + b2^2*p2*q2*(-1 - 4*q1^2 + p2^2*(2 + 3*q1^2) + 
         10*p2*q1^2*q2 + 2*q2^2 + 3*q1^2*q2^2 + 
         p1*q1*(p2^2 + 40*p2*q2 + q2^2) + p1^2*(-4 + 3*p2^2 + 10*p2*q2 + 
           3*q2^2))*sigm^2 - (-14*p1*p2*q1*q2 + p1^2*(-1 + p2^2 + q2^2) + 
         (-1 + q1^2)*(-1 + p2^2 + q2^2))*sigm^4))))/
 (2*(p1*q1*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
     2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2))^(3/2)*
  (p2*q2*(2*b1^2*p1*q1 + 4*b1*b3*p1*q1*(p2 - q2) + 
     2*b3^2*p1*q1*(p2^2 + q2^2) + sigm^2))^(3/2))
)
}


#Get the mean of sqrt(F12) test statistics
meanT12Func = function(b1, b2, b3, p1, p2, sigm, n) { 
q1=1-p1; q2=1-p2; 
return(
sqrt(2/3)*sqrt(n)*sqrt((b1^2*p1*q1 + 2*b1*b3*p1*q1*(p2 - q2) + b2^2*p2*q2 + 
    2*b2*b3*p2*(p1 - q1)*q2 + b3^2*(p1*q1*(p2 - q2)^2 + p1^2*p2*q2 + 
      p2*q1^2*q2))/sigm^2)
)
}


#Get the variance of sqrt(F12) test statistics
varT12Func = function(b1, b2, b3, p1, p2, sigm) { 
q1=1-p1; q2=1-p2; 
return(
(2*b1^4*p1*q1*(p1^2 + 4*p1*q1 + q1^2) + 8*b1^3*b3*p1*q1*
   (p1^2 + 4*p1*q1 + q1^2)*(p2 - q2) + 2*b2^4*p2*q2*(p2^2 + 4*p2*q2 + q2^2) + 
  8*b2^3*b3*p2*(p1 - q1)*q2*(p2^2 + 4*p2*q2 + q2^2) + 
  4*b1^2*(8*b2^2*p1*p2*q1*q2 - 2*b2*b3*(p1 - q1)*
     (p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2) - 2*p1*p2*q1*q2 + 
      (-1 + p1^2 + 2*p1*q1 + q1^2)*q2^2) + 
    b3^2*(p1^4*(p2^2 + q2^2) + q1^2*(-1 + q1^2)*(p2^2 + q2^2) + 
      p1^2*(-1 + 6*q1^2)*(p2^2 + q2^2) + p1^3*q1*(p2^2 - 4*p2*q2 + q2^2) + 
      p1*q1*(p2^2*(4 + q1^2) - 4*p2*q1^2*q2 + (4 + q1^2)*q2^2)) + 
    2*p1*q1*sigm^2) + 8*b2*b3*(p1 - q1)*
   (b3^2*(-(p2^2*q1^2) + p2^4*q1^2 - 4*p1*p2*q1*(p2 - q2)^2*q2 + 
      p2*(q2 + 2*q1^2*q2) + q1^2*q2^2*(-1 + q2^2) + 
      p1^2*(-p2^2 + p2^4 + 2*p2*q2 - q2^2 + q2^4)) + 2*p2*q2*sigm^2) + 
  4*b2^2*(b3^2*(-4*p1*p2*q1*q2*(p2^2 + q2^2) + 
      p1^2*(p2^4 + p2^3*q2 + q2^2*(-1 + q2^2) + p2*q2*(4 + q2^2) + 
        p2^2*(-1 + 6*q2^2)) + q1^2*(p2^4 + p2^3*q2 + q2^2*(-1 + q2^2) + 
        p2*q2*(4 + q2^2) + p2^2*(-1 + 6*q2^2))) + 2*p2*q2*sigm^2) + 
  b3^2*(b3^2*(-16*p1^3*p2*q1*(p2 - q2)^2*q2 - 8*p1*q1*(p2 - q2)^2*
       (-1 + 2*p2*q1^2*q2) + p1^4*(3*p2^4 - 2*p2^2*q2^2 + 3*q2^4) + 
      q1^2*(3*p2^4*q1^2 + 8*p2*q2 + 3*q2^2*(-1 + q1^2*q2^2) - 
        p2^2*(3 + 2*q1^2*q2^2)) - p1^2*(2*p2^4*q1^2 - 32*p2^3*q1^2*q2 + 
        q2^2*(3 + 2*q1^2*q2^2) + p2^2*(3 + 52*q1^2*q2^2) - 
        8*p2*(q2 + 4*q1^2*q2^3))) + 8*(p1*q1*(p2 - q2)^2 + p1^2*p2*q2 + 
      p2*q1^2*q2)*sigm^2) - 4*b1*b3*(p2 - q2)*
   (2*b2^2*(-2*p1*p2*q1*q2 + p1^2*(-1 + p2^2 + 2*p2*q2 + q2^2) + 
      q1^2*(-1 + p2^2 + 2*p2*q2 + q2^2)) - b2*b3*(p1 - q1)*
     (3 - 2*q1^2 + p2^2*(-2 + q1^2) - 2*p2*q1^2*q2 - 2*q2^2 + q1^2*q2^2 + 
      p1^2*(-2 + p2^2 - 2*p2*q2 + q2^2) - 2*p1*q1*(p2^2 + 8*p2*q2 + q2^2)) - 
    2*(b3^2*(-4*p1^3*p2*q1*q2 + p1^4*(p2^2 + q2^2) + q1^2*(-1 + q1^2)*
         (p2^2 + q2^2) - p1^2*(p2^2 - 8*p2*q1^2*q2 + q2^2) + 
        p1*(q1 + 2*p2^2*q1 - 4*p2*q1^3*q2 + 2*q1*q2^2)) + 2*p1*q1*sigm^2)))/
 (24*(b1^2*p1*q1 + 2*b1*b3*p1*q1*(p2 - q2) + b2^2*p2*q2 + 
   2*b2*b3*p2*(p1 - q1)*q2 + b3^2*(p1*q1*(p2 - q2)^2 + p1^2*p2*q2 + 
     p2*q1^2*q2))*sigm^2)
)
}



#Get the covariance b/w sqrt(F12) and T1 test statistics
covT12T1Func = function(b1, b2, b3, p1, p2, sigm) { 
q1=1-p1; q2=1-p2; 
return(
(p1*q1*(b1^2*b3*p1*q1*(p2 - q2)*
    (2*b2^2*(12*p1*p2*q1*q2 + p1^2*(-2 + 2*p2^2 + 7*p2*q2 + 2*q2^2) + 
       q1^2*(-2 + 2*p2^2 + 7*p2*q2 + 2*q2^2)) - 2*b2*b3*(p1 - q1)*
      (3 - 2*q1^2 + p2^2*(-6 + 5*q1^2) - 8*p2*q1^2*q2 - 6*q2^2 + 
       5*q1^2*q2^2 + 6*p1*q1*(p2^2 - 6*p2*q2 + q2^2) + 
       p1^2*(-2 + 5*p2^2 - 8*p2*q2 + 5*q2^2)) + 
     b3^2*(6*p1^4*p2*q2 + 6*p2*q1^4*q2 - 4*p1^3*q1*(p2^2 - 7*p2*q2 + q2^2) - 
       4*p1^2*q1^2*(2*p2^2 + p2*q2 + 2*q2^2) - 
       4*p1*q1*(1 + p2^2*(-2 + q1^2) - 7*p2*q1^2*q2 + (-2 + q1^2)*q2^2)) + 
     3*(p1^2 + 4*p1*q1 + q1^2)*sigm^2) + 
   b1^3*p1*q1*(4*b2*b3*p2*(p1 - q1)*(p1 + q1)^2*q2 + 
     2*b2^2*p2*(p1^2 + q1^2)*q2 + (p1^2 + 4*p1*q1 + q1^2)*
      (2*b3^2*p2*(p1 - q1)^2*q2 + sigm^2)) - 
   b1*(2*b2^4*p1*p2*q1*q2*(p2^2 - 8*p2*q2 + q2^2) + 
     8*b2^3*b3*p2*(p1 - q1)*q2*(p2^2*(-1 + p1^2 + 3*p1*q1 + q1^2) - 
       3*p1*p2*q1*q2 + (-1 + p1^2 + 3*p1*q1 + q1^2)*q2^2) - 
     b3^4*(4*p1^6*p2*q2*(p2^2 + q2^2) + 4*p2*q1^4*(-1 + q1^2)*q2*
        (p2^2 + q2^2) + p1^5*q1*(p2^4 - 10*p2^3*q2 + 2*p2^2*q2^2 - 
         10*p2*q2^3 + q2^4) + p1^4*(-6*p2^4*q1^2 + 4*p2^3*(-1 + 6*q1^2)*q2 + 
         4*p2^2*q1^2*q2^2 + 4*p2*(-1 + 6*q1^2)*q2^3 - 6*q1^2*q2^4) - 
       2*p1^2*q1^2*(p2^4*(-8 + 3*q1^2) - 4*p2^3*(-5 + 3*q1^2)*q2 - 
         2*p2^2*(-2 + (8 + q1^2)*q2^2) + q2^2*(4 + (-8 + 3*q1^2)*q2^2) - 
         4*p2*q2*(2 + (-5 + 3*q1^2)*q2^2)) + 
       p1*q1^3*(p2^4*(-4 + q1^2) + 2*p2^3*(12 - 5*q1^2)*q2 + 
         q2^2*(3 + (-4 + q1^2)*q2^2) + p2^2*(3 + 2*(-4 + q1^2)*q2^2) - 
         2*p2*q2*(4 + (-12 + 5*q1^2)*q2^2)) + 
       p1^3*q1*(-2*p2^4*(2 + 7*q1^2) + 4*p2^3*(6 + 7*q1^2)*q2 + 3*q2^2 - 
         2*(2 + 7*q1^2)*q2^4 + 4*p2*q2*(-2 + (6 + 7*q1^2)*q2^2) + 
         p2^2*(3 - 4*(2 + 15*q1^2)*q2^2))) + 
     b3^2*(-2*p1^4*(p2^2 + q2^2) - 2*q1^2*(-1 + q1^2)*(p2^2 + q2^2) + 
       p1^3*q1*(p2^2 + 4*p2*q2 + q2^2) + 2*p1^2*(p2^2 - 8*p2*q1^2*q2 + 
         q2^2) + p1*q1*(p2^2*(-8 + q1^2) + 4*p2*q1^2*q2 + (-8 + q1^2)*q2^2))*
      sigm^2 - 4*p1*q1*sigm^4 + 
     2*b2^2*(b3^2*(6*p1^4*p2*q2*(p2^2 + q2^2) + 6*p2*q1^2*(-1 + q1^2)*q2*
          (p2^2 + q2^2) - 2*p1^2*p2*q2*(p2^2*(3 + 14*q1^2) - 4*p2*q1^2*q2 + 
           (3 + 14*q1^2)*q2^2) + p1^3*q1*(2*p2^4 + 3*p2^3*q2 + 
           2*q2^2*(-1 + q2^2) - 2*p2^2*(1 + 2*q2^2) + p2*q2*(8 + 3*q2^2)) + 
         p1*q1*(2*p2^4*q1^2 + p2^3*(8 + 3*q1^2)*q2 + 2*q1^2*q2^2*
            (-1 + q2^2) - 2*p2^2*q1^2*(1 + 2*q2^2) + 
           p2*(8*q2^3 + q1^2*q2*(8 + 3*q2^2)))) - 6*p1*p2*q1*q2*sigm^2) + 
     4*b2*b3*(p1 - q1)*(b3^2*p1*q1*(p2^4*(-2 + 4*q1^2) + 
         p2^3*(-4 + q1^2)*q2 - 2*q1^2*q2^2 - 2*q2^4 + 4*q1^2*q2^4 + 
         2*p1*q1*(2*p2^4 - 5*p2^3*q2 + 2*p2^2*q2^2 - 5*p2*q2^3 + 2*q2^4) + 
         p2*q2*(2 - 4*q2^2 + q1^2*(4 + q2^2)) + 
         p1^2*(4*p2^4 + p2^3*q2 - 2*q2^2 + 4*q2^4 + p2*q2*(4 + q2^2) + 
           p2^2*(-2 + 6*q2^2)) + p2^2*(-4*q2^2 + q1^2*(-2 + 6*q2^2))) + 
       (p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2) - p1*p2*q1*q2 + 
         (-1 + p1^2 + 2*p1*q1 + q1^2)*q2^2)*sigm^2)) - 
   b3*(p2 - q2)*(2*b2^4*p2*q2*(2*p1^2*(-1 + p2^2 + 2*p2*q2 + q2^2) + 
       2*q1^2*(-1 + p2^2 + 2*p2*q2 + q2^2) + p1*q1*(p2^2 + 4*p2*q2 + q2^2)) + 
     2*b2^3*b3*p2*(p1 - q1)*q2*(-3 - 2*q1^2 + p2^2*(2 + 3*q1^2) + 
       10*p2*q1^2*q2 + 2*q2^2 + 3*q1^2*q2^2 + 
       6*p1*q1*(p2^2 + 6*p2*q2 + q2^2) + p1^2*(-2 + 3*p2^2 + 10*p2*q2 + 
         3*q2^2)) + b3^4*(-4*p1^6*p2*q2*(p2^2 + q2^2) - 
       4*p2*q1^4*(-1 + q1^2)*q2*(p2^2 + q2^2) - 
       p1^5*q1*(p2^4 - 10*p2^3*q2 - 6*p2^2*q2^2 - 10*p2*q2^3 + q2^4) + 
       2*p1^4*(p2^4*q1^2 + 2*p2^3*q2 - 14*p2^2*q1^2*q2^2 + 2*p2*q2^3 + 
         q1^2*q2^4) + p1*q1^3*(-(p2^4*(-4 + q1^2)) + 2*p2^3*(-8 + 5*q1^2)*
          q2 - q2^2*(3 + (-4 + q1^2)*q2^2) + 
         2*p2*q2*(2 + (-8 + 5*q1^2)*q2^2) + p2^2*(-3 + (8 + 6*q1^2)*q2^2)) + 
       2*p1^2*q1^2*(p2^4*(-4 + q1^2) + 12*p2^3*q2 + 4*p2*q2*(-1 + 3*q2^2) + 
         q2^2*(2 + (-4 + q1^2)*q2^2) - 2*p2^2*(-1 + (4 + 7*q1^2)*q2^2)) + 
       p1^3*q1*(p2^4*(4 + 6*q1^2) - 4*p2^3*(4 + 3*q1^2)*q2 - 
         4*p2*q2*(-1 + (4 + 3*q1^2)*q2^2) + q2^2*(-3 + (4 + 6*q1^2)*q2^2) + 
         p2^2*(-3 + (8 + 60*q1^2)*q2^2))) + 
     b3^2*(-2*p1^4*(p2^2 + q2^2) - 2*q1^2*(-1 + q1^2)*(p2^2 + q2^2) + 
       p1^3*q1*(p2^2 + 8*p2*q2 + q2^2) + p1*q1*(-2 + p2^2*(-4 + q1^2) + 
         8*p2*q1^2*q2 + (-4 + q1^2)*q2^2) + 2*p1^2*(p2^2*(1 + 2*q1^2) - 
         8*p2*q1^2*q2 + (1 + 2*q1^2)*q2^2))*sigm^2 - 4*p1*q1*sigm^4 + 
     2*b2^2*(b3^2*(-2*p1^4*p2*q2*(-1 + p2^2 - 4*p2*q2 + q2^2) + 
         2*p2*q1^2*q2*(-(p2^2*(-3 + q1^2)) + 4*p2*q1^2*q2 - 
           (-3 + q1^2)*(-1 + q2^2)) - 2*p1^2*p2*q2*(3 + p2^2*(-3 + 6*q1^2) + 
           40*p2*q1^2*q2 - 3*q2^2 + q1^2*(-2 + 6*q2^2)) + 
         p1^3*q1*(4*p2^4 + 11*p2^3*q2 + 4*q2^2*(-1 + q2^2) + 
           p2*q2*(4 + 11*q2^2) + p2^2*(-4 + 48*q2^2)) + 
         p1*q1*(4*p2^4*q1^2 + p2^3*(-12 + 11*q1^2)*q2 + 
           4*q1^2*q2^2*(-1 + q2^2) + 4*p2^2*q1^2*(-1 + 12*q2^2) + 
           p2*q2*(10 - 12*q2^2 + q1^2*(4 + 11*q2^2)))) + 
       (p1^2 + q1^2)*(-1 + p2^2 + 2*p2*q2 + q2^2)*sigm^2) - 
     b2*b3*(p1 - q1)*(2*b3^2*(p1^4*p2*q2*(-2 + 5*p2^2 - 2*p2*q2 + 5*q2^2) + 
         p2*q1^2*q2*(3 + p2^2*(-6 + 5*q1^2) - 2*p2*q1^2*q2 - 6*q2^2 + 
           q1^2*(-2 + 5*q2^2)) - p1^3*q1*(3*p2^4 + 8*p2^3*q2 + 
           q2^2*(-2 + 3*q2^2) + p2^2*(-2 + 30*q2^2) + 4*p2*(q2 + 2*q2^3)) - 
         p1^2*(2*p2^4*q1^2 + 6*p2^3*(1 + q1^2)*q2 - 32*p2^2*q1^2*q2^2 + 
           2*q1^2*q2^4 + p2*q2*(-3 + 6*q2^2 + q1^2*(4 + 6*q2^2))) - 
         p1*q1*(p2^4*(2 + 3*q1^2) + 4*p2^3*(-3 + 2*q1^2)*q2 + 
           q2^2*(-3 + 2*q2^2 + q1^2*(-2 + 3*q2^2)) + 
           2*p2*q2*(3 - 6*q2^2 + q1^2*(2 + 4*q2^2)) + 
           p2^2*(-3 + 4*q2^2 + q1^2*(-2 + 30*q2^2)))) + 
       (3 - 2*q1^2 + p2^2*(-2 + q1^2) - 2*p2*q1^2*q2 - 2*q2^2 + q1^2*q2^2 + 
         p1^2*(-2 + p2^2 - 2*p2*q2 + q2^2) - 2*p1*q1*(p2^2 + 10*p2*q2 + 
           q2^2))*sigm^2))))/
 (4*sqrt(3)*sqrt((b1^2*p1*q1 + 2*b1*b3*p1*q1*(p2 - q2) + b2^2*p2*q2 + 
     2*b2*b3*p2*(p1 - q1)*q2 + b3^2*(p1*q1*(p2 - q2)^2 + p1^2*p2*q2 + 
       p2*q1^2*q2))/sigm^2)*sigm^2*
  (p1*q1*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 2*b3^2*p2*(p1^2 + q1^2)*
      q2 + sigm^2))^(3/2))
)
}


#Get the mean of F3|1 T1 test statistics
meanF3l1Func = function(b1, b2, b3, p1, p2, p3, sigm) { 
q1=1-p1; q2=1-p2; q3=1-p3;
return(
(0.125*(8*b2^2*p2*q2 + 
   (b3*(b2*(p1 - q1)*(p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2)*(p3 - q3)^2 + 
        (-1 + p1^2 + 2*p1*q1 + q1^2)*q2^2*(p3 - q3)^2 - 
        2*p2*q2*(p3^2*(-1 + p1^2 + 2*p1*q1 + q1^2) - 4*p1*p3*q1*q3 + 
          (-1 + p1^2 + 2*p1*q1 + q1^2)*q3^2)) + 
      2*b3*p1*q1*(p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2)*(p3 - q3)^2 + 
        (-1 + p1^2 + 2*p1*q1 + q1^2)*q2^2*(p3 - q3)^2 - 
        2*p2*q2*(p3^2*(-1 + p1^2 + 2*p1*q1 + q1^2) - 2*p3*(p1 + q1)^2*q3 + 
          (-1 + p1^2 + 2*p1*q1 + q1^2)*q3^2))))/(p1*p3*q1*q3) + 
   (8*b2^2*p1*p2*p3*q1*q2*q3 + b2*b3*(p1 - q1)*
      (p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2)*(p3 - q3)^2 + 
       (-1 + p1^2 + 2*p1*q1 + q1^2)*q2^2*(p3 - q3)^2 - 
       2*p2*q2*(p3^2*(-1 + p1^2 + 2*p1*q1 + q1^2) - 4*p1*p3*q1*q3 + 
         (-1 + p1^2 + 2*p1*q1 + q1^2)*q3^2)) + 
     b3^2*(-(p1^4*(p2^2*(p3 - q3)^2 + q2^2*(p3 - q3)^2 - 
          2*p2*q2*(p3^2 + q3^2))) - 2*p1*q1*(p2^2*(p3 - q3)^2 + 
         q2^2*(p3 - q3)^2 - 2*p2*q2*(p3^2 + q3^2)) - 
       q1^2*(-1 + q1^2)*(p2^2*(p3 - q3)^2 + q2^2*(p3 - q3)^2 - 
         2*p2*q2*(p3^2 + q3^2)) + p1^2*(p2^2*(1 + 2*q1^2)*(p3 - q3)^2 + 
         (1 + 2*q1^2)*q2^2*(p3 - q3)^2 - 2*p2*q2*(p3^2*(1 + 2*q1^2) - 
           8*p3*q1^2*q3 + (1 + 2*q1^2)*q3^2))))/(p1*p3*q1*q3) + 8*sigm^2))/
 (2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 2*b3^2*p2*(p1^2 + q1^2)*q2 + 
  sigm^2)
)
}




#Get the var of F3|1 T1 test statistics
varF3l1Func = function(b1, b2, b3, p1, p2, p3, sigm) { 
q1=1-p1; q2=1-p2; q3=1-p3;
return(
(0.03125*(64*b2^4*p2^2*q2^2 + 
   (b3^2*(b2*(p1 - q1)*(p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2)*(p3 - q3)^2 + 
         (-1 + p1^2 + 2*p1*q1 + q1^2)*q2^2*(p3 - q3)^2 - 
         2*p2*q2*(p3^2*(-1 + p1^2 + 2*p1*q1 + q1^2) - 4*p1*p3*q1*q3 + 
           (-1 + p1^2 + 2*p1*q1 + q1^2)*q3^2)) + 
       2*b3*p1*q1*(p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2)*(p3 - q3)^2 + 
         (-1 + p1^2 + 2*p1*q1 + q1^2)*q2^2*(p3 - q3)^2 - 
         2*p2*q2*(p3^2*(-1 + p1^2 + 2*p1*q1 + q1^2) - 2*p3*(p1 + q1)^2*q3 + 
           (-1 + p1^2 + 2*p1*q1 + q1^2)*q3^2)))^2)/(p1^2*p3^2*q1^2*q3^2) + 
   (48*b2^2*b3^2*p2*q2*(-(p1^4*(p2^2*(p3 - q3)^2 + q2^2*(p3 - q3)^2 - 
         2*p2*q2*(p3^2 + q3^2))) - 2*p1*q1*(p2^2*(p3 - q3)^2 + 
        q2^2*(p3 - q3)^2 - 2*p2*q2*(p3^2 + q3^2)) - 
      q1^2*(-1 + q1^2)*(p2^2*(p3 - q3)^2 + q2^2*(p3 - q3)^2 - 
        2*p2*q2*(p3^2 + q3^2)) + p1^2*(p2^2*(1 + 2*q1^2)*(p3 - q3)^2 + 
        (1 + 2*q1^2)*q2^2*(p3 - q3)^2 - 2*p2*q2*(p3^2*(1 + 2*q1^2) - 
          8*p3*q1^2*q3 + (1 + 2*q1^2)*q3^2))))/(p1*p3*q1*q3) + 
   (8*b2^2*p1*p2*p3*q1*q2*q3 + b2*b3*(p1 - q1)*
       (p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2)*(p3 - q3)^2 + 
        (-1 + p1^2 + 2*p1*q1 + q1^2)*q2^2*(p3 - q3)^2 - 
        2*p2*q2*(p3^2*(-1 + p1^2 + 2*p1*q1 + q1^2) - 4*p1*p3*q1*q3 + 
          (-1 + p1^2 + 2*p1*q1 + q1^2)*q3^2)) + 
      b3^2*(-(p1^4*(p2^2*(p3 - q3)^2 + q2^2*(p3 - q3)^2 - 
           2*p2*q2*(p3^2 + q3^2))) - 2*p1*q1*(p2^2*(p3 - q3)^2 + 
          q2^2*(p3 - q3)^2 - 2*p2*q2*(p3^2 + q3^2)) - q1^2*(-1 + q1^2)*
         (p2^2*(p3 - q3)^2 + q2^2*(p3 - q3)^2 - 2*p2*q2*(p3^2 + q3^2)) + 
        p1^2*(p2^2*(1 + 2*q1^2)*(p3 - q3)^2 + (1 + 2*q1^2)*q2^2*(p3 - q3)^2 - 
          2*p2*q2*(p3^2*(1 + 2*q1^2) - 8*p3*q1^2*q3 + (1 + 2*q1^2)*q3^2))))^2/
    (p1^2*p3^2*q1^2*q3^2) + (2*b2*b3*(8*b2*p1*p2*p3*q1*q2*q3 + 
      b3*(p1 - q1)*(p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2)*(p3 - q3)^2 + 
        (-1 + p1^2 + 2*p1*q1 + q1^2)*q2^2*(p3 - q3)^2 - 
        2*p2*q2*(p3^2*(-1 + p1^2 + 2*p1*q1 + q1^2) - 4*p1*p3*q1*q3 + 
          (-1 + p1^2 + 2*p1*q1 + q1^2)*q3^2)))*
     (b2*(p1 - q1)*(p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2)*(p3 - q3)^2 + 
        (-1 + p1^2 + 2*p1*q1 + q1^2)*q2^2*(p3 - q3)^2 - 
        2*p2*q2*(p3^2*(-1 + p1^2 + 2*p1*q1 + q1^2) - 4*p1*p3*q1*q3 + 
          (-1 + p1^2 + 2*p1*q1 + q1^2)*q3^2)) + 
      b3*(2*p1^3*q1*(p2 - q2)^2*(p3 - q3)^2 - 
        p1^4*(p2^2*(p3 - q3)^2 + q2^2*(p3 - q3)^2 - 2*p2*q2*(p3^2 + q3^2)) - 
        q1^2*(-1 + q1^2)*(p2^2*(p3 - q3)^2 + q2^2*(p3 - q3)^2 - 
          2*p2*q2*(p3^2 + q3^2)) + 2*p1*q1*(p2^2*(-2 + q1^2)*(p3 - q3)^2 + 
          (-2 + q1^2)*q2^2*(p3 - q3)^2 - 2*p2*q2*(p3^2*(-2 + q1^2) - 
            2*p3*q1^2*q3 + (-2 + q1^2)*q3^2)) + 
        p1^2*(p2^2*(1 + 6*q1^2)*(p3 - q3)^2 + (1 + 6*q1^2)*q2^2*(p3 - q3)^2 - 
          2*p2*q2*(p3^2*(1 + 6*q1^2) - 16*p3*q1^2*q3 + (1 + 6*q1^2)*q3^2)))))/
    (p1^2*p3^2*q1^2*q3^2) + 64*b2^2*p2*q2*sigm^2 + 
   (8*b2*(8*b2*p1*p2*p3*q1*q2*q3 + b3*(p1 - q1)*
       (p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2)*(p3 - q3)^2 + 
        (-1 + p1^2 + 2*p1*q1 + q1^2)*q2^2*(p3 - q3)^2 - 
        2*p2*q2*(p3^2*(-1 + p1^2 + 2*p1*q1 + q1^2) - 4*p1*p3*q1*q3 + 
          (-1 + p1^2 + 2*p1*q1 + q1^2)*q3^2)))*sigm^2)/(p1*p3*q1*q3) + 
   (8*b3*(b2*(p1 - q1)*(p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2)*(p3 - q3)^2 + 
        (-1 + p1^2 + 2*p1*q1 + q1^2)*q2^2*(p3 - q3)^2 - 
        2*p2*q2*(p3^2*(-1 + p1^2 + 2*p1*q1 + q1^2) - 4*p1*p3*q1*q3 + 
          (-1 + p1^2 + 2*p1*q1 + q1^2)*q3^2)) + 
      2*b3*p1*q1*(p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2)*(p3 - q3)^2 + 
        (-1 + p1^2 + 2*p1*q1 + q1^2)*q2^2*(p3 - q3)^2 - 
        2*p2*q2*(p3^2*(-1 + p1^2 + 2*p1*q1 + q1^2) - 2*p3*(p1 + q1)^2*q3 + 
          (-1 + p1^2 + 2*p1*q1 + q1^2)*q3^2)))*sigm^2)/(p1*p3*q1*q3) + 
   (8*b3^2*(-(p1^4*(p2^2*(p3 - q3)^2 + q2^2*(p3 - q3)^2 - 
         2*p2*q2*(p3^2 + q3^2))) - 2*p1*q1*(p2^2*(p3 - q3)^2 + 
        q2^2*(p3 - q3)^2 - 2*p2*q2*(p3^2 + q3^2)) - 
      q1^2*(-1 + q1^2)*(p2^2*(p3 - q3)^2 + q2^2*(p3 - q3)^2 - 
        2*p2*q2*(p3^2 + q3^2)) + p1^2*(p2^2*(1 + 2*q1^2)*(p3 - q3)^2 + 
        (1 + 2*q1^2)*q2^2*(p3 - q3)^2 - 2*p2*q2*(p3^2*(1 + 2*q1^2) - 
          8*p3*q1^2*q3 + (1 + 2*q1^2)*q3^2)))*sigm^2)/(p1*p3*q1*q3) + 
   32*sigm^4))/(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
   2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2)^2
)
}


#Get the mean of sqrt(F1|3) test statistics
meanT1l3Func = function(b1, b2, b3, p1, p2, sigm, n) { 
q1=1-p1; q2=1-p2; 
return(
sqrt(n)*sqrt((p1*q1*(b1 + b3*(p2 - q2))^2)/(2*b2^2*p2*q2 + 
    4*b2*b3*p2*(p1 - q1)*q2 + 2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2))
)
}



#Get the variance of sqrt(F1|3) test statistics
varT1l3Func = function(b1, b2, b3, p1, p2, sigm) { 
q1=1-p1; q2=1-p2; 
return(
(16*b2*p1*p2*(b2 + b3*(p1 - q1))*q1*(b1 + b3*(p2 - q2))^2*q2*sigm^2 + 
  4*p1*q1*(b1 + b3*(p2 - q2))^2*sigm^4 + 
  4*sigm^2*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
    2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2)*
   (2*b3*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2) + 2*b2^2*p2*q2 + 
    4*b2*b3*p2*(p1 - q1)*q2 + 2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2) + 
  8*b3*p1*q1*(b1 + b3*(p2 - q2))*sigm^2*
   (2*b1*p2*(b2*(p1 - q1) + b3*(p1^2 + q1^2))*q2 - 
    (p2 - q2)*(2*b2^2*p2*q2 + 2*b2*b3*p2*(p1 - q1)*q2 + sigm^2)) - 
  2*(b1 + b3*(p2 - q2))*(2*b2^2*p2*(p1 - q1)*q2 + 4*b2*b3*p2*(p1^2 + q1^2)*
     q2 + (p1 - q1)*(2*b3^2*p2*(p1 + q1)^2*q2 + sigm^2))*
   (-2*b2*b3*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2)^2 - 
    (b1*(p1 - q1) + b2*(p2 - q2))*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
      2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2)) - 
  (4*b2*b3*p2*(p1 - q1)^3*q2 + 2*b3^2*p2*(p1 - q1)^2*(p1^2 + q1^2)*q2 + 
    2*b2^2*p2*(p1^2 - 4*p1*q1 + q1^2)*q2 + (p1^2 - 4*p1*q1 + q1^2)*sigm^2)*
   (-2*b3^2*p1*q1*(b1 + b3*(p2 - q2))^2*(p2 - q2)^2 + 
    (b1^2 - b3^2*(p2 - q2)^2)*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
      2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2)) + 
  4*b2*p2*q2*(2*b2*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2) - 
    (p1 - q1)*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
      2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2))*
   (b1*(b2^2 + 2*b2*b3*(p1 - q1) + b3^2*(p1^2 + q1^2))*(p2 - q2) + 
    b3*(b2^2*(p2 + q2)^2 + 2*b2*b3*(p1 - q1)*(p2 + q2)^2 + 
      b3^2*(p1^2 + q1^2)*(p2 + q2)^2 + 2*sigm^2)) - 
  4*b2^2*p1*p2*q1*(b1 + b3*(p2 - q2))*q2*
   (b1*(b2^2 + 2*b2*b3*(p1 - q1) + b3^2*(p1^2 + q1^2))*
     (p2^2 - 4*p2*q2 + q2^2) + b3*(p2 - q2)*(b2^2*(p2^2 + q2^2) + 
      2*b2*b3*(p1 - q1)*(p2^2 + q2^2) + b3^2*(p1^2 + q1^2)*(p2^2 + q2^2) + 
      2*sigm^2)) - (2*(4*b2*b3*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2) + 
     (b2 + b3*(-p1 + q1))*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
       2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2))*
    (b3^3*(p1 - q1)*(p1^3*q1*(p2^2 - q2^2)^2 + 2*p1^4*p2*q2*(p2^2 + q2^2) + 
       2*p2*q1^2*(-1 + q1^2)*q2*(p2^2 + q2^2) + p1*q1*(p2 - q2)^2*
        (1 + p2^2*(-2 + q1^2) + 2*p2*q1^2*q2 + (-2 + q1^2)*q2^2) + 
       2*p1^2*(p2^4*q1^2 + 6*p2^2*q1^2*q2^2 - p2*(1 + 2*q1^2)*q2^3 + 
         q1^2*q2^4 - p2^3*(q2 + 2*q1^2*q2))) + 
     2*b2*b3^2*(2*p1^4*p2*q2*(p2^2 + q2^2) + 2*p2*q1^2*(-1 + q1^2)*q2*
        (p2^2 + q2^2) - 2*p1^2*p2*q2*(p2^2*(1 + 2*q1^2) - 8*p2*q1^2*q2 + 
         (1 + 2*q1^2)*q2^2) - p1^3*q1*(p2^4 + 2*p2*q2 + q2^2*(-1 + q2^2) + 
         p2^2*(-1 + 10*q2^2)) + p1*(-(p2^4*q1^3) + 4*p2^3*q1*q2 + 
         p2^2*q1^3*(1 - 10*q2^2) - q1^3*q2^2*(-1 + q2^2) + 
         p2*(-2*q1^3*q2 + 4*q1*q2^3))) - 4*b2*p1*p2*q1*q2*
      (2*b2^2*p2*q2 + sigm^2) + b1*p1*q1*(p2 - q2)*
      (-2*b2*b3*(p1^2 + q1^2)*(-1 + p2^2 + 2*p2*q2 + q2^2) + 
       b3^2*(p1 - q1)*(1 + p2^2*(-2 + p1^2 + 2*p1*q1 + q1^2) - 
         2*p2*(p1 + q1)^2*q2 + (-2 + p1^2 + 2*p1*q1 + q1^2)*q2^2) + 
       (p1 - q1)*sigm^2) + b3*(p1 - q1)*
      (2*b2^2*p2*q2*(p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2) - 12*p1*p2*q1*q2 + 
         (-1 + p1^2 + 2*p1*q1 + q1^2)*q2^2) + 
       (p2^2*(-1 + p1^2 + 3*p1*q1 + q1^2) - 6*p1*p2*q1*q2 + 
         (-1 + p1^2 + 3*p1*q1 + q1^2)*q2^2)*sigm^2)))/(p1*q1) + 
  (2*b3*(2*b3*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2) + 2*b2^2*p2*q2 + 
     4*b2*b3*p2*(p1 - q1)*q2 + 2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2)*
    (-2*b2^2*p1*p2*q1*(p1^2 + q1^2)*(b1 + b3*(p2 - q2))^2*(p2 - q2)*q2 + 
     2*b2*b3*p1*(p1 - q1)*q1*(b1 + b3*(p2 - q2))^2*(p2 - q2)*
      (-1 + p2^2*q1^2 + q1^2*q2^2 + p1^2*(p2^2 + q2^2)) + 
     b3^2*p1*q1*(p1^2 + q1^2)*(b1 + b3*(p2 - q2))^2*(p2 - q2)*
      (-1 + p2^2*q1^2 + q1^2*q2^2 + p1^2*(p2^2 + q2^2)) + 
     b3*(p1^2 + q1^2)*(b1 + b3*(p2 - q2))*(p2^2 - (p1^2 + q1^2)*(p2 - q2)^2 + 
       q2^2)*(2*b3*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2) + 2*b2^2*p2*q2 + 
       4*b2*b3*p2*(p1 - q1)*q2 + 2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2) + 
     2*b2*p2*(p1^2 + q1^2)*(b1 + b3*(p2 - q2))*q2*
      (2*b2*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2) - 
       (p1 - q1)*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
         2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2)) + 
     (p1 - q1)*(b1 + b3*(p2 - q2))*(p2^2 - (p1^2 + q1^2)*(p2 - q2)^2 + q2^2)*
      (4*b2*b3*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2) + (b2 + b3*(-p1 + q1))*
        (2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 2*b3^2*p2*(p1^2 + q1^2)*
          q2 + sigm^2)) + 2*p1*(p1 - q1)*q1*(b1 + b3*(p2 - q2))*(p2 - q2)*
      (-2*b2*b3*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2)^2 - 
       (b1*(p1 - q1) + b2*(p2 - q2))*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*
          q2 + 2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2)) + 
     p1*q1*(p1^2 + q1^2)*(p2 - q2)*(-2*b3^2*p1*q1*(b1 + b3*(p2 - q2))^2*
        (p2 - q2)^2 + (b1^2 - b3^2*(p2 - q2)^2)*(2*b2^2*p2*q2 + 
         4*b2*b3*p2*(p1 - q1)*q2 + 2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2))))/
   (p1*q1*(b1 + b3*(p2 - q2))) - 
  4*b2*b3*(-2*b2^2*p1*p2*(p1 - q1)*q1*(b1 + b3*(p2 - q2))^2*q2*
     (p2^2 + q2^2) + b3^2*p1*(p1 - q1)*q1*(b1 + b3*(p2 - q2))^2*(p2^2 + q2^2)*
     (-1 + p2^2*q1^2 + q1^2*q2^2 + p1^2*(p2^2 + q2^2)) + 
    2*b2*b3*p1*q1*(b1 + b3*(p2 - q2))^2*(p2^2 + q2^2)*
     (p1^2*(-1 + p2^2 + q2^2) + q1^2*(-1 + p2^2 + q2^2) - 
      2*p1*q1*(p2^2 + q2^2)) + b3*(p1 - q1)*(b1 + b3*(p2 - q2))*(p2 - q2)*
     (1 - (p1^2 + q1^2)*(p2^2 + q2^2))*(2*b3*p1*q1*(b1 + b3*(p2 - q2))*
       (p2 - q2) + 2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
      2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2) + 2*b2*p2*(p1 - q1)*
     (b1 + b3*(p2 - q2))*(p2 - q2)*q2*(2*b2*p1*q1*(b1 + b3*(p2 - q2))*
       (p2 - q2) - (p1 - q1)*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
        2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2)) + (b1 + b3*(p2 - q2))*(p2 - q2)*
     (p1^2 + q1^2 - (p1 - q1)^2*(p2^2 + q2^2))*
     (4*b2*b3*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2) + (b2 + b3*(-p1 + q1))*
       (2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 2*b3^2*p2*(p1^2 + q1^2)*q2 + 
        sigm^2)) + 2*p1*q1*(b1 + b3*(p2 - q2))*(p2^2 + q2^2)*
     (-2*b2*b3*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2)^2 - 
      (b1*(p1 - q1) + b2*(p2 - q2))*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
        2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2)) + p1*(p1 - q1)*q1*(p2^2 + q2^2)*
     (-2*b3^2*p1*q1*(b1 + b3*(p2 - q2))^2*(p2 - q2)^2 + 
      (b1^2 - b3^2*(p2 - q2)^2)*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
        2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2))) - 
  2*b3^2*(-2*b2^2*p1*p2*q1*(p1^2 + q1^2)*(b1 + b3*(p2 - q2))^2*q2*
     (p2^2 + q2^2) + 2*b2*b3*p1*(p1 - q1)*q1*(b1 + b3*(p2 - q2))^2*
     (p2^2 + q2^2)*(-1 + p2^2*q1^2 + q1^2*q2^2 + p1^2*(p2^2 + q2^2)) + 
    b3^2*p1*q1*(p1^2 + q1^2)*(b1 + b3*(p2 - q2))^2*(p2^2 + q2^2)*
     (-1 + p2^2*q1^2 + q1^2*q2^2 + p1^2*(p2^2 + q2^2)) + 
    b3*(p1^2 + q1^2)*(b1 + b3*(p2 - q2))*(p2 - q2)*
     (1 - (p1^2 + q1^2)*(p2^2 + q2^2))*(2*b3*p1*q1*(b1 + b3*(p2 - q2))*
       (p2 - q2) + 2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
      2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2) + 2*b2*p2*(p1^2 + q1^2)*
     (b1 + b3*(p2 - q2))*(p2 - q2)*q2*(2*b2*p1*q1*(b1 + b3*(p2 - q2))*
       (p2 - q2) - (p1 - q1)*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
        2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2)) + (p1 - q1)*(b1 + b3*(p2 - q2))*
     (p2 - q2)*(1 - (p1^2 + q1^2)*(p2^2 + q2^2))*
     (4*b2*b3*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2) + (b2 + b3*(-p1 + q1))*
       (2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 2*b3^2*p2*(p1^2 + q1^2)*q2 + 
        sigm^2)) + 2*p1*(p1 - q1)*q1*(b1 + b3*(p2 - q2))*(p2^2 + q2^2)*
     (-2*b2*b3*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2)^2 - 
      (b1*(p1 - q1) + b2*(p2 - q2))*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
        2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2)) + p1*q1*(p1^2 + q1^2)*
     (p2^2 + q2^2)*(-2*b3^2*p1*q1*(b1 + b3*(p2 - q2))^2*(p2 - q2)^2 + 
      (b1^2 - b3^2*(p2 - q2)^2)*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
        2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2))))/
 (8*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 2*b3^2*p2*(p1^2 + q1^2)*q2 + 
    sigm^2)^3)
)
}



#Get the variance b/w sqrt(F1|3) and sqrt(F2|3) test statistics
covT1l3T2l3Func = function(b1, b2, b3, p1, p2, sigm) { 
q1=1-p1; q2=1-p2; 
return(
(2*p1*p2*(b2 + b3*(p1 - q1))^2*q1*(b1 + b3*(p2 - q2))^2*q2*sigm^4 - 
  4*b1*p1*p2*(b2 + b3*(p1 - q1))^2*q1*(b1 + b3*(p2 - q2))*q2*sigm^2*
   (2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 2*b3^2*p2*(p1^2 + q1^2)*q2 + 
    sigm^2) - 4*p1*p2*(b2 + b3*(p1 - q1))^2*q1*(b1 + b3*(p2 - q2))^2*q2*
   sigm^2*(2*b1^2*p1*q1 + 4*b1*b3*p1*q1*(p2 - q2) + 2*b2*b3*p2*(p1 - q1)*q2 + 
    2*b3^2*(p1*q1*(p2 - q2)^2 + p1^2*p2*q2 + p2*q1^2*q2) + sigm^2) + 
  b1^2*p1*p2*(b2 + b3*(p1 - q1))^2*q1*(b1 + b3*(p2 - q2))^2*q2*
   (4*b2*b3*p2*(p1 - q1)^3*q2 + 2*b3^2*p2*(p1 - q1)^2*(p1^2 + q1^2)*q2 + 
    2*b2^2*p2*(p1^2 - 4*p1*q1 + q1^2)*q2 + (p1^2 - 4*p1*q1 + q1^2)*sigm^2) + 
  4*b3*p1*p2*(b2 + b3*(p1 - q1))^2*q1*(b1 + b3*(p2 - q2))*q2*sigm^2*
   (2*b1*p2*(b2*(p1 - q1) + b3*(p1^2 + q1^2))*q2 - 
    (p2 - q2)*(2*b2^2*p2*q2 + 2*b2*b3*p2*(p1 - q1)*q2 + sigm^2)) + 
  b1*p1*(b2 + b3*(p1 - q1))*q1*(b1 + b3*(p2 - q2))^2*
   (2*b2^2*p2*(p1 - q1)*q2 + 4*b2*b3*p2*(p1^2 + q1^2)*q2 + 
    (p1 - q1)*(2*b3^2*p2*(p1 + q1)^2*q2 + sigm^2))*
   (2*b1^2*p1*q1*(p2 - q2) - 2*b1*(b2*p2*(p1 - q1)*q2 + 
      b3*(p1^2*p2*q2 + p2*q1^2*q2 - 2*p1*q1*(p2^2 - p2*q2 + q2^2))) + 
    (p2 - q2)*(2*b3^2*p1*q1*(p2^2 + q2^2) + sigm^2)) + 
  2*p1*p2*(b2 + b3*(p1 - q1))*q1*(b1 + b3*(p2 - q2))*q2*
   (-2*b1*b3*p2*(b2 + b3*(p1 - q1))*(p1 - q1)^2*q2 - 
    (b1*(p1 - q1) + b2*(p2 - q2))*(2*b1^2*p1*q1 + 4*b1*b3*p1*q1*(p2 - q2) + 
      2*b3^2*p1*q1*(p2^2 + q2^2) + sigm^2))*
   (b1*(b2^2 + 2*b2*b3*(p1 - q1) + b3^2*(p1^2 + q1^2))*(p2 - q2) + 
    b3*(b2^2*(p2 + q2)^2 + 2*b2*b3*(p1 - q1)*(p2 + q2)^2 + 
      b3^2*(p1^2 + q1^2)*(p2 + q2)^2 + 2*sigm^2)) + 
  p1*p2*q1*(b1 + b3*(p2 - q2))*q2*(-2*b3^2*p2*(b2 + b3*(p1 - q1))^2*
     (p1 - q1)^2*q2 + (b2^2 - b3^2*(p1 - q1)^2)*(2*b1^2*p1*q1 + 
      4*b1*b3*p1*q1*(p2 - q2) + 2*b3^2*p1*q1*(p2^2 + q2^2) + sigm^2))*
   (b1*(b2^2 + 2*b2*b3*(p1 - q1) + b3^2*(p1^2 + q1^2))*
     (p2^2 - 4*p2*q2 + q2^2) + b3*(p2 - q2)*(b2^2*(p2^2 + q2^2) + 
      2*b2*b3*(p1 - q1)*(p2^2 + q2^2) + b3^2*(p1^2 + q1^2)*(p2^2 + q2^2) + 
      2*sigm^2)) - (b2 + b3*(p1 - q1))*(b1 + b3*(p2 - q2))*
   (4*b1*b3*p2*(b2 + b3*(p1 - q1))*(p1 - q1)*q2 + 
    (b1 + b3*(-p2 + q2))*(2*b1^2*p1*q1 + 4*b1*b3*p1*q1*(p2 - q2) + 
      2*b3^2*p1*q1*(p2^2 + q2^2) + sigm^2))*
   (b3^3*(p1 - q1)*(p1^3*q1*(p2^2 - q2^2)^2 + 2*p1^4*p2*q2*(p2^2 + q2^2) + 
      2*p2*q1^2*(-1 + q1^2)*q2*(p2^2 + q2^2) + p1*q1*(p2 - q2)^2*
       (1 + p2^2*(-2 + q1^2) + 2*p2*q1^2*q2 + (-2 + q1^2)*q2^2) + 
      2*p1^2*(p2^4*q1^2 + 6*p2^2*q1^2*q2^2 - p2*(1 + 2*q1^2)*q2^3 + 
        q1^2*q2^4 - p2^3*(q2 + 2*q1^2*q2))) + 
    2*b2*b3^2*(2*p1^4*p2*q2*(p2^2 + q2^2) + 2*p2*q1^2*(-1 + q1^2)*q2*
       (p2^2 + q2^2) - 2*p1^2*p2*q2*(p2^2*(1 + 2*q1^2) - 8*p2*q1^2*q2 + 
        (1 + 2*q1^2)*q2^2) - p1^3*q1*(p2^4 + 2*p2*q2 + q2^2*(-1 + q2^2) + 
        p2^2*(-1 + 10*q2^2)) + p1*(-(p2^4*q1^3) + 4*p2^3*q1*q2 + 
        p2^2*q1^3*(1 - 10*q2^2) - q1^3*q2^2*(-1 + q2^2) + 
        p2*(-2*q1^3*q2 + 4*q1*q2^3))) - 4*b2*p1*p2*q1*q2*
     (2*b2^2*p2*q2 + sigm^2) + b1*p1*q1*(p2 - q2)*
     (-2*b2*b3*(p1^2 + q1^2)*(-1 + p2^2 + 2*p2*q2 + q2^2) + 
      b3^2*(p1 - q1)*(1 + p2^2*(-2 + p1^2 + 2*p1*q1 + q1^2) - 
        2*p2*(p1 + q1)^2*q2 + (-2 + p1^2 + 2*p1*q1 + q1^2)*q2^2) + 
      (p1 - q1)*sigm^2) + b3*(p1 - q1)*
     (2*b2^2*p2*q2*(p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2) - 12*p1*p2*q1*q2 + 
        (-1 + p1^2 + 2*p1*q1 + q1^2)*q2^2) + 
      (p2^2*(-1 + p1^2 + 3*p1*q1 + q1^2) - 6*p1*p2*q1*q2 + 
        (-1 + p1^2 + 3*p1*q1 + q1^2)*q2^2)*sigm^2)) - 
  2*b1*b3*p2*(b2 + b3*(p1 - q1))^2*q2*
   (-2*b2^2*p1*p2*q1*(p1^2 + q1^2)*(b1 + b3*(p2 - q2))^2*(p2 - q2)*q2 + 
    2*b2*b3*p1*(p1 - q1)*q1*(b1 + b3*(p2 - q2))^2*(p2 - q2)*
     (-1 + p2^2*q1^2 + q1^2*q2^2 + p1^2*(p2^2 + q2^2)) + 
    b3^2*p1*q1*(p1^2 + q1^2)*(b1 + b3*(p2 - q2))^2*(p2 - q2)*
     (-1 + p2^2*q1^2 + q1^2*q2^2 + p1^2*(p2^2 + q2^2)) + 
    b3*(p1^2 + q1^2)*(b1 + b3*(p2 - q2))*(p2^2 - (p1^2 + q1^2)*(p2 - q2)^2 + 
      q2^2)*(2*b3*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2) + 2*b2^2*p2*q2 + 
      4*b2*b3*p2*(p1 - q1)*q2 + 2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2) + 
    2*b2*p2*(p1^2 + q1^2)*(b1 + b3*(p2 - q2))*q2*
     (2*b2*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2) - 
      (p1 - q1)*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
        2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2)) + (p1 - q1)*(b1 + b3*(p2 - q2))*
     (p2^2 - (p1^2 + q1^2)*(p2 - q2)^2 + q2^2)*
     (4*b2*b3*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2) + (b2 + b3*(-p1 + q1))*
       (2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 2*b3^2*p2*(p1^2 + q1^2)*q2 + 
        sigm^2)) + 2*p1*(p1 - q1)*q1*(b1 + b3*(p2 - q2))*(p2 - q2)*
     (-2*b2*b3*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2)^2 - 
      (b1*(p1 - q1) + b2*(p2 - q2))*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
        2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2)) + p1*q1*(p1^2 + q1^2)*(p2 - q2)*
     (-2*b3^2*p1*q1*(b1 + b3*(p2 - q2))^2*(p2 - q2)^2 + 
      (b1^2 - b3^2*(p2 - q2)^2)*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
        2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2))) + 
  b3*(b2 + b3*(p1 - q1))*(2*b1^2*p1*q1 + 4*b1*b3*p1*q1*(p2 - q2) + 
    2*b3*p2*(b2 + b3*(p1 - q1))*(p1 - q1)*q2 + 2*b3^2*p1*q1*(p2^2 + q2^2) + 
    sigm^2)*(-2*b2^2*p1*p2*(p1 - q1)*q1*(b1 + b3*(p2 - q2))^2*q2*
     (p2^2 + q2^2) + b3^2*p1*(p1 - q1)*q1*(b1 + b3*(p2 - q2))^2*(p2^2 + q2^2)*
     (-1 + p2^2*q1^2 + q1^2*q2^2 + p1^2*(p2^2 + q2^2)) + 
    2*b2*b3*p1*q1*(b1 + b3*(p2 - q2))^2*(p2^2 + q2^2)*
     (p1^2*(-1 + p2^2 + q2^2) + q1^2*(-1 + p2^2 + q2^2) - 
      2*p1*q1*(p2^2 + q2^2)) + b3*(p1 - q1)*(b1 + b3*(p2 - q2))*(p2 - q2)*
     (1 - (p1^2 + q1^2)*(p2^2 + q2^2))*(2*b3*p1*q1*(b1 + b3*(p2 - q2))*
       (p2 - q2) + 2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
      2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2) + 2*b2*p2*(p1 - q1)*
     (b1 + b3*(p2 - q2))*(p2 - q2)*q2*(2*b2*p1*q1*(b1 + b3*(p2 - q2))*
       (p2 - q2) - (p1 - q1)*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
        2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2)) + (b1 + b3*(p2 - q2))*(p2 - q2)*
     (p1^2 + q1^2 - (p1 - q1)^2*(p2^2 + q2^2))*
     (4*b2*b3*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2) + (b2 + b3*(-p1 + q1))*
       (2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 2*b3^2*p2*(p1^2 + q1^2)*q2 + 
        sigm^2)) + 2*p1*q1*(b1 + b3*(p2 - q2))*(p2^2 + q2^2)*
     (-2*b2*b3*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2)^2 - 
      (b1*(p1 - q1) + b2*(p2 - q2))*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
        2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2)) + p1*(p1 - q1)*q1*(p2^2 + q2^2)*
     (-2*b3^2*p1*q1*(b1 + b3*(p2 - q2))^2*(p2 - q2)^2 + 
      (b1^2 - b3^2*(p2 - q2)^2)*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
        2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2))) - 
  b3^2*p2*(b2 + b3*(p1 - q1))^2*q2*(-2*b2^2*p1*p2*q1*(p1^2 + q1^2)*
     (b1 + b3*(p2 - q2))^2*q2*(p2^2 + q2^2) + 2*b2*b3*p1*(p1 - q1)*q1*
     (b1 + b3*(p2 - q2))^2*(p2^2 + q2^2)*(-1 + p2^2*q1^2 + q1^2*q2^2 + 
      p1^2*(p2^2 + q2^2)) + b3^2*p1*q1*(p1^2 + q1^2)*(b1 + b3*(p2 - q2))^2*
     (p2^2 + q2^2)*(-1 + p2^2*q1^2 + q1^2*q2^2 + p1^2*(p2^2 + q2^2)) + 
    b3*(p1^2 + q1^2)*(b1 + b3*(p2 - q2))*(p2 - q2)*
     (1 - (p1^2 + q1^2)*(p2^2 + q2^2))*(2*b3*p1*q1*(b1 + b3*(p2 - q2))*
       (p2 - q2) + 2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
      2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2) + 2*b2*p2*(p1^2 + q1^2)*
     (b1 + b3*(p2 - q2))*(p2 - q2)*q2*(2*b2*p1*q1*(b1 + b3*(p2 - q2))*
       (p2 - q2) - (p1 - q1)*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
        2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2)) + (p1 - q1)*(b1 + b3*(p2 - q2))*
     (p2 - q2)*(1 - (p1^2 + q1^2)*(p2^2 + q2^2))*
     (4*b2*b3*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2) + (b2 + b3*(-p1 + q1))*
       (2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 2*b3^2*p2*(p1^2 + q1^2)*q2 + 
        sigm^2)) + 2*p1*(p1 - q1)*q1*(b1 + b3*(p2 - q2))*(p2^2 + q2^2)*
     (-2*b2*b3*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2)^2 - 
      (b1*(p1 - q1) + b2*(p2 - q2))*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
        2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2)) + p1*q1*(p1^2 + q1^2)*
     (p2^2 + q2^2)*(-2*b3^2*p1*q1*(b1 + b3*(p2 - q2))^2*(p2 - q2)^2 + 
      (b1^2 - b3^2*(p2 - q2)^2)*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
        2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2))))/
 (4*sqrt((p1*q1*(b1 + b3*(p2 - q2))^2)/(2*b2^2*p2*q2 + 
     4*b2*b3*p2*(p1 - q1)*q2 + 2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2))*
  (2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 2*b3^2*p2*(p1^2 + q1^2)*q2 + 
    sigm^2)^2*sqrt((p2*(b2 + b3*(p1 - q1))^2*q2)/
    (2*b1^2*p1*q1 + 4*b1*b3*p1*q1*(p2 - q2) + 2*b3^2*p1*q1*(p2^2 + q2^2) + 
     sigm^2))*(2*b1^2*p1*q1 + 4*b1*b3*p1*q1*(p2 - q2) + 
    2*b3^2*p1*q1*(p2^2 + q2^2) + sigm^2)^2)
)
}


#Get the covariance b/w T1 and sqrt(F1|3) test statistics
covT1T1l3Func = function(b1, b2, b3, p1, p2, sigm) { 
q1=1-p1; q2=1-p2; 
return(
((b1 + b3*(p2 - q2))*(16*b2*p1^2*p2*(b2 + b3*(p1 - q1))*q1^2*
    (b1 + b3*(p2 - q2))^2*q2*sigm^2 + 4*p1^2*q1^2*(b1 + b3*(p2 - q2))^2*
    sigm^4 + 4*p1*q1*sigm^2*(2*p2*((b2 + b3*p1)^2 - 2*b2*b3*q1 + b3^2*q1^2)*
      q2 + sigm^2)*(2*(b3^2*p1*p2^2*q1 + b1*b3*p1*q1*(p2 - q2) + 
       p2*(b2 + b3*(p1 - q1))^2*q2 + b3^2*p1*q1*q2^2) + sigm^2) - 
   4*b2*p1*p2*q1*q2*(((b2 + b3*p1)^2 - 2*b2*b3*q1 + b3^2*q1^2)*
      (b1*(p2 - q2) + b3*(p2 + q2)^2) + 2*b3*sigm^2)*
    (-2*b2*b3*p1*p2^2*q1 + 2*p2*(b2 + b3*(p1 - q1))*
      (b2*(p1 - q1) + b3*(p1^2 + q1^2))*q2 - 2*b2*b3*p1*q1*q2^2 + 
     2*b1*b2*p1*q1*(-p2 + q2) + (p1 - q1)*sigm^2) - 
   4*b2^2*p1^2*p2*q1^2*(b1 + b3*(p2 - q2))*q2*
    (((b2 + b3*p1)^2 - 2*b2*b3*q1 + b3^2*q1^2)*(b3*(p2 - q2)*(p2^2 + q2^2) + 
       b1*(p2^2 - 4*p2*q2 + q2^2)) + 2*b3*(p2 - q2)*sigm^2) + 
   8*b3*p1^2*q1^2*(b1 + b3*(p2 - q2))*sigm^2*
    (2*b1*p2*(b2*(p1 - q1) + b3*(p1^2 + q1^2))*q2 - 
     (p2 - q2)*(2*b2*p2*(b2 + b3*(p1 - q1))*q2 + sigm^2)) - 
   p1*q1*(b1 + b3*(p2 - q2))*(2*p2*(b2 + b3*(p1 - q1))*
      (b3*(p1 - q1)*(p1^2 + q1^2) + b2*(p1^2 - 4*p1*q1 + q1^2))*q2 + 
     (p1^2 - 4*p1*q1 + q1^2)*sigm^2)*
    (b1*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
       2*b3^2*(-(p1*p2^2*q1) + p2*(p1 + q1)^2*q2 - p1*q1*q2^2) + sigm^2) - 
     b3*(p2 - q2)*(2*(b3^2*p1*p2^2*q1 + p2*(b2 + b3*(p1 - q1))^2*q2 + 
         b3^2*p1*q1*q2^2) + sigm^2)) + 
   2*(4*b2*b3*p1*q1*(b1 + b3*(p2 - q2))*(p2 - q2) + 
     (b2 + b3*(-p1 + q1))*(2*p2*((b2 + b3*p1)^2 - 2*b2*b3*q1 + b3^2*q1^2)*
        q2 + sigm^2))*(2*b2*b3^2*(p1*p2^2*(-1 + p2^2)*q1*(p1^2 + q1^2) - 
       2*p2*(p1^4*p2^2 - p1^3*q1 + p2^2*q1^2*(-1 + q1^2) - 
         p1^2*p2^2*(1 + 2*q1^2) + p1*(2*p2^2*q1 - q1^3))*q2 + 
       p1*q1*(p1^2*(-1 + 10*p2^2) - 16*p1*p2^2*q1 + (-1 + 10*p2^2)*q1^2)*
        q2^2 - 2*p2*(p1 - q1)^2*(-1 + p1 + q1)*(1 + p1 + q1)*q2^3 + 
       p1*q1*(p1^2 + q1^2)*q2^4) - b3^3*(p1 - q1)*
      (p1*p2^2*q1*(1 + p2^2*(-2 + (p1 + q1)^2)) + 
       2*p2*(p1^4*p2^2 + p1*(-1 + 2*p2^2)*q1 + p2^2*q1^2*(-1 + q1^2) - 
         p1^2*p2^2*(1 + 2*q1^2))*q2 - 
       p1*q1*(-1 + 2*p2^2*(2 + p1^2 - 6*p1*q1 + q1^2))*q2^2 + 
       2*p2*(p1 - q1)^2*(-1 + p1 + q1)*(1 + p1 + q1)*q2^3 + 
       p1*q1*(-2 + (p1 + q1)^2)*q2^4) + 4*b2*p1*p2*q1*q2*
      (2*b2^2*p2*q2 + sigm^2) - b1*p1*q1*(p2 - q2)*
      (-2*b2*b3*(p1^2 + q1^2)*(-1 + p2 + q2)*(1 + p2 + q2) + 
       b3^2*(p1 - q1)*(1 + p2^2*(-2 + (p1 + q1)^2) - 2*p2*(p1 + q1)^2*q2 + 
         (-2 + (p1 + q1)^2)*q2^2) + (p1 - q1)*sigm^2) - 
     b3*(p1 - q1)*(2*b2^2*p2*q2*(p2^2*(-1 + p1 + q1)*(1 + p1 + q1) - 
         12*p1*p2*q1*q2 + (-1 + p1 + q1)*(1 + p1 + q1)*q2^2) + 
       (p2^2*(-1 + p1^2 + 3*p1*q1 + q1^2) - 6*p1*p2*q1*q2 + 
         (-1 + p1^2 + 3*p1*q1 + q1^2)*q2^2)*sigm^2)) + 
   2*b3*(2*(b3^2*p1*p2^2*q1 + b1*b3*p1*q1*(p2 - q2) + 
       p2*(b2 + b3*(p1 - q1))^2*q2 + b3^2*p1*q1*q2^2) + sigm^2)*
    (-2*b2^3*p2*(p1 - q1)*q2*(p2^2*(-1 + p1 + q1)*(1 + p1 + q1) - 
       4*p1*p2*q1*q2 + (-1 + p1 + q1)*(1 + p1 + q1)*q2^2) - 
     4*b2^2*b3*p2*q2*(p2^2*(p1^4 + p1^3*q1 - q1^2 + q1^4 - 
         p1^2*(1 + 2*q1^2) + p1*(q1 + q1^3)) - 
       2*p1*p2*q1*(3*p1^2 - 4*p1*q1 + 3*q1^2)*q2 + 
       (p1^4 + p1^3*q1 - q1^2 + q1^4 - p1^2*(1 + 2*q1^2) + p1*(q1 + q1^3))*
        q2^2) - b1*p1*q1*(p2 - q2)*(-8*b2^2*p1*p2*q1*q2 + 
       2*b2*b3*(p1 - q1)*(1 + p2^2*(-2 + (p1 + q1)^2) - 
         2*p2*(p1^2 + 6*p1*q1 + q1^2)*q2 + (-2 + (p1 + q1)^2)*q2^2) + 
       b3^2*(p1^2 + q1^2)*(1 + p2^2*(-2 + (p1 + q1)^2) - 
         2*p2*(p1^2 + 6*p1*q1 + q1^2)*q2 + (-2 + (p1 + q1)^2)*q2^2) + 
       (p1^2 - 4*p1*q1 + q1^2)*sigm^2) - b2*(p1 - q1)*
      (2*b3^2*(p1*p2^2*q1*(1 + p2^2*(-2 + (p1 + q1)^2)) + 
         p2*(p1^2*(-1 + p1^2)*p2^2 + 2*p1*(-1 + p1^2*p2^2)*q1 - 
           (1 + 6*p1^2)*p2^2*q1^2 + 2*p1*p2^2*q1^3 + p2^2*q1^4)*q2 - 
         p1*(-1 + 2*p2^2*(2 + 3*(p1 - q1)^2))*q1*q2^2 + 
         p2*(p1^4 + 2*p1^3*q1 - q1^2 + 2*p1*q1^3 + q1^4 - p1^2*(1 + 6*q1^2))*
          q2^3 + p1*q1*(-2 + (p1 + q1)^2)*q2^4) + 
       (p2^2*(-1 + p1 + q1)*(1 + p1 + q1) - 4*p1*p2*q1*q2 + 
         (-1 + p1 + q1)*(1 + p1 + q1)*q2^2)*sigm^2) - 
     b3*p1*q1*(b3^2*(p1^2 + q1^2)*(p2^4*(-2 + (p1 + q1)^2) - 
         8*p1*p2^3*q1*q2 + q2^2 + (-2 + (p1 + q1)^2)*q2^4 + 
         p2^2*(1 - 2*(2 + p1^2 - 6*p1*q1 + q1^2)*q2^2) - 
         2*p2*(q2 + 4*p1*q1*q2^3)) + (p2^2*(-2 + 3*p1^2 + 3*q1^2) - 
         6*p2*(p1^2 + q1^2)*q2 + (-2 + 3*p1^2 + 3*q1^2)*q2^2)*sigm^2)) + 
   2*p1*q1*(b1 + b3*(p2 - q2))*(2*p2*(b2 + b3*(p1 - q1))*
      (b2*(p1 - q1) + b3*(p1 + q1)^2)*q2 + (p1 - q1)*sigm^2)*
    (b2*(p2 - q2)*(2*(b3^2*p1*p2^2*q1 + p2*(b2 + b3*(p1 - q1))^2*q2 + 
         b3^2*p1*q1*q2^2) + sigm^2) + b1*(2*b2^2*p2*(p1 - q1)*q2 + 
       2*b2*b3*(p1*p2^2*q1 + 2*p2*(p1^2 - 3*p1*q1 + q1^2)*q2 + p1*q1*q2^2) + 
       (p1 - q1)*(2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2))) - 
   4*b2*b3*p1*q1*(b1 + b3*(p2 - q2))*
    (b1*p1*q1*(-2*b2*b3*(p1^2 + q1^2)*(p2^4 - 2*p2^3*q2 - q2^2 + q2^4 - 
         2*p2*q2*(-2 + q2^2) + p2^2*(-1 + 2*q2^2)) + 
       b3^2*(p1 - q1)*(-(p2^4*(p1 + q1)^2) + 2*p2^3*(p1 + q1)^2*q2 + q2^2 - 
         (p1 + q1)^2*q2^4 + p2^2*(1 - 2*(p1 + q1)^2*q2^2) + 
         2*p2*q2*(-2 + (p1 + q1)^2*q2^2)) - (p1 - q1)*(8*b2^2*p2^2*q2^2 + 
         (p2^2 + q2^2)*sigm^2)) - (p2 - q2)*
      (b3^3*(p1 - q1)*(p1*p2^2*q1*(-1 + p2*(p1 + q1))*(1 + p2*(p1 + q1)) + 
         2*p2*(p1^4 + p1^3*p2^2*q1 - q1^2 + q1^4 + p1*q1*(2 + p2^2*q1^2) + 
           p1^2*(-1 - 2*(-1 + p2^2)*q1^2))*q2 + 
         p1*q1*(-1 + 2*p2^2*(p1 + q1)^2)*q2^2 + 2*p1*p2*(p1 - q1)^2*q1*q2^3 + 
         p1*q1*(p1 + q1)^2*q2^4) + 2*b2*b3^2*(p1*p2^2*(-1 + p2^2)*q1*
          (p1^2 + q1^2) + p2*(p1^4*(1 + p2^2) + 2*p1^3*p2^2*q1 - 2*q1^2 + 
           (1 + p2^2)*q1^4 + 2*p1*q1*(2 + p2^2*q1^2) + 
           p1^2*(-2 + (2 - 6*p2^2)*q1^2))*q2 + (p1^2 + q1^2)*
          (2*p1^2*p2^2 + 2*p2^2*q1^2 - p1*(q1 + 2*p2^2*q1))*q2^2 + 
         p2*(p1 - q1)^2*(p1^2 + 4*p1*q1 + q1^2)*q2^3 + p1*q1*(p1^2 + q1^2)*
          q2^4) + b2*(-4*p1*p2*q1*q2 + p1^2*(-1 + p2 + q2)*(1 + p2 + q2) + 
         q1^2*(-1 + p2 + q2)*(1 + p2 + q2))*(2*b2^2*p2*q2 + sigm^2) + 
       b3*(p1 - q1)*(2*b2^2*p2*q2*(-1 + 2*p1*q1*(p2 - q2)^2 + 
           p1^2*(-1 + 2*(p2 + q2)^2) + q1^2*(-1 + 2*(p2 + q2)^2)) + 
         (-1 + p1^2 + q1^2 + 3*p1*q1*(p2^2 + q2^2))*sigm^2))) + 
   2*b3^2*p1*q1*(b1 + b3*(p2 - q2))*
    (b1*p1*q1*(8*b2^2*p2*q2*(-(p2*q1) + p1*q2)*(p1*p2 - q1*q2) + 
       2*b2*b3*(p1 - q1)*(p2^4*(p1 + q1)^2 + 4*p2*q2 - 
         2*p2^3*(p1^2 + 6*p1*q1 + q1^2)*q2 - q2^2 - 
         2*p2*(p1^2 + 6*p1*q1 + q1^2)*q2^3 + (p1 + q1)^2*q2^4 + 
         p2^2*(-1 + 2*(p1 + q1)^2*q2^2)) + b3^2*(p1^2 + q1^2)*
        (p2^4*(p1 + q1)^2 + 4*p2*q2 - 2*p2^3*(p1^2 + 6*p1*q1 + q1^2)*q2 - 
         q2^2 - 2*p2*(p1^2 + 6*p1*q1 + q1^2)*q2^3 + (p1 + q1)^2*q2^4 + 
         p2^2*(-1 + 2*(p1 + q1)^2*q2^2)) + (p1^2 - 4*p1*q1 + q1^2)*
        (p2^2 + q2^2)*sigm^2) + (p2 - q2)*
      (2*b2^3*p2*(p1 - q1)*q2*(-1 + p1^2*(p2 + q2)^2 + q1^2*(p2 + q2)^2 + 
         2*p1*q1*(p2^2 + q2^2)) + 4*b2^2*b3*p2*q2*
        (-(p1^2*(1 + 2*q1^2*(p2 - q2)^2)) + p1*(q1 + q1^3*(p2 - q2)^2) + 
         p1^3*q1*(p2 - q2)^2 + p1^4*(p2 + q2)^2 + q1^2*(-1 + q1*(p2 + q2))*
          (1 + q1*(p2 + q2))) + b2*(p1 - q1)*
        (2*b3^2*(p1^4*p2*q2*(p2 + q2)^2 + p1^3*q1*(p2^2 + q2^2)*
            (p2^2 + 4*p2*q2 + q2^2) + p2*q1^2*q2*(-1 + q1*(p2 + q2))*
            (1 + q1*(p2 + q2)) + p1*q1*(p2^2 + q2^2)*
            (-1 + q1^2*(p2^2 + 4*p2*q2 + q2^2)) + 
           p1^2*(2*p2^4*q1^2 - 2*p2^3*q1^2*q2 + 8*p2^2*q1^2*q2^2 + 
             2*q1^2*q2^4 - p2*(q2 + 2*q1^2*q2^3))) + 
         (-1 + p2^2*(p1 + q1)^2 + 2*p2*(p1^2 + q1^2)*q2 + (p1 + q1)^2*q2^2)*
          sigm^2) + b3*p1*q1*(-2*sigm^2 + (p1^2 + q1^2)*(p2^2 + q2^2)*
          (b3^2*(-1 + p2^2*(p1 + q1)^2 + 2*p2*(p1 - q1)^2*q2 + 
             (p1 + q1)^2*q2^2) + 3*sigm^2))))))/
 (4*sqrt(2)*sqrt((p1*q1*(b1 + b3*(p2 - q2))^2)/
    (2*p2*((b2 + b3*p1)^2 - 2*b2*b3*q1 + b3^2*q1^2)*q2 + sigm^2))*
  (2*p2*((b2 + b3*p1)^2 - 2*b2*b3*q1 + b3^2*q1^2)*q2 + sigm^2)^3*
  sqrt(p1*q1*(2*p2*((b2 + b3*p1)^2 - 2*b2*b3*q1 + b3^2*q1^2)*q2 + sigm^2)))

)
}


#Get the covariance b/w T1 and sqrt(F2|3) test statistics
covT1T2l3Func = function(b1, b2, b3, p1, p2, sigm) { 
q1=1-p1; q2=1-p2; 
return(
(p1*q1*(2*p1*p2*(b2 + b3*(p1 - q1))^2*q1*(b1 + b3*(p2 - q2))*q2*sigm^4 - 
   4*b1*p1*p2*(b2 + b3*(p1 - q1))^2*q1*q2*sigm^2*
    (2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 2*b3^2*p2*(p1^2 + q1^2)*q2 + 
     sigm^2) - 4*p1*p2*(b2 + b3*(p1 - q1))^2*q1*(b1 + b3*(p2 - q2))*q2*sigm^2*
    (2*b1^2*p1*q1 + 4*b1*b3*p1*q1*(p2 - q2) + 2*b2*b3*p2*(p1 - q1)*q2 + 
     2*b3^2*(p1*q1*(p2 - q2)^2 + p1^2*p2*q2 + p2*q1^2*q2) + sigm^2) + 
   b1^2*p1*p2*(b2 + b3*(p1 - q1))^2*q1*(b1 + b3*(p2 - q2))*q2*
    (4*b2*b3*p2*(p1 - q1)^3*q2 + 2*b3^2*p2*(p1 - q1)^2*(p1^2 + q1^2)*q2 + 
     2*b2^2*p2*(p1^2 - 4*p1*q1 + q1^2)*q2 + (p1^2 - 4*p1*q1 + q1^2)*sigm^2) + 
   4*b3*p1*p2*(b2 + b3*(p1 - q1))^2*q1*q2*sigm^2*
    (2*b1*p2*(b2*(p1 - q1) + b3*(p1^2 + q1^2))*q2 - 
     (p2 - q2)*(2*b2^2*p2*q2 + 2*b2*b3*p2*(p1 - q1)*q2 + sigm^2)) - 
   b1*p1*(b2 + b3*(p1 - q1))*q1*(b1 + b3*(p2 - q2))*
    (2*b2^2*p2*(p1 - q1)*q2 + 4*b2*b3*p2*(p1^2 + q1^2)*q2 + 
     (p1 - q1)*(2*b3^2*p2*(p1 + q1)^2*q2 + sigm^2))*
    (2*b1*p2*(b2 + b3*(p1 - q1))*(p1 - q1)*q2 - 
     (p2 - q2)*(2*b1^2*p1*q1 + 4*b1*b3*p1*q1*(p2 - q2) + 
       2*b3^2*p1*q1*(p2^2 + q2^2) + sigm^2)) + 2*p1*p2*(b2 + b3*(p1 - q1))*q1*
    q2*(-2*b1*b3*p2*(b2 + b3*(p1 - q1))*(p1 - q1)^2*q2 - 
     (b1*(p1 - q1) + b2*(p2 - q2))*(2*b1^2*p1*q1 + 4*b1*b3*p1*q1*(p2 - q2) + 
       2*b3^2*p1*q1*(p2^2 + q2^2) + sigm^2))*
    (b1*(b2^2 + 2*b2*b3*(p1 - q1) + b3^2*(p1^2 + q1^2))*(p2 - q2) + 
     b3*(b2^2*(p2 + q2)^2 + 2*b2*b3*(p1 - q1)*(p2 + q2)^2 + 
       b3^2*(p1^2 + q1^2)*(p2 + q2)^2 + 2*sigm^2)) + 
   p1*p2*q1*q2*(-2*b3^2*p2*(b2 + b3*(p1 - q1))^2*(p1 - q1)^2*q2 + 
     (b2^2 - b3^2*(p1 - q1)^2)*(2*b1^2*p1*q1 + 4*b1*b3*p1*q1*(p2 - q2) + 
       2*b3^2*p1*q1*(p2^2 + q2^2) + sigm^2))*
    (b1*(b2^2 + 2*b2*b3*(p1 - q1) + b3^2*(p1^2 + q1^2))*
      (p2^2 - 4*p2*q2 + q2^2) + b3*(p2 - q2)*(b2^2*(p2^2 + q2^2) + 
       2*b2*b3*(p1 - q1)*(p2^2 + q2^2) + b3^2*(p1^2 + q1^2)*(p2^2 + q2^2) + 
       2*sigm^2)) + (b2 + b3*(p1 - q1))*
    (4*b1*b3*p2*(b2 + b3*(p1 - q1))*(p1 - q1)*q2 + (b1 + b3*(-p2 + q2))*
      (2*b1^2*p1*q1 + 4*b1*b3*p1*q1*(p2 - q2) + 2*b3^2*p1*q1*(p2^2 + q2^2) + 
       sigm^2))*(-(b3^3*(p1 - q1)*(p1^3*q1*(p2^2 - q2^2)^2 + 
        2*p1^4*p2*q2*(p2^2 + q2^2) + 2*p2*q1^2*(-1 + q1^2)*q2*(p2^2 + q2^2) + 
        p1*q1*(p2 - q2)^2*(1 + p2^2*(-2 + q1^2) + 2*p2*q1^2*q2 + 
          (-2 + q1^2)*q2^2) + 2*p1^2*(p2^4*q1^2 + 6*p2^2*q1^2*q2^2 - 
          p2*(1 + 2*q1^2)*q2^3 + q1^2*q2^4 - p2^3*(q2 + 2*q1^2*q2)))) + 
     2*b2*b3^2*(-2*p1^4*p2*q2*(p2^2 + q2^2) - 2*p2*q1^2*(-1 + q1^2)*q2*
        (p2^2 + q2^2) + 2*p1^2*p2*q2*(p2^2*(1 + 2*q1^2) - 8*p2*q1^2*q2 + 
         (1 + 2*q1^2)*q2^2) + p1^3*q1*(p2^4 + 2*p2*q2 + q2^2*(-1 + q2^2) + 
         p2^2*(-1 + 10*q2^2)) + p1*q1*(p2^4*q1^2 - 4*p2^3*q2 + 
         2*p2*q2*(q1^2 - 2*q2^2) + q1^2*q2^2*(-1 + q2^2) + 
         p2^2*q1^2*(-1 + 10*q2^2))) + 4*b2*p1*p2*q1*q2*
      (2*b2^2*p2*q2 + sigm^2) - b1*p1*q1*(p2 - q2)*
      (-2*b2*b3*(p1^2 + q1^2)*(-1 + p2^2 + 2*p2*q2 + q2^2) + 
       b3^2*(p1 - q1)*(1 + p2^2*(-2 + p1^2 + 2*p1*q1 + q1^2) - 
         2*p2*(p1 + q1)^2*q2 + (-2 + p1^2 + 2*p1*q1 + q1^2)*q2^2) + 
       (p1 - q1)*sigm^2) - b3*(p1 - q1)*
      (2*b2^2*p2*q2*(p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2) - 12*p1*p2*q1*q2 + 
         (-1 + p1^2 + 2*p1*q1 + q1^2)*q2^2) + 
       (p2^2*(-1 + p1^2 + 3*p1*q1 + q1^2) - 6*p1*p2*q1*q2 + 
         (-1 + p1^2 + 3*p1*q1 + q1^2)*q2^2)*sigm^2)) - 
   2*b1*b3*p2*(b2 + b3*(p1 - q1))^2*q2*
    (-2*b2^3*p2*(p1 - q1)*q2*(p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2) - 
       4*p1*p2*q1*q2 + (-1 + p1^2 + 2*p1*q1 + q1^2)*q2^2) - 
     4*b2^2*b3*p2*q2*(p1^4*(p2^2 + q2^2) + q1^2*(-1 + q1^2)*(p2^2 + q2^2) + 
       p1^3*q1*(p2^2 - 6*p2*q2 + q2^2) + p1*q1*(p2^2*(1 + q1^2) - 
         6*p2*q1^2*q2 + (1 + q1^2)*q2^2) - p1^2*(p2^2*(1 + 2*q1^2) - 
         8*p2*q1^2*q2 + (1 + 2*q1^2)*q2^2)) - b1*p1*q1*(p2 - q2)*
      (-8*b2^2*p1*p2*q1*q2 + 2*b2*b3*(p1 - q1)*
        (1 + p2^2*(-2 + p1^2 + 2*p1*q1 + q1^2) - 2*p2*(p1^2 + 6*p1*q1 + q1^2)*
          q2 + (-2 + p1^2 + 2*p1*q1 + q1^2)*q2^2) + b3^2*(p1^2 + q1^2)*
        (1 + p2^2*(-2 + p1^2 + 2*p1*q1 + q1^2) - 2*p2*(p1^2 + 6*p1*q1 + q1^2)*
          q2 + (-2 + p1^2 + 2*p1*q1 + q1^2)*q2^2) + (p1^2 - 4*p1*q1 + q1^2)*
        sigm^2) - b2*(p1 - q1)*(2*b3^2*(p1^4*p2*q2*(p2^2 + q2^2) + 
         p2*q1^2*(-1 + q1^2)*q2*(p2^2 + q2^2) + p1^3*q1*(p2 - q2)^2*
          (p2^2 + 4*p2*q2 + q2^2) + p1^2*(2*p2^4*q1^2 + 12*p2^2*q1^2*q2^2 - 
           p2*(1 + 6*q1^2)*q2^3 + 2*q1^2*q2^4 - p2^3*(q2 + 6*q1^2*q2)) + 
         p1*q1*(p2^4*(-2 + q1^2) + 2*p2^3*q1^2*q2 + q2^2 + (-2 + q1^2)*q2^4 + 
           2*p2*q2*(-1 + q1^2*q2^2) + p2^2*(1 - 2*(2 + 3*q1^2)*q2^2))) + 
       (p2^2*(-1 + p1^2 + 2*p1*q1 + q1^2) - 4*p1*p2*q1*q2 + 
         (-1 + p1^2 + 2*p1*q1 + q1^2)*q2^2)*sigm^2) - 
     b3*p1*q1*(b3^2*(2*p1^3*q1*(p2 - q2)^4 + 2*p1*q1^3*(p2 - q2)^4 + 
         p1^4*(p2^2 - q2^2)^2 + p1^2*(2*p2^4*(-1 + q1^2) - 2*p2*q2 + q2^2 + 
           2*(-1 + q1^2)*q2^4 + p2^2*(1 - 4*(1 + q1^2)*q2^2)) + 
         q1^2*(p2^4*(-2 + q1^2) - 2*p2*q2 + q2^2 + (-2 + q1^2)*q2^4 + 
           p2^2*(1 - 2*(2 + q1^2)*q2^2))) + (p2^2*(-2 + 3*p1^2 + 3*q1^2) - 
         6*p2*(p1^2 + q1^2)*q2 + (-2 + 3*p1^2 + 3*q1^2)*q2^2)*sigm^2)) + 
   b3^2*p2*(b2 + b3*(p1 - q1))^2*q2*
    (b1*p1*q1*(8*b2^2*p2*q2*(p1^2*p2*q2 + p2*q1^2*q2 - p1*q1*(p2^2 + q2^2)) + 
       b3^2*(p1^4*(p2 - q2)^2*(p2^2 + q2^2) + 2*p1^3*q1*(p2^4 - 6*p2^3*q2 + 
           2*p2^2*q2^2 - 6*p2*q2^3 + q2^4) + 2*p1*q1^3*(p2^4 - 6*p2^3*q2 + 
           2*p2^2*q2^2 - 6*p2*q2^3 + q2^4) + q1^2*(-p2^2 + p2^4*q1^2 + 
           4*p2*q2 - 2*p2^3*q1^2*q2 - q2^2 + 2*p2^2*q1^2*q2^2 - 
           2*p2*q1^2*q2^3 + q1^2*q2^4) + p1^2*(-p2^2 + 2*p2^4*q1^2 + 
           4*p2*q2 - 4*p2^3*q1^2*q2 - q2^2 + 4*p2^2*q1^2*q2^2 - 
           4*p2*q1^2*q2^3 + 2*q1^2*q2^4)) + 
       2*b2*b3*(p1^3*(p2 - q2)^2*(p2^2 + q2^2) + 
         p1^2*q1*(p2^4 - 10*p2^3*q2 + 2*p2^2*q2^2 - 10*p2*q2^3 + q2^4) + 
         q1*(-(p2^4*q1^2) + 2*p2^3*q1^2*q2 + q2^2 - q1^2*q2^4 + 
           p2^2*(1 - 2*q1^2*q2^2) + 2*p2*q2*(-2 + q1^2*q2^2)) - 
         p1*(p2^4*q1^2 - 10*p2^3*q1^2*q2 + q2^2 + q1^2*q2^4 + 
           p2^2*(1 + 2*q1^2*q2^2) - 2*p2*q2*(2 + 5*q1^2*q2^2))) + 
       (p1^2 - 4*p1*q1 + q1^2)*(p2^2 + q2^2)*sigm^2) + 
     (p2 - q2)*(2*b2^3*p2*(p1 - q1)*q2*(-1 + p2^2*q1^2 + 2*p2*q1^2*q2 + 
         q1^2*q2^2 + p1^2*(p2 + q2)^2 + 2*p1*q1*(p2^2 + q2^2)) + 
       4*b2^2*b3*p2*q2*(p1*(q1 + q1^3*(p2 - q2)^2) + p1^3*q1*(p2 - q2)^2 + 
         p1^4*(p2 + q2)^2 + q1^2*(-1 + p2^2*q1^2 + 2*p2*q1^2*q2 + 
           q1^2*q2^2) - p1^2*(1 + 2*p2^2*q1^2 - 4*p2*q1^2*q2 + 
           2*q1^2*q2^2)) + b3*p1*q1*(b3^2*(p1^2 + q1^2)*(p2^2 + q2^2)*
          (-1 + p2^2*q1^2 + 2*p1*q1*(p2 - q2)^2 + 2*p2*q1^2*q2 + q1^2*q2^2 + 
           p1^2*(p2 + q2)^2) + (-2 + 3*p2^2*q1^2 + 3*q1^2*q2^2 + 
           3*p1^2*(p2^2 + q2^2))*sigm^2) + b2*(p1 - q1)*
        (2*b3^2*(p1^4*p2*q2*(p2 + q2)^2 + p2*q1^2*q2*(-1 + p2^2*q1^2 + 
             2*p2*q1^2*q2 + q1^2*q2^2) + p1*q1*(p2^2 + q2^2)*
            (-1 + p2^2*q1^2 + 4*p2*q1^2*q2 + q1^2*q2^2) + 
           p1^3*q1*(p2^4 + 4*p2^3*q2 + 2*p2^2*q2^2 + 4*p2*q2^3 + q2^4) + 
           p1^2*(2*p2^4*q1^2 - 2*p2^3*q1^2*q2 + 8*p2^2*q1^2*q2^2 + 
             2*q1^2*q2^4 - p2*(q2 + 2*q1^2*q2^3))) + 
         (-1 + p2^2*q1^2 + 2*p2*q1^2*q2 + q1^2*q2^2 + p1^2*(p2 + q2)^2 + 
           2*p1*q1*(p2^2 + q2^2))*sigm^2))) + b3*(b2 + b3*(p1 - q1))*
    (2*b1^2*p1*q1 + 4*b1*b3*p1*q1*(p2 - q2) + 2*b3*p2*(b2 + b3*(p1 - q1))*
      (p1 - q1)*q2 + 2*b3^2*p1*q1*(p2^2 + q2^2) + sigm^2)*
    (b1*p1*q1*(-2*b2*b3*(p1^2 + q1^2)*(p2^4 - 2*p2^3*q2 - 
         2*p2*q2*(-2 + q2^2) + q2^2*(-1 + q2^2) + p2^2*(-1 + 2*q2^2)) + 
       b3^2*(-(p1^3*(p2 - q2)^2*(p2^2 + q2^2)) - p1^2*q1*(p2 - q2)^2*
          (p2^2 + q2^2) + p1*(p2^4*q1^2 - 2*p2^3*q1^2*q2 + q2^2 + q1^2*q2^4 - 
           2*p2*q2*(2 + q1^2*q2^2) + p2^2*(1 + 2*q1^2*q2^2)) + 
         q1*(p2^4*q1^2 - 2*p2^3*q1^2*q2 + q2^2*(-1 + q1^2*q2^2) + 
           p2^2*(-1 + 2*q1^2*q2^2) + p2*(4*q2 - 2*q1^2*q2^3))) - 
       (p1 - q1)*(8*b2^2*p2^2*q2^2 + (p2^2 + q2^2)*sigm^2)) - 
     (p2 - q2)*(2*b2*b3^2*(p1^4*p2*q2*(1 + p2^2 + 2*p2*q2 + q2^2) + 
         p2*q1^2*q2*(-2 + q1^2*(1 + p2^2 + 2*p2*q2 + q2^2)) + 
         p1^3*q1*(p2^4 + 2*p2^3*q2 + 2*p2*q2^3 + q2^2*(-1 + q2^2) - 
           p2^2*(1 + 2*q2^2)) - 2*p1^2*p2*q2*
          (1 + q1^2*(-1 + 3*p2^2 - 2*p2*q2 + 3*q2^2)) + 
         p1*(p2^4*q1^3 + 2*p2^3*q1^3*q2 + q1^3*q2^2*(-1 + q2^2) - 
           p2^2*q1^3*(1 + 2*q2^2) + 2*p2*q1*q2*(2 + q1^2*q2^2))) + 
       b3^3*(p1 - q1)*(2*p1^4*p2*q2 + 2*p2*q1^2*(-1 + q1^2)*q2 + 
         p1^3*q1*(p2 + q2)^2*(p2^2 + q2^2) + 
         p1*q1*(p2^4*q1^2 + 2*p2^3*q1^2*q2 + q2^2*(-1 + q1^2*q2^2) + 
           2*p2*q2*(2 + q1^2*q2^2) + p2^2*(-1 + 2*q1^2*q2^2)) + 
         2*p1^2*(p2^4*q1^2 - 2*p2^3*q1^2*q2 + 2*p2^2*q1^2*q2^2 + q1^2*q2^4 - 
           p2*(q2 - 2*q1^2*q2 + 2*q1^2*q2^3))) + 
       b2*(-4*p1*p2*q1*q2 + p1^2*(-1 + p2^2 + 2*p2*q2 + q2^2) + 
         q1^2*(-1 + p2^2 + 2*p2*q2 + q2^2))*(2*b2^2*p2*q2 + sigm^2) + 
       b3*(p1 - q1)*(2*b2^2*p2*q2*(-1 + 2*p1*q1*(p2 - q2)^2 + 
           p1^2*(-1 + 2*p2^2 + 4*p2*q2 + 2*q2^2) + 
           q1^2*(-1 + 2*p2^2 + 4*p2*q2 + 2*q2^2)) + 
         (-1 + p1^2 + q1^2 + 3*p1*q1*(p2^2 + q2^2))*sigm^2)))))/
 (2*sqrt(2)*(p1*q1*(2*b2^2*p2*q2 + 4*b2*b3*p2*(p1 - q1)*q2 + 
     2*b3^2*p2*(p1^2 + q1^2)*q2 + sigm^2))^(3/2)*
  sqrt((p2*(b2 + b3*(p1 - q1))^2*q2)/(2*b1^2*p1*q1 + 
     4*b1*b3*p1*q1*(p2 - q2) + 2*b3^2*p1*q1*(p2^2 + q2^2) + sigm^2))*
  (2*b1^2*p1*q1 + 4*b1*b3*p1*q1*(p2 - q2) + 2*b3^2*p1*q1*(p2^2 + q2^2) + 
    sigm^2)^2)
)
}

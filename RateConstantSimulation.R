
#
#  Rovshan G. Sadygov, Ph.D.
#
#  Department of Biochemistry and Molecular Biology
#  
#  The University of Texas Medical Branch
#
#  Galveston, TX
#
#  throughout the code:
#  RA - relative abundance
#
#  M0, M1, M2, M3 - monoisotope, first, second, and third
#  mass isotopomers, correspondingly.
#
#  pX - excess (compared to natural) deuterium enrichment.
#
#  NEH - the number of exchangeable hydrogens in a peptide.
#
#


#
#
#   A function to sample from the Laplace
#   distribution.
#
#   n is the number of sample points (default is set equal to 1);
#
rlaplace <- function (n = 1, m = 0, s = 1)
{
    if (any(s <= 0)) 
        stop("s must be positive")
    q <- runif(n)
    ifelse(q < 0.5, s * log(2 * q) + m, -s * log(2 * (1 - q)) + 
        m)
}



#
#   A function to compute the
#   the RA of the monoisotope,
#   given the NEH, pX, and the
#   RA of the natural monoisotope (I0_0)
#
#   pH is the relative abundance of deuterium in
#   nature.
#

I0_t <- function(I0_0, pX, NEH)
{

     pH = 1.5574*10^(-4)

     I0t = I0_0 * (1 - pX / (1 - pH))^NEH;

     return(I0t)

}

#
#   A function to compute the
#   the RA of the M1,
#   given the NEH, pX, and the
#   RAs of the natural monoisotope (I0_0) and
#   natural M1 (I1_0).
#
#   pH is the relative abundance of deuterium in
#   nature.
#


I1_t <- function(I0_0, I1_0, NEH, pX)
{
     pH = 1.5574*10^(-4)

     I1t = (1 - pX / (1 - pH))^(NEH - 1) * pX * NEH * I0_0 / (1 - pH)^2;

     I1t = I1t + I1_0 * (1 - pX / (1 - pH))^NEH ;

     return(I1t);

}




#
#----------- A function to compute the  the RA of the second   ---------------------
#----------- heavy mass isotoper (M2) at the labeling duration ---------------------
#----------- time t.
#
#   
#            NEH, pX, and the natural RAs: I0_0, I1_0, and I2_0 are inputs.
#
#----------- pX is the deuterium excess labeling at labeling   ------------------
#----------- duration time t.                                  ------------------
#----------- NEH is the number of exchangeable Hydrogens in    ------------------
#----------- the peptide.                                      ------------------
#----------- I0_0, I1_0, and I2_0 are the RIAs of the          ------------------
#----------- monoisotope(M0), first heavy mass isotopomer (M1), ------------------
#----------- and the second heavy mass isotopomer (M3),        ------------------
#----------- respectively, before the start of the labeling.   ------------------
#
#   pH is the relative abundance of deuterium in
#   nature.
#
#

I2_t <- function(I0_0, I1_0, I2_0, pX, NEH)
{

   pH = 1.5574*10^(-4)

   I0t = I0_t(I0_0, pX, NEH);

   I1t = I1_t(I0_0, I1_0, NEH, pX);

   I2t = I0t * (I2_0 / I0_0 - I1_0 * pH * NEH / (I0_0 * (1 - pH)));

   I2t = I2t + I0t * (pH / (1 - pH))^2 * NEH * (NEH + 1) / 2;

   I2t = I2t - I0t * ((pX + pH)/(1 - pH - pX))^2 * NEH * (NEH + 1) / 2;

   I2t = I2t + NEH * (pX + pH) * I1t / (1 - pH - pX);

   return(I2t);
}

#
#  computes the I2_t_tilde - the transformed form
#  which is an exponential decay model
#  Inputs:
#  I0_0, I1_0, I2_0 - are the relative abundances of
#  M0, M1, and M2 mass isotopomers, correspondingly.
#  NEH is the number of exchangeables hydrogens in a peptide
#  pX is the excess enrichment with deuterium of a peptide
#  end of Inputs:
#
#  I2_t_orginal is the original I2_t, that is the RA of the
#  M2 at labeling time t.
#
#   pH is the relative abundance of deuterium in
#   nature.
#
#  I2t_tilde is calculated one term at a time (Eq. 5 of the text).
#
#
I2_t_tilde <- function(I0_0, I1_0, I2_0, pX, NEH)
{
    pH = 1.5574*10^(-4);

    I1t = I1_t(I0_0, I1_0, NEH, pX);

    I2_t_original = I2_t(I0_0, I1_0, I2_0, pX, NEH);

#    print(c(I2_t_original, I1t, I0_0));

    I2t_tilde = I2_t_original + pH * NEH / (1 - pH) * I1_0 * (1 - pX / (1 - pH))^NEH ;


    temp = NEH * (NEH + 1) / 2 * I0_0 * ( (pH / (1 - pH))^2 - ((pX + pH) / (1 - pX - pH))^2 );

    I2t_tilde = I2t_tilde - temp * (1 - pX / (1 - pH))^NEH;

    #print(c(I2t_tilde, temp));

    I2t_tilde = I2t_tilde -  NEH * (pX + pH) * I1t / (1 - pH - pX);

    return (I2t_tilde);
}


#
#
#------------ A function to profile RAs of the first four mass  --------------
#------------ isotopomers for a peptide sequence.               --------------
#------------ inputs: Sequence - peptide sequence               --------------
#------------         pw_limit - the limit for deuterium enrich --------------
#------------                    ment.                          --------------
#
#------------ Nit is the grid length for deuterium enrichment   --------------
#------------ I1, I2, I3 are the orginal RAs.                   --------------
#------------ I1_tilde, I2_tilde are the transformed RAs        --------------
#
#------------ Profile_Peptide("TVLMNPNIASVQTNEVGLK", 0.35)      --------------
#

Profile_Peptide <- function(Sequence, pW_limit=0.35)
{

   Nit = 100;

   Sequence_Elements = ElementalComposition(Sequence);  # call isotope calculator

   NEH = Sequence_Elements$NEH;

   NAll_H =  Sequence_Elements$nH;

   I0_0 = Sequence_Elements$Isotopes[1]; I1_0 = Sequence_Elements$Isotopes[2]; 

   I2_0 = Sequence_Elements$Isotopes[3]; I3_0 = Sequence_Elements$Isotopes[4];


   
   I0 = c(1:Nit);

   I1 = c(1:Nit);

   I2 = c(1:Nit);

   I3 = c(1:Nit);

   I1_tilde = c(1:Nit);

   I2_tilde = c(1:Nit);

   x = c(1:Nit);

   pH = 1.5574*10^(-4);

   for(i in 1:Nit)
   {
      x[i] = 0.000 + (i-1) * pW_limit / Nit;

      I0[i] =  I0_t(I0_0, x[i], NEH);

      I1[i] = I1_t(I0_0, I1_0, NEH, x[i]);

      I2[i] = I2_t(I0_0, I1_0, I2_0, x[i], NEH);

      a = ((pH + x[i])/(1 - pH - x[i]))^3 * NEH * (NEH - 1) * (NEH - 2) / 6;

      b = ((pH + x[i])/(1 - pH - x[i]))^2 * NEH * (NEH - 1) / 2;

      c = ((pH + x[i])/(1 - pH - x[i])) * NEH;

      I3[i] = a * I0[i] + b * (I1[i] - c*I0[i]) + c*(I2[i] - c *(I1[i] - c*I0[i]) -
      	           b * I0[i]);

      I3[i] = I3[i] + I0[i]* (I3_0 / I0_0 - NEH*pH/(1 - pH)* (I2_0/I0_0 -
      (NEH + 1)*pH/(2 * (1 - pH))* (I1_0/I0_0 - NAll_H*pH/(1 - pH)) -
        (3*NAll_H * NEH + 3*NAll_H - NEH^2 - 3 * NEH - 2) * (pH/(1-pH))^2 / 6 ) )


      I1_tilde[i] = I1[i] - (1 - x[i]/(1 - pH))^(NEH - 1) *
	  		 x[i] * NEH * I0_0 / (1 - pH)^2;

      I2_tilde[i] = I2_t_tilde(I0_0, I1_0, I2_0, x[i], NEH);
   }
   

   plot(x, I0, type="l", lwd=2, ylim = c(0, max(I0, I1, I2, I3)),
   	   xlab = "Deuterium Enrichment", ylab="Relative Abundance");

   lines(x, I1, col="blue", lwd=2);

   lines(x, I2, col="purple", lwd=2);

   lines(x, I3, col="gold", lwd=2);

   lines(x, I1_tilde, col="blue", lty=2, lwd=2);

   lines(x, I2_tilde, col="purple",lty = 2, lwd=2);

   legend("topright",
   c(expression(paste("I"[0], "(t)")), expression(paste("I"[1], "(t)")),
   expression(paste(tilde("I")[1], "(t)")), 
   expression(paste("I"[2], "(t)")), expression(paste(tilde("I")[2], "(t)")), 
   expression(paste("I"[3], "(t)"))),
   col=c("black", "blue", "blue", "purple", "purple", "gold"), lwd=c(2,2,2,2,2,2),
   lty=c(1,1,2,1,2,1));

   legend("top", Sequence, cex=0.8, bty="n");

}



#   Rate_Constant_Simulation - function to simulate rate constants   ------------------

#   parameters: deg_rate - degradation rate constant used in         ------------------
#               simulations.
#               noise - the parameters of Laplace distribution:      ------------------
#-------------- noise  is a two-component vector. The first component------------------
#-------------- is the location, the second component is the         ------------------
#-------------- scale parameter of the Laplace random variable.      ------------------
#               BWE - enrichment level of deuterium in diet water
#
#   example: Rate_Constant_Simulation("EPLFGISTGNIITGLAAGAK", deg_rate = 0.110, dates=c(0, 3, 5, 7, 14, 21),noise= c(0, exp(-8)), BWE = 0.03)
#
#
Rate_Constant_Simulation <- function(Sequence, deg_rate, dates,
 			       noise = c(0, exp(-8)), BWE = 0.03, Plot_Figure = FALSE)
{
   # do a number of iteration
   # with different starting values
   # of I_0, and I_asymp
   #

   tt = dates;

   pH = 0.00015574 / (0.99984426 + 0.00015574);

   nDays = length(tt);

   laplace_noise = rlaplace(nDays, noise[1], noise[2]);

   Elements = ElementalComposition(Sequence);

#   print(Elements$mPeptide);

   NEH = Elements$NEH;

   IsotopeDist = Isotopes(Elements$nC, (Elements$nH-NEH), Elements$nN, Elements$nO,
   	       Elements$nS);

   probTemp0 = Isotopes(Elements$nC, Elements$nH, Elements$nN, Elements$nO,
   	       Elements$nS);

   probTemp0 = Elements$Isotopes;

 #  print(Sequence);

   I0_0 = probTemp0[1];   I1_0 = probTemp0[2]; I2_0 = probTemp0[3];

   I0_asymp = probTemp0[1] * (1 - (BWE - pH)/(1 - pH))^NEH;

   I0_t = I0_asymp + (probTemp0[1] - I0_asymp) * exp(-deg_rate * tt);

   pxt = (1 - (I0_t/I0_0)^(1/NEH)) * (1- pH);

   I1_t = c(1:length(tt));

   I1_t_tilde = c(1:length(tt));

   I2_t = c(1:length(tt));

   I2t_tilde = c(1:length(tt));

   for(j in 1:length(tt))
   {
          pX = pH + pxt[j];

          probX = rep(0, 32);

          probX[1] = 1 - pX; probX[2] = pX;

   # The X element with pxt enrichment

          pTemp = probX;
   
          for(i in 1:(NEH-1))
          {
            sum_all = sum(Re(fft(fft(pTemp)*fft(probX)/length(probX), inverse=TRUE)))
	
   	    pTemp = Re(fft(fft(pTemp)*fft(probX)/length(probX), inverse=TRUE)) / sum_all;  
          }

          sum_all = sum(Re(fft(fft(pTemp)*fft(IsotopeDist)/length(probX), inverse=TRUE)))

          probTemp = Re(fft(fft(pTemp)*fft(IsotopeDist)/length(probX), inverse=TRUE)) / sum_all;  

          temp = (1 - pH)^2;

          temp = temp * (probTemp[2]*probTemp0[1] - probTemp0[2]*probTemp[1]);

          pxt_Computed = temp / (NEH*probTemp0[1]*probTemp[1] +
	              (1 - pH)*(probTemp[2]*probTemp0[1] - probTemp0[2]*probTemp[1]))  ;

          I1_t[j] = probTemp[2];

	  I1_t_tilde[j] = I1_t[j] - (1 - pxt_Computed/(1 - pH))^(NEH - 1) *
	  		 pxt_Computed * NEH * I0_0 / (1 - pH)^2;

          I2_t[j] = probTemp[3];

          I2t_tilde[j] = I2_t_tilde(I0_0, I1_0, I2_0, pxt_Computed, NEH);

   }

#
# The next ten lines set up and calculate the isotope distribution
#   a peptide at the plateau of labeling (when the enrichment reaches
#   that of the diet water). The computed RA values are in probTemp.

   pX = BWE;

   probX = rep(0, 32);

   probX[1] = 1 - pX; probX[2] = pX;

   pTemp = probX;
   
   for(i in 1:(NEH-1))
   {
      sum_all = sum(Re(fft(fft(pTemp)*fft(probX)/length(probX), inverse=TRUE)))
	
      pTemp = Re(fft(fft(pTemp)*fft(probX)/length(probX), inverse=TRUE)) / sum_all;  
   }

   sum_all = sum(Re(fft(fft(pTemp)*fft(IsotopeDist)/length(probX), inverse=TRUE)))

   probTemp = Re(fft(fft(pTemp)*fft(IsotopeDist)/length(probX), inverse=TRUE)) / sum_all;  

   I1_asymp = probTemp[2]; 

   I1_asymp_2 = probTemp[2] - (1 - (BWE - pH)/(1 - pH))^(NEH - 1) *
	  		 (BWE - pH) * NEH * I0_0 / (1 - pH)^2;

   I2_asymp = probTemp[3];
   
   I2t_tilde_asymptote = I2_t_tilde(I0_0, I1_0, I2_0, BWE - pH, NEH);

   mod1 <- NULL; mod2 <- NULL; mod3 <- NULL; mod4 <- NULL; mod5 <-NULL;
	     
   try (mod1 <- nls(I1_t + laplace_noise ~ I1_asymp + (I1_0 - I1_asymp) * exp(-k * tt),
           start = list(k = deg_rate)), silent=TRUE );

#   mod1 <- nls(I1_t + laplace_noise ~ I1_asymp + (I1_0 - I1_asymp) * exp(-k * tt),
 #          start = list(k = 0.001));


   try (mod2 <- nls(I0_t + laplace_noise ~ I0_asymp + (I0_0 - I0_asymp) * exp( - k *tt),
   	   start = list(k = deg_rate)), silent=TRUE );
	   
 #  mod2 <- nls(I0_t + laplace_noise ~ I0_asymp + (I0_0 - I0_asymp) * exp( - k *tt),
 #  	   start = list(k = 0.001));

  laplace_noise3 = rlaplace(nDays, noise[1], noise[2]);

   try (mod3 <- nls(I1_t_tilde + laplace_noise3 ~ I1_asymp_2 + (I1_0 - I1_asymp_2) * exp(-k * tt),
           start = list(k = deg_rate)), silent=TRUE );
	   
 #  mod3 <- nls(I1_t_tilde + laplace_noise ~ I1_asymp_2 + (I1_0 - I1_asymp_2) * exp(-k * tt),
 #          start = list(k = deg_rate));

   try (mod4 <- nls(I2t_tilde + laplace_noise ~ I2t_tilde_asymptote +
          (I2_0 - I2t_tilde_asymptote) * exp(-k * tt),  start = list(k = deg_rate)), silent=TRUE );

 #   mod4 <- nls(I2t_tilde + laplace_noise ~ I2t_tilde_asymptote +
 #         (I2_0 - I2t_tilde_asymptote) * exp(-k * tt),  start = list(k = 0.1));

   try (mod5 <- nls(I2_t + laplace_noise ~ I2_asymp +
          (I2_0 - I2_asymp) * exp(-k * tt),  start = list(k = deg_rate)), silent=TRUE );

   #mod5 <- nls(I2_t + laplace_noise ~ I2_asymp +
    #      (I2_0 - I2_asymp) * exp(-k * tt),  start = list(k = 0.1));

  # print(I2_t); print(c(I2_0, I2_asymp))

   a1 = a2 = a3 = a4 = a5 = -1;
   
   if(!is.null(mod1)) a1 = coef(mod1)[[1]];

   if(!is.null(mod2)) a2 = coef(mod2)[[1]];

   if(!is.null(mod3)) a3 = coef(mod3)[[1]];

   if(!is.null(mod4)) a4 = coef(mod4)[[1]];

   if(!is.null(mod5)) a5 = coef(mod5)[[1]];
   
   if(!is.null(mod1) & !is.null(mod5))
   {
      print(c(coef(mod1)[[1]], coef(mod2)[[1]],coef(mod3)[[1]],
              coef(mod4)[[1]], coef(mod5)[[1]])  );
   }
   else
   {
      print("I1 or I2 did not converge");

      print(Sequence);

      return;
   }

   results <- list();


if(Plot_Figure)
{

   ymin = min(I0_t, I1_t, I1_t_tilde);  ymax = max(I0_t, I1_t, I1_t_tilde);

   ymin = min(I0_t, I2t_tilde, I1_t_tilde);  ymax = max(I0_t, I2t_tilde, I1_t_tilde);


   if(!is.null(mod2))
   {
     Plot_Time_Course(tt, I0_t, I0_0, I0_asymp, coef(mod2)[[1]],
                    c(ymin, ymax), FALSE, col="black"); }
   else
      Plot_Time_Course(tt, I0_t, I0_0, I0_asymp, deg_rate,
                    c(ymin, ymax), FALSE, col="black");

   if(!is.null(mod1))
     Plot_Time_Course(tt, I1_t, I1_0, I1_asymp, coef(mod1)[[1]], 
   			c(ymin, ymax), TRUE, color = "blue");

   if(!is.null(mod4))
   Plot_Time_Course(tt, I2t_tilde, I2_0, I2t_tilde_asymptote, coef(mod4)[[1]],
                      c(ymin, ymax), TRUE, color = "cyan");
   if(!is.null(mod3))
   Plot_Time_Course(tt, I1_t_tilde, I1_0, I1_asymp_2, coef(mod3)[[1]], 
   			c(ymin, ymax), TRUE, color="purple");

   if(!is.null(mod5))
   Plot_Time_Course(tt, I2_t, I2_0, I2_asymp, coef(mod5)[[1]], 
   			c(ymin, ymax), TRUE, color="magenta");

   legend("topright",
   c(expression(paste("I"[0], "(t)")),expression(paste("I"[1], "(t)")),
   expression(paste(tilde("I")[1], "(t)")), expression(paste("I"[2], "(t)")),
   expression(paste(tilde("I")[2], "(t)")) ),
   col=c("black", "blue", "purple", "magenta", "cyan"), lwd=c(2,2,2,2,2), cex=0.7);
}

   #print(c(deg_rate, a1, a2, a3, a4, a5) );


 # readline(prompt="Press [enter] to continue")

   results$mPeptide = Elements$mPeptide;

   results$Rates     = c(deg_rate, a2, a1, a3, a5, a4);

   results$NEH = Elements$NEH;

   results$IsoDist = Elements$Isotopes;


   return(results);
   
}

#
#
#-------------------  This function plots multiple RAs on the same plot. It is   -----------------
#------------------   only used inside the function Rate_Constant_Simulation     -----------------
#------------------   to avoid repeatative code lines.                           -----------------
#
Plot_Time_Course <- function(tt, yt, In_0, In_asympt, kdeg = 0.113,
		 ylims=c(0, 1), lines=FALSE, color="black")
{
   if(lines == FALSE)
   {
      #print(ylimits);
      
      plot(tt, In_asympt + (In_0 - In_asympt) * exp( - kdeg *tt), type="l", lwd=2,
       ylim=c(ylims[1], ylims[2]), xlab="Labling Duration (days)",
       			ylab="Relative Abundances");

       points(tt, yt);
    }
    else
    {
        lines(tt, In_asympt + (In_0 - In_asympt) * exp(-tt * kdeg), col=color, lwd=2);

        points(tt, yt, col=color);
    }

}



#
#
#------------- The function calls rlnorm(n, -1.81, 0.914) to simulate  ---------
#------------- the rate constants from mouse liver.                    ---------
#
#
#           Rates_Experim_Noise <- All_Peptides(Peptides_Rates,c(0, 3, 5, 7, 14, 21), c(-0.0035, 0.0159), pW = 0.03)
#
#           All_Rates stores vectors each with six components: true rate, rate from I1(t), rate from I0(t),
#                                                              rate from I1(t) transformed, rate from 
#                                                              transformed I2(t), and the rate from I2(t)
#
#
#

All_Peptides <- function(Peptides_Vector, dates= c(0, 3, 5, 7, 14, 21), noise = c(0, exp(-8)), pW = 0.03)
{
   k = 0;

   All_Rates = matrix(0, nrow = length(Peptides_Vector$Rates), ncol=6);

   Peptide_Masses = c(1:length(Peptides_Vector$Rates) );


   krate = rlnorm(length(Peptides_Vector$Rates),-1.81, 0.914);


   for(i in 1:(length(Peptides_Vector$Rates)))
 #  for(i in 1:10)
   {

      if(i > 1)
      {

        if(toupper(Peptides_Vector$Peptides[i]) != toupper(Peptides_Vector$Peptides[i-1]) )
	{

              k = k + 1;

    	}
      }
      else
      {
          k = k + 1;

      }

   }

   All_Rates = matrix(0, nrow = k, ncol=6);

   All_Isotopes = matrix(0, nrow = k, ncol=32);

   Peptide_Masses = c(1:k);

   Peptide_NEHs = c(1:k);


   krate = rlnorm(length(Peptides_Vector$Rates),-1.81, 0.914);

   print(k);

   k = 0;
   
   for(i in 1:(length(Peptides_Vector$Rates)))
#   for(i in 1:10)
   {

      if(i > 1)
      {

        if(toupper(Peptides_Vector$Peptides[i]) != toupper(Peptides_Vector$Peptides[i-1]) )
	{
	      print(toupper(Peptides_Vector$Peptides[i]));


              k = k + 1;

              temp_results = Rate_Constant_Simulation(toupper(Peptides_Vector$Peptides[i]), krate[i],
 			       dates, noise, BWE = pW);

              All_Rates[k,] = temp_results$Rates;

	      All_Isotopes[k,] = temp_results$IsoDist;

              Peptide_Masses[k] = temp_results$mPeptide;

              Peptide_NEHs[k]   = temp_results$NEH;
         }
      }
      else
      {
          print(toupper(Peptides_Vector$Peptides[i]));

          k = k + 1;


          temp_results = Rate_Constant_Simulation(toupper(Peptides_Vector$Peptides[i]), krate[i], dates,
 			       noise, BWE = pW);

           All_Rates[k,] = temp_results$Rates;

           All_Isotopes[k,] = temp_results$IsoDist;

	   Peptide_Masses[k] = temp_results$mPeptide;

           Peptide_NEHs[k]   = temp_results$NEH;

	   print(c("k = ", i, k));
      }

   }

   results <- list();

   results$Rates = All_Rates[1:k,];

   results$IsoDist = All_Isotopes[1:k,];

   results$noise     = noise;

   results$Peptide_Masses = Peptide_Masses[1:k];

   results$Peptide_NEHs = Peptide_NEHs[1:k];

   print(k);

   return (results)
}
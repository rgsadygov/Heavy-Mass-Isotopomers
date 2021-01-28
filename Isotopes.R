#
#
#
#   A method to compute isotope distribution of a 
#   peptide sequence given its amino acid sequence.
#
#
#

#
#
#------- Letter2num function transforms a letter to a number     -------
#------- It is used to transform the amino acid letter to a      -------
#-------- number.                                                -------
#
Letter2num <- function(Letter)
{
	utf8ToInt(Letter) - utf8ToInt("A") + 1L
}

#
#
#-------- A function - ElementalComposition. Given an amino acid   ------
#---------sequence it returns elemental composition, the number    ------
#-------- of exchangeable hydrogens, and the isotope distribution  ------
#-------- of the amino acid sequence. The results are returned in  ------
#-------- a list. For the isotope distribution, it calls the       ------
#-------- another function, Isotopes.
#
#
#-------- The function first creates the elemental composition     ------
#-------- matrix, ElementalMatrix, which is 30x5 matrix. The       ------
#-------- raws correspond to the amino acids, the columns to the   ------
#-------- atoms.
#-------- #1 Carbon, #2 Hydrogen, #3 Nitrogen, #4 Oxygen, #5 Sulfur
#-------- the NEH values for amino acids are those from the mouse   ------
#          NEH values are stored in the vector NEH_Vector           ------
#
#
#

ElementalComposition <- function(Sequence)
{


    szPeptide = as.character(Sequence);


    mH2OH = 19.01784113;

    mProton = 1.00727642;

    ElementalMatrix = matrix(0, nrow = 30, ncol=5);   # stores the numbers of atom types

    Mass_Vector    = c(rep(0, 30));                   # stores the mass of each amino acid

    NEH_Vector      = c(rep(0, 30));                  # stores the number NEH for each amino acid


    k = Letter2num("A");  ### Alanine

    ElementalMatrix[k, 1]  = 3;  ElementalMatrix[k, 2]  = 5;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 1;

    NEH_Vector[k] = 4;    ### Number of Exchangeable Hydrogens for Ala.

    Mass_Vector[k] = 71.03711378;


    k = Letter2num("G");  ### Glycine

    ElementalMatrix[k, 1]  = 2;  ElementalMatrix[k, 2]  = 3;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 1;

    NEH_Vector[k] = 2.06;

    Mass_Vector[k] = 57.02146372;



    k = Letter2num("S");  ### Serine

    ElementalMatrix[k, 1]  = 3;  ElementalMatrix[k, 2]  = 5;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 2;

    NEH_Vector[k] = 2.61;

    Mass_Vector[k] = 87.03202840;



    k = Letter2num("P");   ### Proline

    ElementalMatrix[k, 1]  = 5;  ElementalMatrix[k, 2]  = 7;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 1;

    NEH_Vector[k] = 2.59;

    Mass_Vector[k] = 97.05276384;



    k = Letter2num("V");     ### Valine

    ElementalMatrix[k, 1]  = 5;  ElementalMatrix[k, 2]  = 9;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 1;

    NEH_Vector[k] = 0.56;

    Mass_Vector[k] = 99.06841390;



    k = Letter2num("T");     ###  Threonine

    ElementalMatrix[k, 1]  = 4;  ElementalMatrix[k, 2]  = 7;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 2;

    NEH_Vector[k] = 0.2;

    Mass_Vector[k] = 101.04767846;



    k = Letter2num("C");      ### Cysteine

    ElementalMatrix[k, 1]  = 3;  ElementalMatrix[k, 2]  = 5;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 1;
    ElementalMatrix[k, 5]  = 1;

    NEH_Vector[k] = 1.62;

    Mass_Vector[k] = 103.00918451;


    k = Letter2num("L");       ### Leucine

    ElementalMatrix[k, 1]  = 6;  ElementalMatrix[k, 2]  = 11;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 1;

    NEH_Vector[k] = 0.6;

    Mass_Vector[k] = 113.08406396;



    k = Letter2num("I");       ### Isoleucine

    ElementalMatrix[k, 1]  = 6;  ElementalMatrix[k, 2]  = 11;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 1;

    NEH_Vector[k]  = 1.0;

    Mass_Vector[k] = 113.08406396;



    k = Letter2num("N");        ### Asparagine, Asn 

    ElementalMatrix[k, 1]  = 4;  ElementalMatrix[k, 2]  = 6;  ElementalMatrix[k, 3]  = 2;  ElementalMatrix[k, 4]  = 2;

    NEH_Vector[k]  = 1.89;

    Mass_Vector[k] = 114.04292744;



    k = Letter2num("D");         ### Aspartate, Asp

    ElementalMatrix[k, 1]  = 4;  ElementalMatrix[k, 2]  = 5;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 3;

    NEH_Vector[k]  = 1.89;

    Mass_Vector[k] = 115.02694302;



    k = Letter2num("Q");         ### Glutamine, Gln

    ElementalMatrix[k, 1]  = 5;  ElementalMatrix[k, 2]  = 8;  ElementalMatrix[k, 3]  = 2;  ElementalMatrix[k, 4]  = 2;

    NEH_Vector[k]  = 3.95;

    Mass_Vector[k] = 128.05857750;



    k = Letter2num("K");          ### Lysine, Lys

    ElementalMatrix[k, 1]  = 6;  ElementalMatrix[k, 2]  = 12;  ElementalMatrix[k, 3]  = 2;  ElementalMatrix[k, 4]  = 1;

    NEH_Vector[k] = 0.54;

    Mass_Vector[k] = 128.09496300;



    k = Letter2num("E");          ### Glutamate, Glu

    ElementalMatrix[k, 1]  = 5;  ElementalMatrix[k, 2]  = 7;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 3;

    NEH_Vector[k] = 3.95;

    Mass_Vector[k] = 129.04259308;



    k = Letter2num("M");          ### Methionine, Met

    ElementalMatrix[k, 1]  = 5;  ElementalMatrix[k, 2]  = 9;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 1;

    ElementalMatrix[k, 5]  = 1;

    NEH_Vector[k] = 1.12;

    Mass_Vector[k] = 131.04048463;



    k = Letter2num("H");           ### Histidine, His

    ElementalMatrix[k, 1]  = 6;  ElementalMatrix[k, 2]  = 7;  ElementalMatrix[k, 3]  = 3;  ElementalMatrix[k, 4]  = 1;

    NEH_Vector[k] = 2.88;

    Mass_Vector[k] = 137.05891186;



    k = Letter2num("F");           ### Phenylalanine, Phe

    ElementalMatrix[k, 1]  = 9;  ElementalMatrix[k, 2]  = 9;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 1;

    NEH_Vector[k]  = 0.32;

    Mass_Vector[k] = 147.06841390;



    k = Letter2num("R");           ### Arginine, Arg

    ElementalMatrix[k, 1]  = 6;  ElementalMatrix[k, 2]  = 12;  ElementalMatrix[k, 3]  = 4;  ElementalMatrix[k, 4]  = 3;

    NEH_Vector[k]  = 3.43;

    Mass_Vector[k] = 156.10111102;



    k = Letter2num("Y");           ### Tyrosine, Tyr

    ElementalMatrix[k, 1]  = 9;  ElementalMatrix[k, 2]  = 9;  ElementalMatrix[k, 3]  = 1;  ElementalMatrix[k, 4]  = 2;

    NEH_Vector[k]  = 0.42;

    Mass_Vector[k] = 163.06332852;


    k = Letter2num("W");          ### Tryptophan, Trp

    ElementalMatrix[k, 1]  = 11;  ElementalMatrix[k, 2]  = 10;  ElementalMatrix[k, 3]  = 2;  ElementalMatrix[k, 4]  = 1;

    NEH_Vector[k] = 0.08;

    Mass_Vector[k] = 186.07931294;
    

    nH = nO = nC = nN = nS = 0;

    nO = 1; nN = 0;

    nH = 3;

    NEH = 0;

    mPeptide = mH2OH;


    for(i in 1:nchar(szPeptide))
    {
	 j = Letter2num(substring(szPeptide,i,i));

         nC = nC + as.numeric(ElementalMatrix[j, 1]);

	 nH = nH + ElementalMatrix[j, 2];
	 
         nN = nN + ElementalMatrix[j, 3];

	 nO = nO + ElementalMatrix[j, 4];

	 nS = nS + ElementalMatrix[j, 5];

	 NEH = NEH + NEH_Vector[j];

	 mPeptide = mPeptide + Mass_Vector[j];
    }

    NEH =as.integer(NEH);

    Composition = list ();

    Composition$nC = nC;   Composition$nH = nH;
    Composition$nN = nN;   Composition$nO = nO;
    Composition$nS = nS;   Composition$NEH = NEH;
    Composition$Sequence = szPeptide;
    Composition$Isotopes = Isotopes(nC, nH, nN, nO, nS);
    Composition$mPeptide = mPeptide;

    Composition
    
}


#
#
#-------- Function - Isotopes. Given the numbers of atoms of       ------
#-------- each type, the function computes the isotope             ------
#-------- distribution of a molecule using FFTs.                   ------
#-------- the results are returned in a 1:32 vector.               ------
#
#

Isotopes <- function(nC, nH, nN, nO, nS)
{

   probC = seq(1:32); probH = seq(1:32); probN = seq(1:32); probO = seq(1:32); probS = seq(1:32);

   for(i in 1:32)
     probC[i] = probH[i] = probN[i] = probO[i] = probS[i] = 0;

#   pH = 0.00015574 / (0.99984426 + 0.00015574);

   # Carbon atom isotopes
   probC[1] = 0.988922 / (0.988922 + 0.011078); probC[2] = 0.011078 / (0.988922 + 0.011078);


   # Hydrogen atom isotopes
   probH[1] = 0.99984426 / (0.99984426 + 0.00015574); probH[2] = 0.00015574 / (0.99984426 + 0.00015574);

 
   # Nitrogen atom isotopes
   probN[1] = 0.996337 / (0.996337 + 0.003663); probN[2] =  0.003663 / (0.996337 + 0.003663);


   # Oxygen isotopes
   probO[1] = 0.9976206 / (0.9976206 + 0.0003790 + 0.0020004);

   probO[2] = 0.000379 /(0.9976206 + 0.0003790 + 0.0020004);

   probO[3] = 0.0020004 /(0.9976206 + 0.0003790 + 0.0020004);


   # Sulfor atom isotopes
   probS[1] = 0.9504074 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

   probS[2] = 0.0074869 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

   probS[3] = 0.0419599 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);

   probS[5] = 0.0001458 / (0.9504074 + 0.0074869 + 0.0419599 + 0.0001458);


   probTemp = probC;
   
   for(i in 1:(nC-1))
   {
	sum_all = sum(Re(fft(fft(probC)*fft(probTemp)/length(probC), inverse=TRUE)))

	probTemp = Re(fft(fft(probC)*fft(probTemp)/length(probC), inverse=TRUE)) / sum_all;    
   }

   probC = probTemp;


   probTemp = probH;
   
   for(i in 1:(nH-1))
   {
	sum_all = sum(Re(fft(fft(probH)*fft(probTemp)/length(probH), inverse=TRUE)))

	probTemp = Re(fft(fft(probH)*fft(probTemp)/length(probH), inverse=TRUE)) / sum_all;    
   }

   probH = probTemp;


   probTemp = probN;
   
   for(i in 1:(nN-1))
   {
	sum_all = sum(Re(fft(fft(probN)*fft(probTemp)/length(probN), inverse=TRUE)))

	probTemp = Re(fft(fft(probN)*fft(probTemp)/length(probN), inverse=TRUE)) / sum_all;    
   }

   probN = probTemp;


   probTemp = probO;
   
   for(i in 1:(nO-1))
   {
	sum_all = sum(Re(fft(fft(probO)*fft(probTemp)/length(probO), inverse=TRUE)))

	probTemp = Re(fft(fft(probO)*fft(probTemp)/length(probO), inverse=TRUE)) / sum_all;    
   }

   probO = probTemp;



   probTemp = probS;

   for(i in 1:(nS-1))
   {
	sum_all = sum(Re(fft(fft(probS)*fft(probTemp)/length(probS), inverse=TRUE)))

	probTemp = Re(fft(fft(probS)*fft(probTemp)/length(probS), inverse=TRUE)) / sum_all;    
   }

   sum_all = sum(Re(fft(fft(probC)*fft(probH)/length(probH), inverse=TRUE)))

   probTemp = Re(fft(fft(probC)*fft(probH)/length(probH), inverse=TRUE)) / sum_all;


   sum_all = sum(Re(fft(fft(probN)*fft(probTemp)/length(probN), inverse=TRUE)))

   probTemp = Re(fft(fft(probTemp)*fft(probN)/length(probN), inverse=TRUE)) / sum_all;    


   sum_all = sum(Re(fft(fft(probO)*fft(probTemp)/length(probO), inverse=TRUE)))

   probTemp = Re(fft(fft(probTemp)*fft(probO)/length(probO), inverse=TRUE)) / sum_all;    


   if(nS > 0)
   {
          sum_all = sum(Re(fft(fft(probS)*fft(probTemp)/length(probS), inverse=TRUE)))

	  probTemp = Re(fft(fft(probTemp)*fft(probS)/length(probS), inverse=TRUE)) / sum_all;    
   }


   probTemp
}


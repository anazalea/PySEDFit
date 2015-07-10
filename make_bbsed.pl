#! /usr/bin/perl 

# version 2008 Aug 28

use PDL;
#use PGPLOT;
#use PDL::Graphics::PGPLOT;
use PDL::AutoLoader;
use PDL::Math;
#use PDL::Astro::Cosmology;
use Astro::Cosmology;

$version = 20090309; 
print "make_bbsed version $version \n"; 




# READ IN AND PARSE THE PARAMETERS THAT WILL CONTROL THE PROCEDURE

$paramfile =@ARGV[0];
#print "Parsing parameter file $paramfile ... \n";

open (PARAMETERS , $paramfile) or die "USAGE: make_bbased parameterfile [-PARAMETER1 value[s] -PARAMETER2 value[s]...]\n";
while ($line = <PARAMETERS>){

  ($line, $junk) = split '#', $line;
  (@objects) = split ' ', $line;
  $parameter = shift @objects;

  if ((length $parameter) > 0){
    $parameter = "-".$parameter; 
    push @paramfileparams, $parameter;
    push @paramfileparams, @objects;
  }
}

@commandlineparams = @ARGV;
shift @commandlineparams;
push @paramfileparams, @commandlineparams;
@params = @paramfileparams;


$lp = -1;
for ($elem=0; $elem < (scalar  @params); $elem++){
  $argv = @params[$elem];
  if ((substr ($argv,0,1)) =~ "-"){   # if first character of this string is a "-" then this string is a parameter name
    $lp++; 
    @parameter[$lp] = "";
  }
  @parameter[$lp] = @parameter[$lp] . $argv. " ";
}


$parameters = "";
$Nlp = $lp;
for ($lp=0; $lp<=$Nlp; $lp++){
  $parameters = $parameters . " " . @parameter[$lp];
  $line = @parameter[$lp]; 
  (@objects) = split ' ', $line;
  $parameter = uc ( shift @objects);
  if ($parameter =~ "MODEL_DIR")       {@MODEL_DIR = @objects;}
  if ($parameter =~ "MODEL_FILE")      {@MODEL_FILE = @objects;}
  if ($parameter =~ "MODEL_NUMBERS")   {@MODEL_NUMBERS = @objects;}
  if ($parameter =~ "BC_PARAM_FILE")   {@BC_PARAM_FILE = @objects;}
  if ($parameter =~ "BC_PARAM_DIR")    {@BC_PARAM_DIR = @objects;}
  if ($parameter =~ "OUTPUT_FILE")     {@OUTPUT_FILE = @objects;}
  if ($parameter =~ "OUTPUT_OVERWRITE"){@OUTPUT_OVERWRITE = @objects;}
  if ($parameter =~ "EXTINCTION_LAW")  {@EXTINCTION_LAW = @objects;}
  if ($parameter =~ "EXTINCTION_R")    {@EXTINCTION_R = @objects;}
  if ($parameter =~ "EXTINCTION_EBV")  {@EXTINCTION_EBV = @objects;}
  if ($parameter =~ "REDSHIFTS")       {@REDSHIFTS = @objects;}
  if ($parameter =~ "COSMOLOGY")       {@COSMOLOGY = @objects;}
  if ($parameter =~ "SPECMODE")        {@SPECMODE = @objects;}
  if ($parameter =~ "FILTER_DIR")      {@FILTER_DIR = @objects;}
  if ($parameter =~ "FILTER_NAME")     {push @FILTER_NAMES, @objects;}
  if ($parameter =~ "COSMIC_OPACITY")  {@COSMIC_OPACITY = @objects;}
  if ($parameter =~ "COSMIC_OP_LAW")   {@COSMIC_OP_LAW = @objects;}
}


if (scalar @MODEL_FILE == 0)     { print "ERROR:  No model file specified\n"; exit;}
if (scalar @MODEL_NUMBERS == 0)  { print "ERROR:  No models specified\n";exit;}
if (scalar @EXTINCTION_LAW == 0) { print "ERROR:  No extinction law specified\n";exit;}
if (scalar @EXTINCTION_R == 0)   { print "ERROR:  No extinction R specified\n";exit;}
if (scalar @EXTINCTION_EBV == 0) { print "ERROR:  No extinction E(B-V) specified\n";exit;}
if (scalar @REDSHIFTS == 0)      { print "ERROR:  No redshift range specified\n";exit;}
if (scalar @COSMOLOGY <3 )       { print "ERROR:  No cosmology specified\n";exit;}
if (($SPECMODE =~ "mag") & (scalar @FILTER_NAMES == 0)) { print "ERROR:  Magnitudes requested but no filters specified\n";exit;}
if (scalar @COSMIC_OPACITY == 0) { print "ERROR:  No cosmological opacity specified\n";exit;}


$MODEL_DIR = @MODEL_DIR[0];
$MODEL_FILE = @MODEL_FILE[0];
#print "Reading model array $MODEL_FILE ...\n";
@modelArray = rcols $MODEL_DIR . $MODEL_FILE;
$Nmodels = scalar @modelArray - 1;

print " N models = $Nmodels \n";

$MODELstyle = shift @MODEL_NUMBERS;	
if ($MODELstyle =~ 'value'){
  @MODEL_LIST = @MODEL_NUMBERS;}
elsif ($MODELstyle =~ 'range'){
  $min = @MODEL_NUMBERS[0]; $max = @MODEL_NUMBERS[1]; $step = @MODEL_NUMBERS[2];
  for ($value=$min; $value<=$max; $value=$value+$step){
    push @MODEL_LIST, $value;}
}
elsif ($MODELstyle =~ 'all'){
  $min = 1; $max = $Nmodels;
  @MODEL_LIST = $min..$max;
}

if (scalar @BC_PARAM_FILE != 0){
  $BC_PARAM_FILE = @BC_PARAM_FILE[0];
  $BC_PARAM_DIR = @BC_PARAM_DIR[0];
  ($BCmodelages, $BCmodelMstars, $BCmodelMgas, $BCmodelMgalaxy, $BCmodelSFR)=rcols $BC_PARAM_DIR.$BC_PARAM_FILE, 0, 6,7,8,9;
}


if (scalar @OUTPUT_FILE != 0){
  $OUTPUT_FILE = @OUTPUT_FILE[0];
  $OUTPUT_OVERWRITE = @OUTPUT_OVERWRITE[0];
  if (-e $OUTPUT_FILE){
    if (lc($OUTPUT_OVERWRITE) !~ 'yes'){
      print "ERROR:  Specified output file already exists.  Set OUTPUT_OVERWRITE to 'yes' if you want me to overwrite it. \n";
      exit; }
  }
  open (OUTFILE , "> $OUTPUT_FILE");
}


# parsing the dust law selection
$EXTINCTION_LAW = @EXTINCTION_LAW[0];
if (
    (lc($EXTINCTION_LAW) !~ 'calzetti2000')
    & (lc($EXTINCTION_LAW) !~ 'calzetti1997')
    & (lc($EXTINCTION_LAW) !~ lc('MW_fitz'))
    & (lc($EXTINCTION_LAW) !~ lc('LMC_fitz'))
    & (lc($EXTINCTION_LAW) !~ lc('30Dor_fitz'))
    & (lc($EXTINCTION_LAW) !~ lc('prevot_fitz'))
   )
  {
    print "\n";
    print "ERROR:  The requsted extinction law '$EXTINCTION_LAW' is not supported\n";
    print "        The following are currently supported (please specify in full one of the names in 1st column):\n";
    print "              calzetti2000 		(Calzetti 2000 starburst law)\n";
    print "              calzetti1997		(Calzetti 1997 starburst law)\n";
    print "              LMC_fitzpatrick1986 	(Fitzpatrick 1986 LMC)\n";
    print "              MW_fitzpatrick1986	(Milky Way: Seaton via Fitzpatrick)\n";
    print "              30Dor_fitzpatrick1986	(Fitzpatrick 1986 30 Doradus)\n";
    print "              SMC_prevot1984		(Prevot et al 1984 SMC)\n";
    exit;
 }



$EXTINCTION_R = @EXTINCTION_R[0];

$EBVstyle = shift @EXTINCTION_EBV;	
if ($EBVstyle =~ 'value'){
  @EBV_LIST = @EXTINCTION_EBV;}
elsif ($EBVstyle =~ 'range'){
  $min = @EXTINCTION_EBV[0]; $max = @EXTINCTION_EBV[1]; $step = @EXTINCTION_EBV[2];
  for ($value=$min; $value<=$max; $value=$value+$step){
    push @EBV_LIST, $value;}
}

$REDSHIFTSstyle = shift @REDSHIFTS;	
if ($REDSHIFTSstyle =~ 'value'){
  @REDSHIFTS_LIST = @REDSHIFTS;}
elsif ($REDSHIFTSstyle =~ 'range'){
  $min = @REDSHIFTS[0]; $max = @REDSHIFTS[1]; $step = @REDSHIFTS[2];
  for ($value=$min; $value<=$max; $value=$value+$step){
    push @REDSHIFTS_LIST, $value;}
}

$Omega_matter = shift @COSMOLOGY;
$Omega_lambda = shift @COSMOLOGY;
$Ho = shift @COSMOLOGY;
#$cosmo=PDL::Astro::Cosmology->new(matter=>$Omega_matter, lambda=>$Omega_lambda, H0=>$Ho);
$cosmo=Astro::Cosmology->new(matter=>$Omega_matter, lambda=>$Omega_lambda, H0=>$Ho);




# parsing the cosmic opacity law selection

if (scalar @COSMIC_OP_LAW == 0){  # if none specified then default to Madau 1995
  print "\nWARNING:  No cosmic opacity law specified.  Defaulting to Madau (1995) \n";
  $cosmic_op_law = 'madau';
}
else {
  $cosmic_op_law = @COSMIC_OP_LAW[0];
  $cosmic_op_law = lc $cosmic_op_law;
  if (  ( $cosmic_op_law !~ 'madau')
      & ( $cosmic_op_law !~ 'inoue') ) 
    {
      print "\n";
      print "ERROR:  The requested cosmic opacity law '$cosmic_op_law' is not supported\n";
      print "        The following are currently supported (please specify in full one of the names in 1st column):\n";
      print "              madau 		(Madau 1995)\n";
      print "              inoue 		(Inoue et al. 2005)\n";
      exit;
    }
}




$specmode = shift @SPECMODE;
if ($specmode =~ 'spec'){
  #dev '/xs';
  #dev '/vcps';
  #plotting_setup (6,1);
  #env 0, 30000, 60,0;
}

$filter_dir = shift @FILTER_DIR;
$Nfilters = scalar @FILTER_NAMES;
@filterNames = @FILTER_NAMES;

@COSMIC_OPACITY_LIST = @COSMIC_OPACITY;


$modelLambda = @modelArray[0];

# PRINT OUTPUT FILE HEADER:
print OUTFILE "# parameter file =    $paramfile \n";
print OUTFILE "# model =             $MODEL_FILE \n";
print OUTFILE "# cosmology =         $cosmo \n";
print OUTFILE "# dust law =          $EXTINCTION_LAW \n";
print OUTFILE "# intergal opacity =  $cosmic_op_law \n";
print OUTFILE "# full params list =  $parameters \n";
print OUTFILE "# \n";
if ($specmode =~ 'spec'){
  print OUTFILE "# wavelenghts:\n";
  for ($element=0; $element<nelem $modelLambda; $element++){
    print OUTFILE "     ".$modelLambda->index($element)." ";
  }
  print OUTFILE "\n";
}
print OUTFILE "# Columns:\n";  $columnCounter=0;
print OUTFILE "#  (".$columnCounter.") z \n"; $columnCounter++;
print OUTFILE "#  (".$columnCounter.") model \n"; $columnCounter++;
print OUTFILE "#  (".$columnCounter.") E(B-V) \n"; $columnCounter++;
print OUTFILE "#  (".$columnCounter.") intergalactic attenuation \n"; $columnCounter++;
if (scalar @BC_PARAM_FILE != 0){
  print OUTFILE "#  (".$columnCounter.") model age (log Gyr) \n"; $columnCounter++;
  print OUTFILE "#  (".$columnCounter.") Mstars (Msun) \n"; $columnCounter++;
  print OUTFILE "#  (".$columnCounter.") SFR (Msun/yr) \n"; $columnCounter++;
}
if ($specmode =~ 'mag'){
  for ($line=0; $line<$Nfilters; $line++){
    print OUTFILE "#  (". $columnCounter .") ". @filterNames[$line]. " (AB mag)\n"; $columnCounter++;
  }
}
if ($specmode =~ 'spec'){
  print OUTFILE "#  (".$columnCounter."+) spectrum in AB units \n"; $columnCounter++;
}






# End of preliminaries
##############################################
# Start of the meat of the program









#READ IN THE FTCs AND ALSO RECALCULATE THE SCALE INTO nu UNITS (Hz):
# (filters are stored as filter [0-N] sequentially); 

for ($filter=0; $filter<$Nfilters; $filter++){
  $filterfile = $filter_dir . $filterNames[$filter] ;

  # read in the filters:
  (@filterLambda[$filter], @filterTrans[$filter]) = rcols $filterfile,0,1;

  # translate wavelengths into corresonding frequencies (in Hz):
  @filterNu[$filter] = 3e18/@filterLambda[$filter]; 

  # resample the nu scale onto a uniform grid in frequency (needed later for integrating):
  ($minNu, $maxNu) = minmax @filterNu[$filter];
  @filterNuResampled[$filter] = $minNu + (sequence(100)/100)*($maxNu-$minNu);
  @filterTransNuResampled[$filter] = interpol @filterNuResampled[$filter], @filterNu[$filter], @filterTrans[$filter];

  # plot the FTCs (this can be commented out):
#  line (@filterLambda[$filter], 40-80*@filterTrans[$filter],{color=>cyan});
}






foreach $z (@REDSHIFTS_LIST){



############################################################################################
#
# CALCULATE THE COSMIC OPACITY ATTENUATION
#
############################################################################################

  foreach $cosmic_opacity_amount (@COSMIC_OPACITY){
    # CALCULATE THE COSMIC OPACITY AT THIS REDSHIFT
    $z_em=$z;
    $wav_em = $modelLambda; 
    $wav_obs = $wav_em*(1+$z_em);

 
    $cosmic_trans = ones nelem $wav_obs;
#    $cosmic_op_law = @COSMIC_OP_LAW[0];
#    if ($cosmic_op_law == 0){
#	$cosmic_op_law = madau;
#    }
    if ($cosmic_op_law =~ "madau") {
      $wavL = 912;  # 912A
      $A_j = pdl(0.0036, 0.0017, 0.0012, 0.00093); # A_j coefficients for lines
      $wav_j = pdl (1216, 1026, 973, 950); # wavelenghts of LyA, B, G, D

	# line blanketing:
	foreach $line (0,1,2,3){
	    $tau_thisline = $A_j->index($line)*($wav_obs/$wav_j->index($line))**3.46;
	    $mask = $wav_em < $wav_j->index($line);
	    $tau_thisline = $tau_thisline * $mask;
	    $cosmic_trans = $cosmic_trans * 2.71**-$tau_thisline;
	}
	
	# photoelectric:
	$x_em = 1+$z_em;
	$z_c = $wav_obs/$wavL -1;
	$mask = $wav_obs > $wavL;
	$z_c = $z_c * $mask;
	$x_c = 1+$z_c;
	
	$tau_photoelectric = 0.25*$x_c**3*($x_em**0.46-$x_c**0.46)
	    + 9.4*$x_c**1.5*($x_em**0.18-$x_c**0.18)
	    - 0.7*$x_c**3*($x_c**-1.32-$x_em**-1.32)
	    - 0.023*($x_em**1.68-$x_c**1.68);
	
	$mask = $wav_obs<$wavL*(1+$z_em); 
	$tau_photoelectric = $tau_photoelectric*$mask;
	
	$cosmic_trans = $cosmic_trans * 2.71**-$tau_photoelectric; 
	$cosmic_trans = $cosmic_trans**$cosmic_opacity_amount;
    }

    elsif ($cosmic_op_law =~ "inoue"){
	$cosmo_trans = igm_aki($z_em, $wav_obs);
	$cosmic_trans = $cosmic_trans * $cosmo_trans;
	$cosmic_trans = $cosmic_trans**$cosmic_opacity_amount;
    }

# wcols $wav_obs, $cosmic_trans;
    
    foreach $model (@MODEL_LIST){
      foreach $Ebv (@EBV_LIST){
#	$modelLambda = @modelArray[0];
	$L_lambda = @modelArray[$model];
      	# $modelWav is in Angstroms and
       	# $modelFlux is in Lsun/s


#################################################################################################################
#
# APPLY DUST REDDENING CURVE TO THE MODEL
#
#################################################################################################################


# Apply the Calzetti2000 law as defined in Calzetti et al 2000 (ApJ 533 682) Equations 2,3,4.  This is an extrapolation below 1200A and above 2.2um
	if (lc($EXTINCTION_LAW) =~ 'calzetti2000'){

	  $modelLambda_microns = $modelLambda/10000;
	
	  $k_short= 2.659*(-2.156 + 1.509/$modelLambda_microns -0.198/$modelLambda_microns**2 + 0.011/$modelLambda_microns**3) + $EXTINCTION_R;
	  $k_long = 2.659*(-1.857 + 1.040/$modelLambda_microns) + $EXTINCTION_R;

	  $k_short = where ($k_short, $modelLambda_microns < 0.63);
	  $k_long  = where ($k_long, $modelLambda_microns >= 0.63);
	  $k_calzetti = append $k_short, $k_long; 

	  $Fobs_over_Fintrinsic = 10**(-0.4*$Ebv*$k_calzetti);

	  $L_lambda = $L_lambda * $Fobs_over_Fintrinsic;
	}


# Apply the Calzetti1997 law as defined in Calzetti et al 1997 (astro-ph/9706121) Eqns 1 and 2.  This is an extrapolation below 1200A and above 1um
	elsif (lc($EXTINCTION_LAW) =~ 'calzetti1997'){
	  $modelLambda_microns = $modelLambda/10000;

	  $k_short = 2.656*(-2.156 + 1.509/$modelLambda_microns -0.198/$modelLambda_microns**2 + 0.011/$modelLambda_microns**3)+4.88;
	  $k_long = (((1.86 - 0.48/$modelLambda_microns)/$modelLambda_microns-0.1)/$modelLambda_microns)+1.73;

	  $k_short = where ($k_short, $modelLambda_microns < 0.63);
	  $k_long  = where ($k_long, $modelLambda_microns >= 0.63);
	  $k_calzetti = append $k_short, $k_long; 

	  $Fobs_over_Fintrinsic = 10**(-0.4*$Ebv*$k_calzetti);

	  $L_lambda = $L_lambda * $Fobs_over_Fintrinsic;
	}


# Apply the LMC, MW, or 30 Doradus extinction law as defined in Fitzpatrick 1986
	elsif (lc($EXTINCTION_LAW) =~ lc('fitz'))   {

	  $modelLambda_microns = $modelLambda/10000;
	  $lambda_inv = 1/$modelLambda_microns;

	  if (lc($EXTINCTION_LAW) =~ lc('LMC')){
	    $c1= -0.69;
	    $c2= 0.89;
	    $c3= 2.55;
	    $c4= 0.50;
	    $lambda_naught_inv= 4.608;
	    $gamma= 0.994;
	  }

	  if (lc($EXTINCTION_LAW) =~ lc('MW')){
	    $c1= -0.38;
	    $c2= 0.74;
	    $c3= 3.96;
	    $c4= 0.26;
	    $lambda_naught_inv= 4.595;
	    $gamma= 1.051;
	  }

	  if (lc($EXTINCTION_LAW) =~ lc('30Dor')){
	    $c1= -2.19;
	    $c2= 1.39;
	    $c3= 1.49;
	    $c4= 0.43;
	    $lambda_naught_inv= 4.606;
	    $gamma= 0.894;
	  }

  #	  compute $Elv_over_Ebv where $Elv_over_Ebv = E(lambda-V)/E(B-V) : 
	  $c4_A = where ( (zeroes(nelem $lambda_inv)), $lambda_inv < 5.9);
	  $c4_B = where ( ($c4*ones(nelem $lambda_inv)),   $lambda_inv >= 5.9);
	  $c4_v = append $c4_B, $c4_A;

	  $Elv_over_Ebv = $c1 
	    + $c2*$lambda_inv 
	      + $c3/(($lambda_inv - (($lambda_naught_inv**2)/$lambda_inv))**2 + $gamma**2)
		+ $c4_v*(0.539*($lambda_inv-5.9)**2 + 0.0564*($lambda_inv-5.9)**3);

	  $Av_over_Ebv = 3.1; # Assumes that  R(V) = A(V)/E(B-V) = 3.1
	  $tau = ($Elv_over_Ebv+$Av_over_Ebv) * $Ebv/1.086;
	  $L_lambda = $L_lambda * exp(-1*$tau);

	}		

# Apply the Prevot et al 1984 SMC law:
	elsif (lc($EXTINCTION_LAW) =~ 'prevot1984'){
	  $modelLambda_microns = $modelLambda/10000;
	  $lambda_inv = 1/$modelLambda_microns;

  #	  compute $Elv_over_Ebv where $Elv_over_Ebv = E(lambda-V)/E(B-V) (uses an interpolation of the data in Table 2 of Prevot et al 1984): 

	  $lambda_inv_prevot = pdl (7.84, 7.52, 7.23, 6.98, 6.72,
	  6.48, 6.27, 6.07, 5.88, 5.70, 5.52, 5.38, 5.24, 5.00, 4.73,
	  4.50, 4.28, 4.09, 3.92, 3.75, 3.60, 3.46, 3.34, 3.22, 2.70,
	  2.35, 1.89, 0);

	  $Elv_over_Ebv_prevot = pdl(13.54, 12.52, 11.51, 10.80, 9.84, 9.28,
	  9.06, 8.49, 8.01, 7.71, 7.17, 6.90, 6.76, 6.38, 5.85, 5.30,
	  4.53, 4.24, 3.91, 3.49, 3.15, 3.00, 2.65, 2.29, 1.67, 1.00,
	  0.00, 0.00);

	  $Elv_over_Ebv = interpol ($lambda_inv,  $lambda_inv_prevot, $Elv_over_Ebv_prevot);

	  $Av_over_Ebv = 3.1; # Assumes that  R(V) = A(V)/E(B-V) = 3.1
	  $tau = ($Elv_over_Ebv+$Av_over_Ebv) * $Ebv/1.086;
	  $L_lambda = $L_lambda * exp(-1*$tau);

	}



	# APPLY COSMOLOGICAL OPACITY (COMUTED EARLIER) TO THE MODEL:
	$L_lambda = $L_lambda * $cosmic_trans;


	# CONVERT TO CGS UNITS AND ALSO TO FREQUENCY (RATHER THAN WAVELENGTH) SCALE:
	$L_lambda = $L_lambda * 3.826e33; # converts L_lambda from Lsun/A to erg/s/A
	$L_nu = $L_lambda * $modelLambda**2 / 3e18; # converts L_lambda to L_nu
	$modelNu = 3e18/$modelLambda;  # calculate the nu (frequencies in Hz) from the lambda(A)

	# FOR Z=0, CONVERT THE SPECTRUM INTO ABSOLUTE MAGS (D=10pc):
	if ($z==0){
	  $dist = 3.085e19; # 10 paresecs in cm
	  $f_nu = $L_nu / (4*3.141 * $dist**2);
	  $modelLambdaZ = $modelLambda;
	  $modelNuZ = $modelNu;
	}
	# OTHERWISE, REDSHIFT THE SPECTRUM IN PREP FOR APPARENT MAGS:
	if ($z>0){
	  # the luminosity distance is: 
	  $D_L = $cosmo->lum_dist($z); # luminosity distance in Mpc
	  $d_L = $D_L * 3.085e24; #luminosity distance in cm
	  
          # the redshifted spectrum is: 
	  $f_nu = ((1+$z) * $L_nu / (4*3.141 * $d_L**2));  # redshifted spectrum in erg/s/Hz
	  $modelLambdaZ = $modelLambda*(1+$z);
	  $modelNuZ = $modelNu/(1+$z);
       	}
	
	printf OUTFILE "%7.4f  %3i %5.2f %4.2f    ",$z, $model, $Ebv, $cosmic_opacity_amount;
	if (scalar @BC_PARAM_FILE != 0){
	  printf OUTFILE "%9.6f  %9.4e  %9.4e    ",
	    $BCmodelages->index($model-1),
	      $BCmodelMstars->index($model-1),
		$BCmodelSFR->index($model-1);
	}

	if ($specmode =~ 'spec'){
	  # THIS IS THE REDSHIFTED _SPECTRUM_ in AB magnitudes.  SUITABLE FOR PLOTTING AND ALL KINDS OF CALCULATIONS:
	  $f_nu_AB = -(2.5*log10($f_nu) + 48.60);
	  #line $modelLambdaZ, $f_nu_AB; hold;
#	  wcols $modelLambdaZ, $f_nu_AB; hold;
	  for ($element=0; $element < nelem $f_nu_AB; $element++){
	    print OUTFILE $f_nu_AB->index($element); print OUTFILE " "
	  }
	  print OUTFILE "\n";
	}

	if ($specmode =~ 'mag'){
	  # CONVOLVE WITH FTCs TO MAKE MAGNITUDES:
	  for ($filter=0; $filter<$Nfilters; $filter++){

	    #resample the spectrum onto the same (resampled) frequency scale as the FTC:
	    $f_nu_Resampled = interpol @filterNuResampled[$filter], $modelNuZ, $f_nu;

	    # calculate the integrals ($bottom is just the integral over the FTC, $top is the integral of the FTC*spec):
	    $top = sumover ($f_nu_Resampled*@filterTransNuResampled[$filter]/@filterNuResampled[$filter]);
	    $bottom = sumover (@filterTransNuResampled[$filter]/@filterNuResampled[$filter]);

	    @mAB[$filter] = (-2.5*log10($top/$bottom) -48.60);
	    printf OUTFILE "%7.3f ", @mAB[$filter];
	  }
	}
	print OUTFILE "\n";
      }
    }
  }
}





###############SUB ROUTINES################

sub ldp {
    my ($variable) = @_;
    my ($nelem, $i, $temp_val, @values);
    $nelem = nelem $variable;
    while ($i<$nelem){
	$temp_val = $variable($i);
	chop $temp_val;
	$temp_val = substr $temp_val,1;
	push @values, $temp_val;
	$i++;
    }
    return @values;
}

sub igm_aki($) {
    my ($ca, $gamma1, $gamma2, $ca1, $ca2, $zc,$zs, $zmin);
    my ($beta, $delta, $nlow, $factor, $b, $lam, $tau, $c_l, @waveform);
    my ($nu, $nu1, $f1, $nud, $cs0, $ctau, $index, @trans, $waveform);
    my ($i, $j);
    
    # Parameters of line distribution function by Kim et al.(2002)
    $ca = 40.0;
    $gamma1=0.2;
    $gamma2=2.5;
    $zc=1.1;
    $ca1=$ca*(1.0+$zc)**(-$gamma1); # pow((1.0+$zc),-$gamma1);
    $ca2=$ca*(1.0+$zc)**(-$gamma2); # pow((1.0+$zc),-$gamma2);
    $beta=1.5; # 1 < beta < 2;
    $nlow=4.365e13; # cm-2
    $sign = 0.;
#    $factor=exp(sqrt 3.14)*$nlow**($beta-1.0); # pow($nlow,($beta-1.0));
    $factor=exp(gammln(2.0-$beta))*($nlow**($beta-1.0));  # <-- modified by iwata
    #printf "factor = %e\n", $factor;

    # Mean Doppler parameter [km/s]
    $b=30.0;
    
    # Atomic data input
    $c_l = 2.9979e10; # cm s-1
    
    # hydrogen: Lyman series
    # line frequency [Hz]
    $nu1[2]=1.0e8*$c_l/1215.67;
    $nu1[3]=1.0e8*$c_l/1025.72;
    $nu1[4]=1.0e8*$c_l/972.537;
    $nu1[5]=1.0e8*$c_l/949.743;
    $nu1[6]=1.0e8*$c_l/937.803;
    $nu1[7]=1.0e8*$c_l/930.748;
    $nu1[8]=1.0e8*$c_l/926.226;
    $nu1[9]=1.0e8*$c_l/923.150;
    $nu1[10]=1.0e8*$c_l/920.963;
    $nu1[11]=1.0e8*$c_l/919.352;
    $nu1[12]=1.0e8*$c_l/918.129;
    $nu1[13]=1.0e8*$c_l/917.181;
    $nu1[14]=1.0e8*$c_l/916.429;
    $nu1[15]=1.0e8*$c_l/915.824;
    $nu1[16]=1.0e8*$c_l/915.329;
    $nu1[17]=1.0e8*$c_l/914.919;
    $nu1[18]=1.0e8*$c_l/914.576;
    $nu1[19]=1.0e8*$c_l/914.286;
    $nu1[20]=1.0e8*$c_l/914.039;
    $nu1[21]=1.0e8*$c_l/913.826;
    $nu1[22]=1.0e8*$c_l/913.641;
    $nu1[23]=1.0e8*$c_l/913.480;
    $nu1[24]=1.0e8*$c_l/913.339;
    $nu1[25]=1.0e8*$c_l/913.215;
    $nu1[26]=1.0e8*$c_l/913.104;
    $nu1[27]=1.0e8*$c_l/913.006;
    $nu1[28]=1.0e8*$c_l/912.918;
    $nu1[29]=1.0e8*$c_l/912.839;
    $nu1[30]=1.0e8*$c_l/912.768;
    $nu1[31]=1.0e8*$c_l/912.703;
    $nu1[32]=1.0e8*$c_l/912.645;
    $nu1[33]=1.0e8*$c_l/912.592;
    $nu1[34]=1.0e8*$c_l/912.543;
    $nu1[35]=1.0e8*$c_l/912.499;
    $nu1[36]=1.0e8*$c_l/912.458;
    $nu1[37]=1.0e8*$c_l/912.420;
    $nu1[38]=1.0e8*$c_l/912.385;
    $nu1[39]=1.0e8*$c_l/912.353;
    $nu1[40]=1.0e8*$c_l/912.324;
    
    # absorption oscillator strength
    $f1[2]=0.4162;
    $f1[3]=7.910e-2;
    $f1[4]=2.899e-2;
    $f1[5]=1.394e-2;
    $f1[6]=7.799e-3;
    $f1[7]=4.814e-3;
    $f1[8]=3.183e-3;
    $f1[9]=2.216e-3;
    $f1[10]=1.605e-3;
    $f1[11]=1.201e-3;
    $f1[12]=9.214e-4;
    $f1[13]=7.227e-4;
    $f1[14]=5.774e-4;
    $f1[15]=4.686e-4;
    $f1[16]=3.856e-4;
    $f1[17]=3.211e-4;
    $f1[18]=2.702e-4;
    $f1[19]=2.296e-4;
    $f1[20]=1.967e-4;
    $f1[21]=1.698e-4;
    $f1[22]=1.476e-4;
    $f1[23]=1.291e-4;
    $f1[24]=1.136e-4;
    $f1[25]=1.005e-4;
    $f1[26]=8.923e-5;
    $f1[27]=7.970e-5;
    $f1[28]=7.144e-5;
    $f1[29]=6.429e-5;
    $f1[30]=5.806e-5;
    $f1[31]=5.261e-5;
    $f1[32]=4.782e-5;
    $f1[33]=4.360e-5;
    $f1[34]=3.986e-5;
    $f1[35]=3.653e-5;
    $f1[36]=3.357e-5;
    $f1[37]=3.092e-5;
    $f1[38]=2.854e-5;
    $f1[39]=2.640e-5;
    $f1[40]=2.446e-5;
    
    for($i=2; $i<=40; $i++){
	$nud[$i]=($b*1.0e5/$c_l)*$nu1[$i];
	$cs0[$i]=0.01497*$f1[$i]/$nud[$i]; # cm^2
	
	# 0.01497=sqrt(pi)*e^2/(me c) [cm^2 s^-1]
    }
    
    # for ionizing photons
    $nu1[1]=1.0e8*$c_l/911.50;
    $cs0[1]=6.3e-18; # cm^2
    
    $zs = $_[0];
    $waveform = $_[1]; $snum = nelem $waveform;
    @waveform = ldp($waveform);
        
    for ($i = 0; $i < $snum; $i++) {
	$lam = @waveform[$i]; # observed-frame wavelength in AA
	$nu = 2.9979e18/$lam;
	
	# Lyman serise lines
	$delta = 2.6;
	for($j=2; $j<=40; $j++){
	    if(($nu1[$j]/$nu)>=(1.0+$zs) || ($nu1[$j]/$nu)<=1.0) {
		$tau[$j]=0.0;
	    } elsif (($nu1[$j]/$nu)<(1.0+$zc)) {
		$tau[$j]=$delta*$factor*$ca1*$cs0[$j]**($beta-1.0)*($nud[$j]/$nu1[$j])*($nu1[$j]/$nu)**($gamma1+1.0);
	    } else {
		$tau[$j]=$delta*$factor*$ca2*$cs0[$j]**($beta-1.0)*($nud[$j]/$nu1[$j])*($nu1[$j]/$nu)**($gamma2+1.0);
	    }
	}
	
	# Lyman continuum
	$zmin=$nu1[1]/$nu-1.0;
	if($zmin<0.0) {
	    $zmin=0.0;
	}
	
	if($zmin>=$zs) {
	    $tau[1]=0.0;
	} 
	elsif($zc>=$zs){
	    $index=$gamma1-3.0*$beta+4.0;
	    if ($index==0.0) {
		$tau[1]=$ca1*$factor*$cs0[1]**($beta-1.0)*($nu1[1]/$nu)**(3.0*($beta-1.0))*log((1.0+$zs)/(1.0+$zmin));
	    } 
	    else {
		$tau[1]=$ca1*$factor*$cs0[1]**($beta-1.0)*($nu1[1]/$nu)**(3.0*($beta-1.0))*((1.0+$zs)**$index-(1.0+$zmin)**$index)/$index;
	    }
	} 
	else {
	    if($zmin>=$zc){
		$index=$gamma2-3.0*$beta+4.0;
		if($index==0.0) {
		    $tau[1]=$ca2*$factor*$cs0[1]**($beta-1.0)*($nu1[1]/$nu)**(3.0*($beta-1.0))*log((1.0+$zs)/(1.0+$zmin));
		}
		else {
		    $tau[1]=$ca2*$factor*$cs0[1]**($beta-1.0)*($nu1[1]/$nu)**(3.0*($beta-1.0))*((1.0+$zs)**$index-(1.0+$zmin)**$index)/$index;
		}
	    } 
	    else {
		$index=$gamma2-3.0*$beta+4.0;
		if ($index==0.0) {
		    $tau[1]=$ca2*$factor*$cs0[1]**($beta-1.0)*($nu1[1]/$nu)**(3.0*($beta-1.0))*log((1.0+$zs)/(1.0+$zc));
		} 
		else {
		    $tau[1]=$ca2*$factor*$cs0[1]**($beta-1.0)*($nu1[1]/$nu)**(3.0*($beta-1.0))*((1.0+$zs)**$index-(1.0+$zc)**$index)/$index;
		}
		$index=$gamma1-3.0*$beta+4.0;
		if ($index==0.0) {
		    $tau[1]+= $ca1*$factor*$cs0[1]**($beta-1.0)*($nu1[1]/$nu)**(3.0*($beta-1.0))*log((1.0+$zc)/(1.0+$zmin));
		}
		else {
		    $tau[1]+= $ca1*$factor*$cs0[1]**($beta-1.0)*($nu1[1]/$nu)**(3.0*($beta-1.0))*((1.0+$zc)**$index-(1.0+$zmin)**$index)/$index;
		}
	    }
	}
	
	# total igm opacity
	$ctau = 0.0;
	for ($j=1; $j<=40; $j++) {
	    $ctau += $tau[$j];
	}
	
	$trans = exp(-$ctau);
	push @trans, $trans;
    }
    $trans = pdl @trans;
    return  $trans;
}

sub gammln($) {
  my ($x, $tmp, $ser, $j);
#double x,tmp,ser;
  my @cof;
#       static double cof[6]={76.18009173,-86.50532033,24.01409822,
#               -1.231739516,0.120858003e-2,-0.536382e-5};
#       int j;
  $cof[0]= 76.18009173;
  $cof[1] =-86.50532033;
  $cof[2] = 24.01409822;
  $cof[3] = -1.231739516;
  $cof[4] = 0.120858003e-2;
  $cof[5] = -0.536382e-5;
#       x=xx-1.0;
#       tmp=x+5.5;
  $xx = $_[0];
  $x = $xx - 1.0;
  $tmp = $x + 5.5;
  $tmp -= ($x+0.5)*log($tmp);
  $ser = 1.0;
#       tmp -= (x+0.5)*log(tmp);
#       ser=1.0;
  for ($j = 0; $j <= 5; $j++) {
    $x += 1.0;
    $ser += $cof[$j]/$x;
  }
#       for (j=0;j<=5;j++) {
       #       x += 1.0;
       #       ser += cof[j]/x;
#       }
#       return -tmp+log(2.50662827465*ser);

  return -$tmp+log(2.50662827465*$ser);
}

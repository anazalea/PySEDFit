#! /usr/bin/perl
#
# Program sed_fit.pl for maximum-likelihood fitting of galaxy data to spectral energy distribution models. 
#

use PDL;
use PGPLOT;
use PDL::Graphics::PGPLOT;
use PDL::AutoLoader;
use PDL::Math;
use PDL::Primitive;
use PDL::Astro::Cosmology;

$version = 20090802; 
print "sed_fit version $version \n"; 

$TINY = 1e-30;
$HUGE = 1e+30;

# READ IN AND PARSE THE PARAMETERS THAT WILL CONTROL THE PROCEDURE


$paramfile =@ARGV[0];
#print "Parsing parameter file $paramfile ... \n";

open (PARAMETERS , $paramfile) or die "USAGE: sed_fit parameterfile [-PARAMETER1 value[s] -PARAMETER2 value[s]...]\n";
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
  if (((substr ($argv,0,1)) =~ "-") & (substr ($argv,1,1) !~ '[0123456789]')){ # if first character of this string is a "-" and the 2nd character is NOT a number then this string is a parameter name
    $lp++; 
    @parameter[$lp] = "";
  }
  @parameter[$lp] = @parameter[$lp] . $argv. " ";

}


$parameters = "";
$Nrestricted_parameters=0;
$Nlp = $lp;
for ($lp=0; $lp<=$Nlp; $lp++){
  $parameters = $parameters . " " . @parameter[$lp];
  $line = @parameter[$lp]; 
  (@objects) = split ' ', $line;
  $parameter = uc ( shift @objects);
#  print @objects; print "\n";

  if ($parameter =~ "MODEL_DIR")            {@MODEL_DIR = @objects;}
  if ($parameter =~ "MODEL_BBSED_FILE")     {@MODEL_BBSED_FILE = @objects;}
  if ($parameter =~ "MODEL_MAGS")           {@MODEL_MAGS = @objects;}
  if ($parameter =~ "MODEL_PARAM")          {@MODEL_PARAMETER = @objects;					     
                                             if (scalar @objects <=3){$temp = "no";}
					     else {$temp = pop @objects;}
					     push @MODEL_PARAM_FOURTHCOL, $temp;

    					     $temp = pop @objects; push @MODEL_PARAM_OUTPUTFORMAT, $temp;
					     $temp = pop @objects; push @MODEL_PARAM_NAME, $temp;
					     $temp = pop @objects; push @MODEL_PARAM_COLUMN, $temp;
					     }
  if ($parameter =~ "RESTRICT_PARAM")       {@RESTRICT_PARAMETER = @objects;
					     $nelem = scalar @RESTRICT_PARAMETER;
					     if (($nelem<3) | ($nelem >4)){
					       print "ERROR: a RESTRICT_PARAMETER not specified properly\n"; exit;
					     }
					     push @restrict_parameter_column, $RESTRICT_PARAMETER[0];
					     $restrict_parameter_style = lc($RESTRICT_PARAMETER[1]);
					     if (($restrict_parameter_style =~ 'range') & ($restrict_parameter_style =~ 'close')){
					       print "ERROR: a RESTRICT_PARAMETER not specified properly\n"; exit;}
					     push @restrict_parameter_style, $restrict_parameter_style;
					     $temp3 = @RESTRICT_PARAMETER[2];
					     if ($nelem == 4){$temp4 = @RESTRICT_PARAMETER[3];} #if there are 4 parameter, read in the 4th
					     else {$temp4 = $temp3;}                            #otherwise, assign the value of the third to the 4th
					     ($restrict_parameter_minvalue, $restrict_parameter_maxvalue) = minmax (append($temp3, $temp4)); # find which one is min and which one is max
					     push @restrict_parameter_minvalue, $restrict_parameter_minvalue;
					     push @restrict_parameter_maxvalue, $restrict_parameter_maxvalue;
					   }
				
###
  if ($parameter =~ "OUTPUT_FILE")          {@OUTPUT_FILE = @objects;}
  if ($parameter =~ "OUTPUT_OVERWRITE")     {@OUTPUT_OVERWRITE = @objects;}

  if ($parameter =~ "BESTFIT_SPECTRA_YN")   {@BESTFIT_SPECTRA_YN = @objects;}
  if ($parameter =~ "BESTFIT_SPECTRA_MKFILE"){@BESTFIT_SPECTRA_MKFILE = @objects;}

  if ($parameter =~ "CHISQ_MATRIX_YN")      {@CHISQ_MATRIX_YN = @objects;}
  if ($parameter =~ "CHISQ_MATRIX_DIR")     {@CHISQ_MATRIX_DIR = @objects;}

  if ($parameter =~ "FIT_ERRORBARS")        {@FIT_ERRORBARS = @objects;}
  if ($parameter =~ "SAVE_MC_RESULTS_YN")   {@SAVE_MC_RESULTS_YN = @objects;}

  if ($parameter =~ "VERBOSE")              {@VERBOSE = @objects;}

###
  if ($parameter =~ "DATA_FILE")            {@DATA_FILE = @objects;}
  if ($parameter =~ "DATA_DIR")             {@DATA_DIR = @objects;}
  if ($parameter =~ "DATA_PARAM")           {
                                             if (scalar @objects <=3){$temp = "no";}
					     else {$temp = pop @objects;}
					     push @DATA_PARAM_FOURTHCOL, $temp;

                                             $temp = pop @objects; push @DATA_PARAM_OUTPUTFORMAT, $temp;
					     $temp = pop @objects; push @DATA_PARAM_NAME, $temp;
					     $temp = pop @objects; push @DATA_PARAM_COLUMN, $temp;}
  if ($parameter =~ "DATA_MAGS")            {@DATA_MAGS = @objects;}
  if ($parameter =~ "DATA_UNCERTAINTIES")   {@DATA_UNCERTAINTIES = @objects;}


  if ($parameter =~ "DATA_UPPERLIM_YN")        {@DATA_UPPERLIM_YN = @objects;}
  if ($parameter =~ "DATA_UPPERLIM_FLAG_COLS") {@DATA_UPPERLIM_FLAG_COLS = @objects;}
  if ($parameter =~ "DATA_UPPERLIM_FLAG_VALS") {@DATA_UPPERLIM_FLAG_VALS = @objects;}
  if ($parameter =~ "DATA_UPPERLIM_FLAG_OPER") {@DATA_UPPERLIM_FLAG_OPER = @objects;}
  if ($parameter =~ "DATA_UPPERLIM_LIM_COLS")  {@DATA_UPPERLIM_LIM_COLS  = @objects;}
  if ($parameter =~ "DATA_UPPERLIM_NSIGMA")    {@DATA_UPPERLIM_NSIGMA    = @objects;}
  

  if ($parameter =~ "DATA_MAG_OFFSETS")     {@DATA_MAG_OFFSETS   = @objects;}
  if ($parameter =~ "DATA_WAVELENGTHS")     {@DATA_WAVELENGTHS   = @objects;}
###
  if ($parameter =~ "MAG_SOFTENING")        {@MAG_SOFTENING = @objects;}
  if ($parameter =~ "FITMETHOD_DETECTED")   {@FITMETHOD_DETECTED = @objects;}
  if ($parameter =~ "FITMETHOD_UPPERLIM")   {@FITMETHOD_UPPERLIM = @objects;}

###
  if ($parameter =~ "PLOTTING_SPECFILE")    {@PLOTTING_SPECFILE = @objects;}
}




if (scalar @MODEL_BBSED_FILE == 0)          { print "ERROR: No bbsed model specified\n";exit;}
if (scalar @MODEL_MAGS != scalar @DATA_MAGS){ print "ERROR: Number of data filters does not match number of model filters\n";exit;}
if (scalar @DATA_MAGS != scalar @DATA_UNCERTAINTIES){ print "ERROR: Number of data filters does not match number of data uncertainties\n";exit;}
if ((scalar @DATA_MAG_OFFSETS > 0) & (scalar @DATA_MAG_OFFSETS != scalar @DATA_MAGS))
  { print "ERROR: you have requested mag zeropoint offsets, but the number supplied does not match the number of data mags\n";exit;}

$MODEL_DIR = @MODEL_DIR[0];
$MODEL_BBSED_FILE = @MODEL_BBSED_FILE[0];
$DATA_DIR =  @DATA_DIR[0];
$DATA_FILE = @DATA_FILE[0];


if (scalar @OUTPUT_FILE != 0){
  $OUTPUT_FILE = @OUTPUT_FILE[0];
  $OUTPUT_OVERWRITE = @OUTPUT_OVERWRITE[0];
  if (-e ($OUTPUT_FILE)){
    if (lc($OUTPUT_OVERWRITE) !~ 'yes'){
      print "ERROR:  Specified output file already exists.  Set OUTPUT_OVERWRITE to 'yes' if you want me to overwrite it. \n";
      exit; }
  }
  open (OUTPUT_FILE , "> $OUTPUT_FILE");
}

if (@VERBOSE[0] =~ 'y'){
  $VERBOSE=1;
  print "Verbose mode on.\n";
}
else {$VERBOSE=0;}


$BESTFIT_SPECTRA_YN = @BESTFIT_SPECTRA_YN[0];
if ($BESTFIT_SPECTRA_YN =~ 'yes'){
  $BESTFIT_SPECTRA_MKFILE = @BESTFIT_SPECTRA_MKFILE[0];
  open (BESTFIT_SPECTRA_MKFILE, "> $BESTFIT_SPECTRA_MKFILE");
}

$CHISQ_MATRIX_YN = @CHISQ_MATRIX_YN[0];
$CHISQ_MATRIX_DIR = @CHISQ_MATRIX_DIR[0];

$ERRORBARS_TYPE= @FIT_ERRORBARS[0];
if (length $ERRORBARS_TYPE <= 0){$ERRORBARS_TYPE='none';}
if (lc($ERRORBARS_TYPE) =~ 'dchisq'){
  $Dchisq_level = @FIT_ERRORBARS[1];
  if ($Dchisq_level <= 0){print "ERROR: Must specify Dchisq that's > 0 in parameter FIT_ERRORBARS\n"; exit;}
}
elsif (lc($ERRORBARS_TYPE) =~ 'frac') {
  $ERRORBARS_LEVEL = @FIT_ERRORBARS[1];
  $MC_PHOT_NITER = @FIT_ERRORBARS[2];
  if (($ERRORBARS_LEVEL <= 0) | ($ERRORBARS_LEVEL >= 1)){print "ERROR: Must specify Dchisq that's between 0 and 1 in parameter FIT_ERRORBARS\n"; exit;}
  if ($MC_PHOT_NITER <= 1) {print "ERROR: Must specify number of MC iterations that's greater than 1 in parameter FIT_ERRORBARS\n"; exit;}
  if ($MC_PHOT_NITER <= 50) {print "WARNING: You have specified $MC_PHOT_NITER Monte Carlo iterations in parameter FIT_ERRORBARS for estimating parameter uncertainties.  Recommend that you use a larger number, but proceeding with the fit anyway. \n";}
}

$SAVE_MC_RESULTS_YN = @SAVE_MC_RESULTS_YN[0];
if (lc($ERRORBARS_TYPE) !~ 'frac') {$SAVE_MC_RESULTS_YN = 'no';}  # no point saving the mc results if we are not doing a mc simulation!

#open the output file for storing MC iteration results
if (lc ($SAVE_MC_RESULTS_YN) !~ 'n'){
  $OUTPUT_MC_FILE = $OUTPUT_FILE . ".mc";
  open (OUTPUT_MC_FILE , "> $OUTPUT_MC_FILE");
}


# The three separate parsers (DATA_UPPERLIM_FLAG, DATA_UPPERLIM_NSIGMA, MAG_SOFTENING) that follow could potentially be replaced by one function that handles this kind of parsing of single vs. multiple parameter values



$DATA_UPPERLIM_YN = @DATA_UPPERLIM_YN[0];
if ($DATA_UPPERLIM_YN =~ 'yes'){
  # the flux-upper-limit parsers:
  if (scalar @DATA_UPPERLIM_FLAG_COLS == 0){
    @DATA_UPPERLIM_FLAG_COLS = @DATA_MAGS;  # this sets the default columns for upper limit flags
  }
  elsif (scalar @DATA_UPPERLIM_FLAG_COLS > 0){
    if ((scalar @DATA_UPPERLIM_FLAG_COLS) != (scalar @DATA_MAGS)){
      print "ERROR: the specified number of DATA_UPPERLIM_FLAG_COLS is inconsistent with that expected from the number of DATA_MAGS. \n"; exit;
    }
  }
  
  if (scalar @DATA_UPPERLIM_FLAG_VALS >0){
    if (scalar @DATA_UPPERLIM_FLAG_VALS == 1){
      $data_upperlim_flag_vals = @DATA_UPPERLIM_FLAG_VALS[0];
      $data_upperlim_flag_vals =  $data_upperlim_flag_vals * ones(scalar @DATA_MAGS);
    }
    elsif (scalar @DATA_UPPERLIM_FLAG_VALS != scalar @DATA_MAGS){
      print "ERROR: Number of upper limit flag values specified in DATA_UPPERLIM_FLAG_VALS is > 1 but does not match the number of data columns specified in DATA_MAGS\n"; exit;}
    else{
      $data_upperlim_flag_vals = pdl(0);
      foreach $datamag (@DATA_UPPERLIM_FLAG_VALS){
	$data_upperlim_flag_vals = append ($data_upperlim_flag_vals, $datamag);
      }
      $nelem = nelem $data_upperlim_flag_vals;
      $slice = '1:' . ($nelem-1);
      $data_upperlim_flag_vals=$data_upperlim_flag_vals->slice($slice);
    }
  }

  if (scalar @DATA_UPPERLIM_FLAG_OPER == 0){
    $data_upperlim_flag_oper = "=";
  }
  elsif (scalar @DATA_UPPERLIM_FLAG_OPER == 1){
    $data_upperlim_flag_oper = @DATA_UPPERLIM_FLAG_OPER[0]; 
  }
  elsif (scalar @DATA_UPERLIM_FLAG_OPER >1) {
    print "ERROR: Number of upper limit flag values specified in DATA_UPPERLIM_FLAG_OPER is > 1 \n"; exit;
  }
  
  if (scalar @DATA_UPPERLIM_LIM_COLS == 0){
    @DATA_UPPERLIM_LIM_COLS = @DATA_UNCERTAINTIES;  # this sets the default columns for upper limit flags
  }
  elsif (scalar @DATA_UPPERLIM_LIM_COLS > 0){
    if ((scalar @DATA_UPPERLIM_LIM_COLS) != (scalar @DATA_MAGS)){
      print "ERROR: the specified number of DATA_UPPERLIM_LIM_COLS is inconsistent with that expected from the number of DATA_MAGS. \n"; exit;
    }
  }
  
  if ((scalar @DATA_UPPERLIM_FLAG_VALS > 0) & (scalar @DATA_UPPERLIM_NSIGMA == 0 )){
    print "ERROR: upper limit fitting requested but no sigma value specified in DATA_UPPERLIM_NSIGMA. \n"; exit;}
  elsif (scalar @DATA_UPPERLIM_NSIGMA >0){
    if (scalar @DATA_UPPERLIM_NSIGMA == 1){
      $data_upperlim_nsigma = @DATA_UPPERLIM_NSIGMA[0];
      $data_upperlim_nsigma =  $data_upperlim_nsigma * ones(scalar @DATA_MAGS);
    }
    elsif (scalar @DATA_UPPERLIM_NSIGMA != scalar @DATA_MAGS){
      print "ERROR: Number of upper limit flag values specified in DATA_UPPERLIM_NSIGMA is > 1 but does not match the number of data columns specified in DATA_MAGS\n"; exit;}
    else{
      $data_upperlim_nsigma = pdl(0);
      foreach $datamag (@DATA_UPPERLIM_NSIGMA){
	$data_upperlim_nsigma = append ($data_upperlim_nsigma, $datamag);
      }
      $nelem = nelem $data_upperlim_nsigma;
      $slice = '1:' . ($nelem-1);
      $data_upperlim_nsigma=$data_upperlim_nsigma->slice($slice);
    }
  }
}
  
  
# end of flux-upper-limits parsers 

#print @DATA_UPPERLIM_FLAG_COLS, "\n";
#print @DATA_UPPERLIM_LIM_COLS, "\n";
#print @DATA_UPPERLIM_FLAG_OPER, "\n";
#print $data_upperlim_flag_vals; 
#print $data_upperlim_nsigma; 




#### XXX




# the mag-softening parser:
if (scalar @MAG_SOFTENING > 0){
  if (scalar @MAG_SOFTENING == 1){
    $mag_softening = @MAG_SOFTENING[0];
    $mag_softening = $mag_softening * ones(scalar @DATA_MAGS);
  }
  elsif(scalar @MAG_SOFTENING != scalar @DATA_MAGS){
    print "ERROR: Number of magnitude softening parameters specified in MAG_SOFTENING is > 1 but does not match the number of data columns specified in DATA_MAGS\n"; exit;}
  else{
    $mag_softening = pdl(0);
    foreach $datamag (@MAG_SOFTENING){
      $mag_softening = append ($mag_softening, $datamag);
    }
    $nelem = nelem $mag_softening;
    $slice = '1:' . ($nelem-1);
    $mag_softening=$mag_softening->slice($slice);
  }
}

# the fitting method parser:
#    for detected objects:
if (scalar @FITMETHOD_DETECTED == 0) {$FITMETHOD_DETECTED = 'brute';} #default to brute force if no method is specified
elsif (scalar @FITMETHOD_DETECTED > 0) {
  $FITMETHOD_DETECTED = lc @FITMETHOD_DETECTED[0];
  if ($FITMETHOD_DETECTED =~ 'brute'){
    $FITMETHOD_DETECTED = 'brute';}
  elsif ($FITMETHOD_DETECTED =~ 'downhill'){
    $FITMETHOD_DETECTED = 'downhill';}
  elsif ($FITMETHOD_DETECTED =~ 'anneal'){
    $FITMETHOD_DETECTED = 'brute';
    print "WARNING: Simulated annealing search not implemented yet for FITMETHOD_DETECTED; defaulting to brute force search\n";
  }
  else {
    $FITMETHOD_DETECTED = 'brute';
    print "WARNING: unrecognized FITMETHOD_DETECTED was specified ; defaulting to brute force search\n";}
}
if ($VERBOSE>0){print "Fitting method for detected objects: $FITMETHOD_DETECTED \n";}

#    for undetected objects:
if (scalar @FITMETHOD_UPPERLIM == 0) {$FITMETHOD_UPPERLIM = 'brute';} #default to brute force if no method is specified
elsif (scalar @FITMETHOD_UPPERLIM > 0) {
  $FITMETHOD_UPPERLIM = lc @FITMETHOD_UPPERLIM[0];
  if (scalar @FITMETHOD_UPPERLIM > 1) {
    $FITMETHOD_UPPERLIM_OPTION1 = lc  @FITMETHOD_UPPERLIM[1];
  }
  if ($FITMETHOD_UPPERLIM =~ 'brute'){
    $FITMETHOD_UPPERLIM = 'brute';}
  elsif ($FITMETHOD_UPPERLIM =~ 'downhill'){
    $FITMETHOD_UPPERLIM = 'downhill';
    $Ndownhill_repeats = 1; 
    if ($FITMETHOD_UPPERLIM_OPTION1 >= 1){
      $Ndownhill_repeats = floor $FITMETHOD_UPPERLIM_OPTION1; 
    }
  }
  elsif ($FITMETHOD_UPPERLIM =~ 'anneal'){
    $FITMETHOD_UPPERLIM = 'brute';
    print "WARNING: Simulated annealing search not implemented yet for FITMETHOD_UPPERLIM; defaulting to brute force search\n";
  }
  else {
    $FITMETHOD_UPPERLIM = 'brute';
    print "WARNING: unrecognized FITMETHOD_UPPERLIM was specified ; defaulting to brute force search\n";}
}
if ($VERBOSE>0){print "Fitting method for undetected objects: $FITMETHOD_UPPERLIM \n";}





if (scalar @RESTRICT_PARAMETER > 0){
  $restrict_parameter_column =  pdl @restrict_parameter_column;
  $restrict_parameter_minvalue =  pdl @restrict_parameter_minvalue;
  $restrict_parameter_maxvalue =  pdl @restrict_parameter_maxvalue;
}


$Nmodel_params = scalar @MODEL_PARAM_COLUMN;
$Ndata_params = scalar @DATA_PARAM_COLUMN;

##  READ IN DATA AND MAKE INTO PERL-FRIENDLY PIDDLES:

#print "### Reading data array $DATA_FILE ...\n";
@dataArray = rcols $DATA_DIR . $DATA_FILE;
$Ndata = nelem @dataArray[0];

# slice out the object parameters:
if ($Ndata_params >0){
  $param = pdl(0);
  foreach $data_param (@DATA_PARAM_COLUMN){
    $param = append ($param, @dataArray[$data_param]);
  }
  $nelem = nelem $param;
  $slice = '1:' . ($nelem-1);
  $data_param=$param->slice($slice);
}
# this lists which column in $model_param contains which column of the original bbsed file:
$data_param_column =  pdl(@DATA_PARAM_COLUMN);
## bale out if some of the parameters are read in more than once: 
#if ((nelem $data_param_column) > (nelem uniq $data_param_column)){
#  print "ERROR: some DATA_PARAM parameter columns are specified more than once.  This can cause problems. Please edit your DATA_PARAMs\n"; exit;}



# slice out the requested filters:
$md_all = pdl(0);
foreach $datamag (@DATA_MAGS){
  $md_all = append ($md_all, @dataArray[$datamag]);
}
$mag_sigma_all = pdl(0);
foreach $datamag (@DATA_UNCERTAINTIES){
  $mag_sigma_all = append ($mag_sigma_all, @dataArray[$datamag]);
}

# get the data upper limit flags: 
if ($DATA_UPPERLIM_YN =~ 'yes'){
  $data_upperlim_flag_all = pdl(0);
  foreach $datamag (@DATA_UPPERLIM_FLAG_COLS){
    $data_upperlim_flag_all = append ($data_upperlim_flag_all, @dataArray[$datamag]);
  }
  $data_upperlim_lim_all = pdl(0);
  foreach $datamag (@DATA_UPPERLIM_LIM_COLS){
    $data_upperlim_lim_all = append ($data_upperlim_lim_all, @dataArray[$datamag]);
  }
}
else {
  $data_upperlim_flag_all = zeroes $md_all;
  $data_upperlim_lim_all = zeroes $md_all; 
}

$nelem = nelem $md_all;
$slice = '1:' . ($nelem-1);
$md_all=$md_all->slice($slice);
$mag_sigma_all=$mag_sigma_all->slice($slice);
$data_upperlim_flag_all=$data_upperlim_flag_all->slice($slice);
$data_upperlim_lim_all=$data_upperlim_lim_all->slice($slice); 



# reshape magnitudes into 2d piddles:
reshape $md_all, $Ndata, (nelem $md_all)/$Ndata;
reshape $mag_sigma_all, $Ndata, (nelem $md_all)/$Ndata;
reshape $data_param, $Ndata, (nelem $data_param)/$Ndata;
reshape $data_upperlim_flag_all, $Ndata, (nelem $md_all)/$Ndata;
reshape $data_upperlim_lim_all, $Ndata, (nelem $md_all)/$Ndata;
$data_upperlim_nsgima = outer ((ones ((nelem $md_all)/(nelem $data_upperlim_nsigma))), $data_upperlim_nsigma);



# make a matrix that identifies which elements in the data matrix are upper limits (1=upper limit, 0=not)

if ($DATA_UPPERLIM_YN =~ 'yes'){
  $data_upperlim_flag_vals = outer (ones ($Ndata), $data_upperlim_flag_vals); 
  $data_upperlim_flag_yn = ($data_upperlim_flag_all -$data_upperlim_flag_vals);
  if ($data_upperlim_flag_oper =~'>='){$data_upperlim_flag_yn = $data_upperlim_flag_yn >= 0;}
  elsif ($data_upperlim_flag_oper =~'<='){$data_upperlim_flag_yn = $data_upperlim_flag_yn <= 0;}
  elsif ($data_upperlim_flag_oper =~ '='){$data_upperlim_flag_yn = $data_upperlim_flag_yn == 0;}
  elsif ($data_upperlim_flag_oper =~ '>'){$data_upperlim_flag_yn = $data_upperlim_flag_yn > 0;}
  elsif ($data_upperlim_flag_oper =~ '<'){$data_upperlim_flag_yn = $data_upperlim_flag_yn < 0;}
  else {print "ERROR: The operator specified in DATA_UPPERLIM_FLAG_OPER is not valid\n"; exit;}
}
else {
  $data_upperlim_flag_yn = zeroes $md_all;
}
### XXXX I am here in my development of UL fitting:
#print $data_upperlim_flag_yn;  # <- this variable contains the flag that identifies upper limits in the data
#print $md_all; 


# add magnitude zeropoint offsets (if requested):
if (scalar @DATA_MAG_OFFSETS > 0){
  $data_mag_offsets = @DATA_MAG_OFFSETS[0];
  for ($datamag=1; $datamag < (scalar @DATA_MAG_OFFSETS); $datamag++){
    $data_mag_offsets = append ($data_mag_offsets, @DATA_MAG_OFFSETS[$datamag]);
  }
  $data_mag_offsets =  outer (ones((nelem $md_all)/ (nelem $data_mag_offsets)), $data_mag_offsets);
  $md_all = $md_all + $data_mag_offsets; # <- adding the offsets to the magnitudes
  $data_upperlim_lim_all = $data_upperlim_lim_all + $data_mag_offsets # <- adding the offsets to the upper limits
}


# convert to fluxes:
$fd_all = 10**(-0.4*$md_all);   # <- this is the observed flux (in flux units, of course)
$fd_upperlim_lim_all = 10**(-0.4*$data_upperlim_lim_all);   # <- this is the upper limit in flux units

# figure out what is 1 sigma in the upper limit case:
$data_upperlim_nsigma_all = outer (ones($Ndata),$data_upperlim_nsigma);
$upperlim_sigma_all = $fd_upperlim_lim_all / $data_upperlim_nsigma_all;   #<- this stores, in flux units, the 1 sigma value of the upper limit


# add softening to the uncertainties: 
$mag_softening = outer ((ones ($Ndata)), $mag_softening);
$mag_sigma_all = (($mag_sigma_all**2 + $mag_softening**2)**0.5);

##  READ IN MODELS AND MAKE INTO PERL-FRIENDLY PIDDLES:

#print "### Reading model array $MODEL_BBSED_FILE ...\n";

$modelfile = $MODEL_DIR . $MODEL_BBSED_FILE;

open (MODELFILE, $modelfile);
$linecount=0;
while (($line = <MODELFILE>) && ($linecount<1000)){
  $linecount++;
  if ($line =~ "# full params list"){
    ($junk, $BBSEDFILE_PARAMS) = split '=', $line;
    chomp $BBSEDFILE_PARAMS;
  }
}

close (MODELFILE);

@modelArray = rcols $MODEL_DIR . $MODEL_BBSED_FILE;
$Nmodels = nelem @modelArray[0];


$mm = pdl(0);
$Nmodel_mags = scalar @MODEL_MAGS;
foreach $modelmag (@MODEL_MAGS){
  $mm = append ($mm, @modelArray[$modelmag]);
}
$nelem = nelem $mm;
$slice = '1:' . ($nelem-1);
$mm=$mm->slice($slice);


# slice out the requested model parameters
if ($Nmodel_params > 0){
  $param = pdl(0);
  foreach $model_param (@MODEL_PARAM_COLUMN){
    $param = append ($param, @modelArray[$model_param]);
  }
  $nelem = nelem $param;
  $slice = '1:' . ($nelem-1);
  $model_param=$param->slice($slice);
}
# this lists which column in $model_param contains which column of the original bbsed file:
$model_param_column =  pdl(@MODEL_PARAM_COLUMN);
# bale out if some of the parameters are read in more than once: 
if ((nelem $model_param_column) > (nelem uniq $model_param_column)){
  print "ERROR: some MODEL_PARAM parameter columns are specified more than once.  This can cause problems. Please edit your MODEL_PARAMs\n"; exit;}

# reshape magnitudes into 2d piddles:
reshape $mm, $Nmodels, (nelem $mm)/$Nmodels;
reshape $model_param, $Nmodels, (nelem $model_param)/$Nmodels;

#convert to fluxes;
$fm = 10**(-0.4*$mm);
$fm = $fm->xchg(0,1);
$model_param = $model_param->xchg(0,1);

#and save the original for future use and reference:
$fm_orig = $fm;
$model_param_orig = $model_param;




## PRINT THE HEADER:
print OUTPUT_FILE "# parameter file =    $paramfile \n";
print OUTPUT_FILE "# model =             $MODEL_BBSED_FILE \n";
print OUTPUT_FILE "# data  =             $DATA_FILE \n";
if ((scalar @FIT_ERRORBARS > 0) & (lc(@FIT_ERRORBARS[0]) !~ 'none')){
  print OUTPUT_FILE "# uncertainties =     @FIT_ERRORBARS[0] = @FIT_ERRORBARS[1] \n";}
if (scalar @DATA_MAG_OFFSETS > 0){
  print OUTPUT_FILE "# data mag offsets =  @DATA_MAG_OFFSETS \n";}
#print OUTPUT_FILE "# softening (mags) =  $mag_softening \n";
print OUTPUT_FILE "# full params list =  $parameters\n";
print OUTPUT_FILE "#\n";
print OUTPUT_FILE "# Columns:\n";
# print the data parameter column names:
$columnCounter = 0;
for ($parameter = 0; $parameter < $Ndata_params; $parameter++){
  print OUTPUT_FILE "#   (", $columnCounter, ") ".@DATA_PARAM_NAME[$parameter]." \n";  $columnCounter++;
}
# print the best-fit model column names:
for ($parameter = 0; $parameter < $Nmodel_params; $parameter++){
  print OUTPUT_FILE "#   (", $columnCounter, ") ".@MODEL_PARAM_NAME[$parameter]." \n"; $columnCounter++;
  if (lc($ERRORBARS_TYPE) !~ 'none'){
    print OUTPUT_FILE "#   (", $columnCounter, ")     lower allowable value of ".@MODEL_PARAM_NAME[$parameter]." \n"; $columnCounter++;
    print OUTPUT_FILE "#   (", $columnCounter, ")     upper allowable value of ".@MODEL_PARAM_NAME[$parameter]." \n"; $columnCounter++;
  }
}
# print the column names for the chi^2 and flux normalization columns:
print OUTPUT_FILE "#   (", $columnCounter, ") flux normalization \n"; $columnCounter++;
if (lc($ERRORBARS_TYPE) !~ 'none'){
  print OUTPUT_FILE "#   (", $columnCounter, ")     lower allowable value of flux normalization \n"; $columnCounter++;
  print OUTPUT_FILE "#   (", $columnCounter, ")     upper allowable value of flux normalization \n"; $columnCounter++;
}
print OUTPUT_FILE "#   (", $columnCounter, ") chi^2 \n"; $columnCounter++;
if (lc($ERRORBARS_TYPE) !~ 'none'){
  print OUTPUT_FILE "#   (", $columnCounter, ") Delta_chi^2 for errorbars \n"; $columnCounter++;
}



## AND PRINT THE HEADER FOR THE MC OUTPUT FILE:

print OUTPUT_MC_FILE "# parameter file =    $paramfile \n";
print OUTPUT_MC_FILE "# model =             $MODEL_BBSED_FILE \n";
print OUTPUT_MC_FILE "# data  =             $DATA_FILE \n";
if ((scalar @FIT_ERRORBARS > 0) & (lc(@FIT_ERRORBARS[0]) !~ 'none')){
  print OUTPUT_MC_FILE "# uncertainties =     @FIT_ERRORBARS[0] = @FIT_ERRORBARS[1] \n";}
if (scalar @DATA_MAG_OFFSETS > 0){
  print OUTPUT_MC_FILE "# data mag offsets =  @DATA_MAG_OFFSETS \n";}
#print OUTPUT_MC_FILE "# softening (mags) =  $mag_softening \n";
print OUTPUT_MC_FILE "# full params list =  $parameters\n";
print OUTPUT_MC_FILE "#\n";
print OUTPUT_MC_FILE "# Columns:\n";
# print the data parameter column names:
$columnCounter = 0;
for ($parameter = 0; $parameter < $Ndata_params; $parameter++){
  print OUTPUT_MC_FILE "#   (", $columnCounter, ") ".@DATA_PARAM_NAME[$parameter]." \n";  $columnCounter++;
}
# print the best-fit model column names:
for ($parameter = 0; $parameter < $Nmodel_params; $parameter++){
  print OUTPUT_MC_FILE "#   (", $columnCounter, ") ".@MODEL_PARAM_NAME[$parameter]." \n"; $columnCounter++;
  if (lc($ERRORBARS_TYPE) !~ 'none'){
    print OUTPUT_MC_FILE "#   (", $columnCounter, ")     upper allowable value of ".@MODEL_PARAM_NAME[$parameter]." \n"; $columnCounter++;
    print OUTPUT_MC_FILE "#   (", $columnCounter, ")     lower allowable value of ".@MODEL_PARAM_NAME[$parameter]." \n"; $columnCounter++;
  }
}
# print the column names for the chi^2 and flux normalization columns:
print OUTPUT_MC_FILE "#   (", $columnCounter, ") flux normalization \n"; $columnCounter++;
if (lc($ERRORBARS_TYPE) !~ 'none'){
  print OUTPUT_MC_FILE "#   (", $columnCounter, ")     upper allowable value of flux normalization \n"; $columnCounter++;
  print OUTPUT_MC_FILE "#   (", $columnCounter, ")     lower allowable value of flux normalization \n"; $columnCounter++;
}
print OUTPUT_MC_FILE "#   (", $columnCounter, ") chi^2 \n"; $columnCounter++;
if (lc($ERRORBARS_TYPE) !~ 'none'){
  print OUTPUT_MC_FILE "#   (", $columnCounter, ") Delta_chi^2 for errorbars \n"; $columnCounter++;
}




## AND NOW THAT ALL THE PRELIMINARIES ARE TAKEN CARE OF, DO THE FITTING:


#print "### Fitting..."; print " \n";
for ($object=0; $object < $Ndata; $object++){
  if ($VERBOSE>0){print "\n\nFitting object $object \n";}
  # copy the model data into working arrays:
  $fm = $fm_orig;
  $model_param = $model_param_orig;

  # and now (if requested) select only those models that have been requested:
  $Nrestrict_parameters = scalar @RESTRICT_PARAMETER;
  if ($Nrestrict_parameters> 0){
    for ($restrict_parameter = 0; $restrict_parameter<$Nrestrict_parameters; $restrict_parameter++){

      if ((@restrict_parameter_style[$restrict_parameter] =~ 'range')
	  | (@restrict_parameter_style[$restrict_parameter] =~ 'closeval')
	  | (@restrict_parameter_style[$restrict_parameter] =~ 'closecol')){

	$this_restrict_parameter_column = $restrict_parameter_column->index($restrict_parameter);
	$this_restrict_parameter_minvalue = $restrict_parameter_minvalue->index($restrict_parameter);
	$this_restrict_parameter_maxvalue = $restrict_parameter_maxvalue->index($restrict_parameter);

	# this is the parameter we will be working on: 
	$this_column = which ($model_param_column == $this_restrict_parameter_column);
	$this_param = $model_param->index($this_column);



	# if the parameter matching is requested as 'close', then do the following.  This finds the model parameter
	#  that's closest to the requested value and then sets the min and max allowed values to that value:

	if (@restrict_parameter_style[$restrict_parameter] =~ 'close'){

	  # what is the requested value?:
	
	  #if we are being asked for the closest match to a specified value:
	  if (@restrict_parameter_style[$restrict_parameter] =~ 'closeval'){
	    $this_restrict_parameter_closest_requestedvalue =  $this_restrict_parameter_maxvalue;
	  }
	  # else if, instead,  we are being asked for a closest match to the value in a specified column:
	  elsif (@restrict_parameter_style[$restrict_parameter] =~ 'closecol'){
	    $this_restrict_parameter_closest_requestedvalue =  $this_restrict_parameter_maxvalue;
#	     print $this_restrict_parameter_closest_requestedvalue; print "\n";
#	     print $data_param_column->index(31); print "\n";
	    $this_data_param_column =  ($data_param_column->index($this_restrict_parameter_closest_requestedvalue));
	    $this_restrict_parameter_closest_requestedvalue = (($data_param->index($object))->index($this_data_param_column));
	  }

	

	  # compile the unique values of the parameter of interest:
	  $this_parameter_uniq_values = uniq $this_param;

	  # find the closest match between the model parameters and the requested parameter:
	  $this_restrict_parameter_closest_requestedvalue;
	  $minimum_ind = minimum_ind(abs ($this_parameter_uniq_values - $this_restrict_parameter_closest_requestedvalue));
	  $this_restrict_parameter_closest_modelvalue = $this_parameter_uniq_values->index($minimum_ind);

	  # set the max and min allowed values to that closest model value:
	  $this_restrict_parameter_minvalue = $this_restrict_parameter_closest_modelvalue;
	  $this_restrict_parameter_maxvalue = $this_restrict_parameter_closest_modelvalue;
	}


	# note that we must do everything on the $model_param array separately from the $fm array because 
	#   they have differnet number of elements. So here we go...

	# do some prep work with arrays:
	$ones_param = ones ((nelem $model_param) / (nelem $this_param));
	$this_param_param = outer $ones_param, $this_param;
	$ones_fm = ones ((nelem $fm) / (nelem $this_param));
	$this_param_fm = outer $ones_fm, $this_param;

	#      print $this_param_param;

	# select the models within the right range of requested restricted parameter values
	$model_param  = where ($model_param,
			       ($this_param_param <= $this_restrict_parameter_maxvalue) &
			       ($this_param_param >= $this_restrict_parameter_minvalue));
	$fm           = where ($fm,
			       ($this_param_fm    <= $this_restrict_parameter_maxvalue) &
			       ($this_param_fm    >= $this_restrict_parameter_minvalue));

	# reshape the arrays back into the right format:
	reshape $model_param, $Nmodel_params, (nelem $model_param)/$Nmodel_params;
	reshape $fm, $Nmodel_mags, (nelem $fm)/$Nmodel_mags;
      }

      if ((nelem $model_param) == 0) {print "ERROR: there are no models to fit; you may have been too stringent in using RESTRICT_PARAM to restrict the allowable parameter values\n"; exit;}

    }
  }
  $Nmodels = $fm->dim(1);    # <- number of models remaining after the required ones have been sliced out
  

  $sigma_all = $fd_all*$mag_sigma_all;

# go through the Monte Carlo iterations for perturbing the model photometry (0th iteration is the actual fit using actual values): 
  $bestchisq_MC  = zeroes ($MC_PHOT_NITER+1);
  for ($MC_fit_iteration = 0; $MC_fit_iteration <= $MC_PHOT_NITER; $MC_fit_iteration++){
    if ($VERBOSE>0) {print " doing MC iteration $MC_fit_iteration"; `date`; print "\n"};

    $fd = $fd_all -> index($object);
    $fd_upperlim_lim = $fd_upperlim_lim_all -> index($object);

    $sigma = $sigma_all -> index ($object);
    $upperlim_sigma_all;
    $upperlim_sigma = $upperlim_sigma_all -> index($object);
    
    
   #perturb the model fluxes by the gaussian errorbars (assumes that the phot errors are given in sigmas in the catalog):
    if ( $MC_fit_iteration > 0 ){ $fd = $fd+$sigma*grandom(nelem($sigma)); }


# Thus end the preliminaries. Here comes the meat and guts of the fitting program:
############################################################
    # DO THE CHI-SQUARE FITTING:


#    $fd_upperlim_lim = $fd_upperlim_lim*$ones; 
#    $upperlim_sigma = $upperlim_sigma*$ones; 

    $thisobject_data_upperlim_flag_yn = $data_upperlim_flag_yn -> index($object);   # <- identifies filters that are upper lims for this object
    $thisobject_data_detection_flag_yn = abs($thisobject_data_upperlim_flag_yn-1);  # <- identifies filters that are detections for this object


# the chisq fitter for the case when all bands are detections:
    if (sumover($thisobject_data_upperlim_flag_yn) == 0){   
#      print "not an upper limit \n";

      $ones = ones $fm;
      $fd = $fd*$ones;
      $sigma = $sigma*$ones;
      $sigma2 = $sigma**2;

      # this here is the actual fitting (for the case of detections in all bands): 
      $s = ((sumover($fd*$fm/$sigma2))/(sumover($fm**2/$sigma2)));  # the scale factor
      $s = ((($ones->xchg(0,1))*$s)->xchg(0,1));
      $chisq = sumover ((($fd-$s*$fm)/$sigma)**2);
     
      $ss = medover $s; 

      $ix = qsorti($chisq);
      $bestix = $ix->at(0);
      $bestchisq = $chisq->index($bestix);
      $flux_scale=$ss->index($bestix);

    }




# the chisq fitter for the case when at least some of the bands are UPPER LIMITS:
    else {  

#      print $fd_upperlim_lim; 
#      print $upperlim_sigma;
#      print $fd; 
#      print $thisobject_data_upperlim_flag_yn;
#      print $sigma; 



# brute force method for upper limit case:
      
      if ($FITMETHOD_UPPERLIM =~ 'brute'){

	# just simply call the integrator routine that calculates chisq in the upper-limit case for all models in the model grid:

	($chisq, $ss) = find_allmodelsChisq_UL ($fm, $Nmodels, $fd, $sigma, $upperlim_sigma, 
						$fd_upperlim_lim, $thisobject_data_detection_flag_yn, 
						$thisobjec_data_upperlim_flag_yn, $fd_upperlim_lim);

	# identify the minimum chisq and the corresponding model index:

	$ix = qsorti($chisq);
        $bestix = $ix->at(0);
	$bestchisq = $chisq->index($bestix);
	$flux_scale=$ss->index($bestix);

      }


# simulated annealing OR simple downhill for upper limit case:

#      elsif ($FITMETHOD_UPPERLIM =~ 'downhill'){
      else{
	$bestchisq = $HUGE;
      	
	if ($MC_fit_iteration == 0){$Ndownhill_i = $Ndownhill_repeats;}	else {$Ndownhill_i = 1;}  # <- do downhill repeats only if this is the 0th (i.e. unperturbed) MC iteration
	if ($MC_fit_iteration > 0){ $current_model_ix = $bestix;}  # <- makes sure the starting point for the MC iterations is at the best-ever found position, not just the best position from the last downhill repeat

	for ($downhill_i=0; $downhill_i < $Ndownhill_i; $downhill_i ++){
	  if ($VERBOSE > 0) { print  " doing repeat number $downhill_i of the upper limit downhill search\n";}

	  #initialize some temp storage variables: 
	  $best_chisq_sofar = $HUGE;
	  $best_fluxscale_sofar = 0; 
	  $allmodels_ix = sequence($Nmodels);
	  $chisq = $HUGE * ones ($Nmodels); 
	  $ss = zeroes($Nmodels); 
       
	  #select a random model as starting point: 

	  if ($MC_fit_iteration == 0){ # if this is the unperturbed (0th) iteration, then start with a random model:
	    $current_model_ix = sumover (floor ($Nmodels*(random (1))));        
	  }
	  else { # else, start with the model that we found was the best model in the last iteration:
	    $current_model_ix = $current_model_ix; 
	  }


	  if ($VERBOSE>0){print "  Doing upper limit search. Trying models: \n";}

	  $Nneighbours_submodelgrid=1; # the size of subgrid (i.e., how many neighbours to include in each parameter). 
	  # loop until you find the minimum::
	  $newmodel_yn = 1;
	  while ($newmodel_yn > 0) {

	    if ($VERBOSE > 0) { print  "     $current_model_ix \n";}


	    # go through the parameters and slice out only those models that are neighbors of the current favorite model:
	    ($model_param_subgrid, $fm_subgrid, $allmodels_ix_subgrid) 
	      = select_neighbour_models ($model_param, $fm, $allmodels_ix, $Nmodel_params, $current_model_ix, $Nneighbours_submodelgrid);


	    # if we are doing a simple downhill search then get the chisq for all the models in the neighbour subgrid...	   
	    if ($FITMETHOD_UPPERLIM =~ 'downhill'){
	      $Nmodels_subgrid = (nelem $model_param_subgrid ) / $Nmodel_params;
	      # call the UL chisq fitter routine and pass it the subgrid models (rather than the full grid of models as is the case in the brute force method):
	      ($chisq_subgrid, $ss_subgrid) = find_allmodelsChisq_UL ($fm_subgrid, $Nmodels_subgrid, $fd, $sigma, $upperlim_sigma, 
								      $fd_upperlim_lim, $thisobject_data_detection_flag_yn, 
								      $thisobjec_data_upperlim_flag_yn, $fd_upperlim_lim);
	      
	      # identify the model within the subgrid that has the best-fitting chisq value:
	      $ix_subgrid = qsorti($chisq_subgrid);
	      $bestix_subgrid = $ix_subgrid->at(0);
	      $bestchisq_subgrid = $chisq_subgrid->index($bestix_subgrid);
	      $flux_scale_subgrid=$ss_subgrid->index($bestix_subgrid);
	      $bestix_allmodels_ix_subgrid = $allmodels_ix_subgrid->index($bestix_subgrid);  # <- this is the model number (in the full-model grid) that the best-fitting chisq of this SUBgrid corresponds to
	      
	      # populate the (full-grid) chisq piddle (and also the flux scale piddle) with the newly computed chisq values that were computed just for the current subgrid. 
	      #   (this is done because chances are we may return to some of these models in which case we don't want to invest the expensive time to calculate their chisq again!)
	      $chisq -> index($allmodels_ix_subgrid) .= $chisq_subgrid; 
	      $ss -> index($allmodels_ix_subgrid) .= $ss_subgrid; 
	      ##	  print "\n",$chisq; print "\n"; 
	      
	      # test to see if there is a new best-fitting model and if yes then update the current best-fit model and keep going:
	      $newmodel_yn = 0;
	      if ($bestchisq_subgrid < $best_chisq_sofar) {
		$best_chisq_sofar = $bestchisq_subgrid;
		$best_fluxscale_sofar = $flux_scale_subgrid;
		$current_model_ix = $bestix_allmodels_ix_subgrid;
		$newmodel_yn = 1;
	      }

	    }

	    # else if we are doing a simulated annealing search, then pick just one random model from the neighbours subgrid and get just its chisq:
	    elsif ($FITMETHOD_UPPERLIM =~ 'anneal'){
	      if ($VERBOSE > 0) { print  " annealing ";}

	      # go through the parameters and slice out only those models that are neighbors of the current favorite model:
	      ($model_param_subgrid, $fm_subgrid, $allmodels_ix_subgrid) 
		= select_neighbour_models ($model_param, $fm, $allmodels_ix, $Nmodel_params, $current_model_ix, $Nneighbours_submodelgrid);
	      
	      # next, pick a random model out of the subgrid. This is the model that we will be testing to see if we want to jump to it.  
	      $annealTry_model_ix = sumover (floor ($Nmodels_subgrid*(random (1))));        
	      $annalTry_ix = $allmodels_ix_subgrid->index($annealTry_model_ix);
	      
	      
	      $Nmodels_subgrid = (nelem $model_param_subgrid ) / $Nmodel_params;
	      # call the UL chisq fitter routine and pass it the subgrid models (rather than the full grid of models as is the case in the brute force method):
	      ($chisq_subgrid, $ss_subgrid) = find_allmodelsChisq_UL ($fm_subgrid, $Nmodels_subgrid, $fd, $sigma, $upperlim_sigma, 
								      $fd_upperlim_lim, $thisobject_data_detection_flag_yn, 
								      $thisobjec_data_upperlim_flag_yn, $fd_upperlim_lim);
	      
	    }




	  } # keep looping until we stop finding lower chisq valies


	  # record the best-fitting parameters for posterity (i.e., for printing to the out file):

	  if ($best_chisq_sofar < $bestchisq){
	    if ($VERBOSE){print "  this downhill repeat found a new best model \n";}
	    $bestix = $current_model_ix;
	    $bestchisq = $best_chisq_sofar;
	    $flux_scale= $best_fluxscale_sofar;
	  }
	  if ($VERBOSE){print "\n"}

	} # end downhill UL method
      } # end the loop for downhill repeats 
    } # end of UL fitting if statement
    
    
#print "$bestix \n"; 
#print "  this iteration's bestchisq and fluxscale:   $bestchisq $flux_scale \n";
#print "\n";    




      
    # CHI-SQUARE FITTING DONE!

############################################################

    #store results if this is the 0th montecarlo iteration - i.e. if this is the actual, unperturbed data:
    if ($MC_fit_iteration == 0){
      $chisq_actual = $chisq;
#      $ix_actual = $ix;
      $bestix_actual = $bestix;
      $bestchisq_actual = $bestchisq;
      $ss_actual = $ss; 
      $flux_scale_actual = $flux_scale;
    }
    # and instead, if this is a perturbed MC iteration, store the best-fitting chisq in the *actual* chisq matrix...
    if ($MC_fit_iteration >= 0){
      set $bestchisq_MC, $MC_fit_iteration, $chisq_actual->index($bestix); 
      # ...and also store the results of MC results:
      $chisq_mc[$MC_fit_iteration] = $chisq;
#      $ix_mc[$MC_fit_iteration] = $ix;
      $bestix_mc[$MC_fit_iteration] = $bestix;
      $bestchisq_mc[$MC_fit_iteration] = $bestchisq;
      $ss_mc[$MC_fit_iteration] = $ss; 
      $flux_scale_mc[$MC_fit_iteration] = $flux_scale;
    }


  }


  if ($MC_PHOT_NITER > 0){
    if (lc($ERRORBARS_TYPE) =~ 'frac'){
      $bestchisq_MC_ordered = qsort ($bestchisq_MC);
      $chisq_level = interpol ($ERRORBARS_LEVEL,
			      (sequence(nelem $bestchisq_MC_ordered))/(nelem $bestchisq_MC_ordered), 
			      $bestchisq_MC_ordered);
      $Dchisq_level = $chisq_level-$bestchisq_actual;
    }
  }





  # output the best-fit results for this object:

  #    first output the object's parameters:
  for ($parameter = 0; $parameter < $Ndata_params; $parameter++){
    printf OUTPUT_FILE @DATA_PARAM_OUTPUTFORMAT[$parameter], ($data_param->index ($object))->index($parameter); print OUTPUT_FILE "  ";
  }
  #    now output the model parameters for the best-fitting model:
  for ($parameter = 0; $parameter < $Nmodel_params; $parameter++){

    if (@MODEL_PARAM_FOURTHCOL[$parameter] =~ "fluxscale"){
      $fluxScaleMultiplier=$flux_scale_actual;
      $fluxScaleMultiplierArray = $ss_actual;
    }
    else {
      $fluxScaleMultiplier = 1;
      $fluxScaleMultiplierArray = $ss_actual*0+1;
    }

    # print the best-fit value:
    printf OUTPUT_FILE @MODEL_PARAM_OUTPUTFORMAT[$parameter], (($model_param->index ($parameter))->index($bestix_actual))*$fluxScaleMultiplier; print OUTPUT_FILE " ";
    # and the allowable range:
    if (lc ($ERRORBARS_TYPE) !~ 'none'){
      $max = max (where (($model_param->index($parameter))*$fluxScaleMultiplierArray, $chisq_actual <= $bestchisq_actual + $Dchisq_level));
      $min = min (where (($model_param->index($parameter))*$fluxScaleMultiplierArray, $chisq_actual <= $bestchisq_actual + $Dchisq_level));
      printf OUTPUT_FILE @MODEL_PARAM_OUTPUTFORMAT[$parameter], $min; print OUTPUT_FILE " ";
      printf OUTPUT_FILE @MODEL_PARAM_OUTPUTFORMAT[$parameter], $max; print OUTPUT_FILE " ";
    }
    print OUTPUT_FILE "  ";

  }




  #   and now output the flux scale normalization :
  printf OUTPUT_FILE "%7.4e ", $flux_scale_actual;

  if (lc($ERRORBARS_TYPE) !~ 'none'){
    $max = max (where ($ss_actual, $chisq_actual <= $bestchisq_actual + $Dchisq_level));
    $min = min (where ($ss_actual, $chisq_actual <= $bestchisq_actual + $Dchisq_level));
    printf OUTPUT_FILE "%7.4e ", $min; print OUTPUT_FILE " ";
    printf OUTPUT_FILE "%7.4e ", $max; print OUTPUT_FILE " ";
  }

	
  # and now output the chisq information: 		
  printf OUTPUT_FILE "  %7.4e ",  $bestchisq_actual;
  if (lc($ERRORBARS_TYPE) !~ 'none'){ printf OUTPUT_FILE "%7.4e ", $Dchisq_level;}

  # end the putput line:
  print OUTPUT_FILE "\n";




  # output the chisq matrices and also make the parameter files for making the bestfit spectra:
  for ($parameter = 0; $parameter < $Nmodel_params; $parameter++){
    $argv = @DATA_PARAM_FOURTHCOL[$parameter];
    if (uc $argv =~ "OBJID"){
      $OBJID = ($data_param->index ($object))->index($parameter);
    }
  }
  $bestfit_makebbsed_outfilename = $MODEL_BBSED_FILE . "." . $OBJID. ".bbsed_spec";
  $chisq_outfilename = $MODEL_BBSED_FILE . "." . $OBJID. ".chisq.fits";





  #outputting chisq matrix:
  if ($CHISQ_MATRIX_YN =~ 'yes'){
    $outchisqmatrix=$chisq_actual;
    for ($parameter = 0; $parameter < $Nmodel_params; $parameter++){
      #    if (@MODEL_PARAM_FOURTHCOL[$parameter] =~ "fluxscale"){$fluxScaleMultiplier=$flux_scale;}
      #    else {$fluxScaleMultiplier = 1;}
      $fluxScaleMultiplier = 1;
      $outchisqmatrix = append ($outchisqmatrix, ($fluxScaleMultiplier*($model_param->index($parameter))));
    }
    reshape ($outchisqmatrix, ((nelem $outchisqmatrix)/($Nmodel_params+1)), ($Nmodel_params+1));
    ($outchisqmatrix) -> wfits($CHISQ_MATRIX_DIR . $chisq_outfilename);
  }




  # outputting instructions file for making best-fit spectra:
  if ($BESTFIT_SPECTRA_YN =~ 'yes'){
    $append_params = "";
    for ($parameter = 0; $parameter < $Nmodel_params; $parameter++){
      $argv = @MODEL_PARAM_FOURTHCOL[$parameter];
      if (uc $argv =~ "REDSHIFT"){
	$append_params = $append_params . " -REDSHIFTS value " . (($model_param->index ($parameter))->index($bestix_actual));
      }
      if (uc $argv =~ "MODEL_NUMBER"){
	$append_params = $append_params . " -MODEL_NUMBERS value " . (($model_param->index ($parameter))->index($bestix_actual));;
      }
      if (uc $argv =~ "EXTINCTION_EBV"){
	$append_params = $append_params . " -EXTINCTION_EBV value " . (($model_param->index ($parameter))->index($bestix_actual));;
      }
    }
    print BESTFIT_SPECTRA_MKFILE "make_bbsed " . $paramfile . " "
      . $BBSEDFILE_PARAMS
	. "-OUTPUT_FILE ". $bestfit_makebbsed_outfilename
	  . $append_params
	    . " -OUTPUT_OVERWRITE no -SPECMODE spectra "
	      . "\n";


  }



  # output the best-fit results for EACH of the MC realizations: 


  for ($MC_fit_iteration =1 ; $MC_fit_iteration <= $MC_PHOT_NITER; $MC_fit_iteration++){
  #    first output the object's parameters:
    for ($parameter = 0; $parameter < $Ndata_params; $parameter++){
      printf OUTPUT_MC_FILE @DATA_PARAM_OUTPUTFORMAT[$parameter], ($data_param->index ($object))->index($parameter); print OUTPUT_MC_FILE "  ";
    }
    #    now output the model parameters for the best-fitting model:
    for ($parameter = 0; $parameter < $Nmodel_params; $parameter++){
      
      if (@MODEL_PARAM_FOURTHCOL[$parameter] =~ "fluxscale"){
	$fluxScaleMultiplier=$flux_scale_mc[$MC_fit_iteration];
	$fluxScaleMultiplierArray = $ss_mc[$MC_fit_iteration];
      }
      else {
	$fluxScaleMultiplier = 1;
	$fluxScaleMultiplierArray = $ss_mc[$MC_fit_iteration]*0+1;
      }
      
      # print the best-fit value:
      printf OUTPUT_MC_FILE @MODEL_PARAM_OUTPUTFORMAT[$parameter], (($model_param->index ($parameter))->index($bestix_mc[$MC_fit_iteration]))*$fluxScaleMultiplier; print OUTPUT_MC_FILE " ";
      # and the allowable range:
      if (lc ($ERRORBARS_TYPE) !~ 'none'){
	$max = max (where (($model_param->index($parameter))*$fluxScaleMultiplierArray, $chisq_mc[$MC_fit_iteration] <= $bestchisq_mc[$MC_fit_iteration] + $Dchisq_level));
	$min = min (where (($model_param->index($parameter))*$fluxScaleMultiplierArray, $chisq_mc[$MC_fit_iteration] <= $bestchisq_mc[$MC_fit_iteration] + $Dchisq_level));
	printf OUTPUT_MC_FILE @MODEL_PARAM_OUTPUTFORMAT[$parameter], $min; print OUTPUT_MC_FILE " ";
	printf OUTPUT_MC_FILE @MODEL_PARAM_OUTPUTFORMAT[$parameter], $max; print OUTPUT_MC_FILE " ";
      }
      print OUTPUT_MC_FILE "  ";
      
    }
    
    #   and now output the flux scale normalization :
    printf OUTPUT_MC_FILE "%7.4e ", $flux_scale_mc[$MC_fit_iteration];
    if (lc($ERRORBARS_TYPE) !~ 'none'){
      $max = max (where ($ss_mc[$MC_fit_iteration], $chisq_mc[$MC_fit_iteration] <= $bestchisq_mc[$MC_fit_iteration] + $Dchisq_level));
      $min = min (where ($ss_mc[$MC_fit_iteration], $chisq_mc[$MC_fit_iteration] <= $bestchisq_mc[$MC_fit_iteration] + $Dchisq_level));
      printf OUTPUT_MC_FILE "%7.4e ", $min; print OUTPUT_MC_FILE " ";
      printf OUTPUT_MC_FILE "%7.4e ", $max; print OUTPUT_MC_FILE " ";
    }
    
    # and now output the chisq information: 		
    printf OUTPUT_MC_FILE "  %7.4e ",  $bestchisq_mc[$MC_fit_iteration];
    if (lc($ERRORBARS_TYPE) !~ 'none'){ printf OUTPUT_MC_FILE "%7.4e ", $Dchisq_level;}
    
    # end the putput line:
    print OUTPUT_MC_FILE "\n";
    


  }
}

close (OUTPUT_FILE);
close (OUTPUT_MC_FILE);
close (BESTFIT_SPECTRA_MKFILE);




















##############################  FUNCTIONS / SUBROUTINES  ####################################



######################
# subroutine to evaluate chisq of all models in the UL case:

sub find_allmodelsChisq_UL {

  my ($fm, $Nmodels, $fd, $sigma, $upperlim_sigma, 
   $fd_upperlim_lim, $thisobject_data_detection_flag_yn, 
   $thisobjec_data_upperlim_flag_yn, $fd_upperlim_lim) = @_;
  
  my ($fm_thismodel, $this_chisqUL, $this_s, $niter, $smin, $smax, $s, $ss, $chisqUL, $fmx, $nelem, $slice, $chisq);

  $chisqUL = pdl(0);
  $s = pdl(0); 
  $fmx = $fm->xchg(0,1);
      
  for ($model=0; $model < $Nmodels; $model++){       # <- go through the models one by one 
    $fm_thismodel = $fmx->index($model);   # <- gets the fluxes for just the model we are working on
    
    # here is the golden-section minimum finder for the upper limit case: 
    #    define the intial minimum and maximum that golden section will start with; 
    #    the max corresponds to the case where we shift all the model points up so 
    #    that all but one are higher than the data, with the lowest one being equivalent to the data: 

    $smin = 0; 
    $smax = max ($fd/$fm_thismodel *$thisobject_data_detection_flag_yn  
		 + $fd_upperlim_lim/$fm_thismodel *$thisobject_data_upperlim_flag_yn);
	  
    ($this_chisqUL, $this_s, $niter) = golden_section ($smin, $smax, $thisobject_data_detection_flag_yn, 
						       $thisobjec_data_upperlim_flag_yn, $fd, $fm_thismodel, 
						       $sigma, $upperlim_sigma, $fd_upperlim_lim); 
    $chisqUL = append $chisqUL, $this_chisqUL;
    $s = append $s, $this_s;
  }	

  $nelem = nelem $s;
  $slice = '1:' . ($nelem-1);
  $chisq = $chisqUL->slice($slice); 
  $ss=$s->slice($slice);	

  return ($chisq, $ss) 
}



##################################################################################
# golden section subtroutine:

sub golden_section {

  my ($smin, $smax, $thisobject_data_detection_flag_yn, $thisobjec_data_upperlim_flag_yn, $fd, $fm_thismodel, $sigma, $upperlim_sigma, $fd_upperlim_lim) = @_;

  $R = 0.6180339; 
  $C = 1-$R;
  $tol = 1e-8;

  $ax = $smin; 
  $cx = $smax; 
  $bx = $R*$smax;

  $x0 = $ax; 
  $x3 = $cx; 

  if (abs($cx-$bx) > abs($bx-$ax)){
    $x1=$bx;
    $x2=$bx+$C*($cx-$bx);}
  else{
    $x2=$bx;
    $x1=$bx-$C*($bx-$ax);
  }

#  print "$x0, $x1, $x2, $x3 \n";

  my $f1 = chisq_UL ($x1, $thisobject_data_detection_flag_yn, $thisobjec_data_upperlim_flag_yn, $fd, $fm_thismodel, $sigma, $upperlim_sigma, $fd_upperlim_lim); 
  my $f2 = chisq_UL ($x2, $thisobject_data_detection_flag_yn, $thisobjec_data_upperlim_flag_yn, $fd, $fm_thismodel, $sigma, $upperlim_sigma, $fd_upperlim_lim);
  
  $iter = 0; 

  while (abs($x3-$x0) > $tol*(abs($x1)+abs($x2))){
  
    $iter ++; 
#    print "$iter, $x0, $x1, $x2    $f0, $f1, $f2\n"; 
    if ($f2 < $f1){
      $x0=$x1; 
      $x1=$x2;
      $x2=$R*$x1+$C*$x3;
      $fo=$f1;
      $f1=$f2;
      $f2=chisq_UL ($x2, $thisobject_data_detection_flag_yn, $thisobjec_data_upperlim_flag_yn, $fd, $fm_thismodel, $sigma, $upperlim_sigma, $fd_upperlim_lim);
    }
    else{
      $x3=$x2;
      $x2=$x1;
      $x1=$R*$x2+$C*$x0;
      $f3=$f2;
      $f2=$f1;
      $f1=chisq_UL ($x1, $thisobject_data_detection_flag_yn, $thisobjec_data_upperlim_flag_yn, $fd, $fm_thismodel, $sigma, $upperlim_sigma, $fd_upperlim_lim)
    }
  }
  if ($f1<$f2){
    $golden=$f1; 
    $xmin = $x1;
  }
  else{
    $golden=$f2; 
    $xmin=$x2; 
  }
#  print "golden, xmin = $golden  $xmin \n";
  return ($golden, $xmin, $iter); 
}


#######################################################
# the upper limit chisq function is evaluated here: 
sub chisq_UL {
  my ($s, $thisobject_data_detection_flag_yn, $thisobjec_data_upperlim_flag_yn, $fd, $fm_thismodel, $sigma, $upperlim_sigma, $fd_upperlim_lim ) = @_;

  # the sum for detected bands:
  my $chisqUL_detections = sumover ($thisobject_data_detection_flag_yn * ((($fd - $s * $fm_thismodel)/$sigma)**2)); 

  # the sum for the undetected bands:
  my $chisqUL_nondetections = -2*log((3.1415926/2)**0.5 
				     * (sumover ( $thisobject_data_upperlim_flag_yn * $upperlim_sigma
						  * (erf(($fd_upperlim_lim - $s*$fm_thismodel)/($upperlim_sigma*2**0.5)) 
						     + erf(($s * $fm_thismodel) / ($upperlim_sigma * 2**0.5))))));  

  # now add the two: 
  my $chisqUL = $chisqUL_detections + $chisqUL_nondetections; 

  return $chisqUL; 
}










############################
# Subroutine that, given a grid of models and a focus model, selects out a subgrid of models that are neighbours (in parameter space) of the focus model. 
# The subroutine returns the model parameter subgrid, model fluxes subgrid, and a list of model index numbers, all for the neigbour models. 
# The "radius" of the neigbourhood (in terms of number of adjacent models) is specified by the user via a variable. 

sub select_neighbour_models {

  # input variables: 
  #       $model_param                the full grid of model parameter values
  #       $fm                         the full grid of model fluxes 
  #       $allmodels_ix               sequential list of model IDs (corresponds to $model_params, $fm, on a 1-to-1 basis)
  #       $current_model_ix           the index number (in the $allmodels_ix list) of the model for which neighbours are to be found
  #       $Nneighbours_submodelgrid   the number of neighbours (in each direction) that are allowed in each dimension of parameter space
  #
  # returned variables:
  #       $model_param_subgrid        the sugbrid of model parameters that contains the neighbours of the model for which neighbours are to be found
  #       $fm_subgrid                 the sugbrid of model fluxes
  #       $allmodels_ix_subgrid       the list of model index numbers of the neighbour models
  

  my ($model_param, $fm, $allmodels_ix, $Nmodel_params, $current_model_ix, $Nneighbours_submodelgrid) = @_; 

  my ($model_param_subgrid, $fm_subgrid, $allmodels_ix_subgrid, $parameter, $this_parameter, 
      $uniq_this_parameter, $Nuniq_this_parameter, $current_model_this_parameter, $uniq_this_parameter_ix, 
      $min_ix, $max_ix, $this_parameter_min_allowed, $this_parameter_max_allowed, 
      $model_param_this_param, $fm_this_param, $model_param_this_param_array);
      

  #make a temporary copy of the model parameters and fluxes that we can then trim down for the neighborhood:
  $model_param_subgrid = $model_param;
  $fm_subgrid = $fm;
  $allmodels_ix_subgrid = $allmodels_ix; 
      

  for ($parameter = 0; $parameter < $Nmodel_params; $parameter++){
    
    $this_parameter = ($model_param_subgrid -> index ($parameter)); # the values of this parameter (in the model subgrid)
    $uniq_this_parameter = qsort (uniq $this_parameter); 
    $Nuniq_this_parameter = nelem $uniq_this_parameter;  
    
    $current_model_this_parameter = (( $model_param->index($parameter))->index($current_model_ix)); # the value of this parameter in the current model
    
    #identify the element that corresponds to the current model value in the uniq list of elements for this subgrid
    $uniq_this_parameter_ix = ((qsorti (abs($uniq_this_parameter - $current_model_this_parameter)))->at(0)); 
    
    #identify the largest and smallest allowable value in the new subgrid we are generating:
    $min_ix = $uniq_this_parameter_ix - $Nneighbours_submodelgrid; if ($min_ix < 0){$min_ix=0;}
    $max_ix = $uniq_this_parameter_ix + $Nneighbours_submodelgrid; if ($max_ix >= $Nuniq_this_parameter) {$max_ix= $Nuniq_this_parameter-1;}
    
    # the allowed min and max values of this parameter: 
    $this_parameter_min_allowed = $uniq_this_parameter->index($min_ix);  
    $this_parameter_max_allowed = $uniq_this_parameter->index($max_ix);   
    
    $model_param_this_param = ($model_param_subgrid -> index ($parameter));
    $fm_this_param          = ((((ones $fm_subgrid)->xchg(0,1))*($model_param_this_param))->xchg(0,1));
    $model_param_this_param_array = ((($model_param_this_param)*((ones $model_param_subgrid)->xchg(0,1)))->xchg(0,1));
    
    ($model_param_subgrid)  = where ($model_param_subgrid, $model_param_this_param_array >= $this_parameter_min_allowed &
				     $model_param_this_param_array <= $this_parameter_max_allowed );
    ($fm_subgrid) =           where ($fm_subgrid, $fm_this_param >= $this_parameter_min_allowed &
				     $fm_this_param <= $this_parameter_max_allowed );
    ($allmodels_ix_subgrid) = where ($allmodels_ix_subgrid, $model_param_this_param >= $this_parameter_min_allowed & 
				     $model_param_this_param <= $this_parameter_max_allowed); 
    
    $model_param_subgrid = $model_param_subgrid->reshape ($Nmodel_params, ((nelem $model_param_subgrid)/$Nmodel_params));
    $fm_subgrid = $fm_subgrid->reshape ($Nmodel_mags, ((nelem $fm_subgrid)/$Nmodel_mags));
    
  } 
  
  return ($model_param_subgrid, $fm_subgrid, $allmodels_ix_subgrid); 

}#end of selection of the subgrid of models
	  

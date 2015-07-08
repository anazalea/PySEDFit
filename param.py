import os
from copy import deepcopy

#FIXME: maybe add model flux capability?
#FIXME: make error bar stuff backwards compatible
#FIXME: add in other optional parameters
fitsedOptParams = {"fitting_method" : "brute",
                   "model_flux_unit" : "mag",
                   "data_flux_unit" : "mag",
                   #restrict_model_param
                   "model_param" : [0, "model_param", "%.4f", "n"],
                   "output_overwrite" : "n",
                   "errorbar_method" : "None",
                   "errorbar_range" : [0.,],
                   #bestfit_spectra_file : ""
                   #data_flux_offsets
                   "montecarlo_iters" : 111,
                   "data_param" : [0, "data_param", "%.4f"],
                   "mag_softening" : 0.03}

fitsedMandParams = {"model_file" : "",
                    "model_flux_columns" : [0,],
                    "output_file" : "",
                    "data_file" : "",
                    "data_flux_columns" : [0,],
                    "data_error_columns" : [0,]}

def SetParams(pfile, mode="fitsed", args=None):
    """Read in parameters from a file and overwrite default params. 
    Parameter file must have lines of the form 'keyword value'.
    Anything after Hashtags are ignored. Whitespace is ignored.
    Parameters can also be input directly from command line with command line
    over-riding text file over-riding default values.

    Inputs:
        pfile: String of file path to parameter text file
        mode: Only sedfit right now. Controls which dictionaries to load

    Outputs:
        Dictionary containing parameter values.
    """
    #FIXME: add individual param value checks (i.e. spaces in file names)
    if mode == "fitsed":
        optDict = deepcopy(fitsedOptParams)
        mandDict = deepcopy(fitsedMandParams)
    else:
        string = "AHHH I don't understand the mode!"
        raise ValueError(string)

    keys = []
    values = []
    params = optDict
    params["model_param_column"] = []
    params["model_param_name"] = []
    params["model_param_output_format"] = []
    params["model_param_fluxscale"] = []
    with open(pfile, "r") as f:
        for line in f:
            line = line.split('#')[0].strip() #remove comments & whitespace
            if len(line) > 0: #ignore blank lines
                stuff = line.split(None, 1) #split into keyword, value
                if len(stuff) != 2:
                    string = ("Missing value or keyword in param file:\n"
                              "%s" % line)
                    raise ValueError(string)
                key, value = stuff
                key = key.lower() 
                key, value = CompatabilityCheck(key, value, params, mode)
                if key not in mandDict and key not in optDict: 
                    string = "In %s, %s is not a valid keyword" % (pfile,key)
                    raise KeyError(string)
                keys.append(key)
                values.append(value)
                

    if args is not None:
        clKeys = []
        clValues = []
        keyDex = []
        i = 0
        for x in args:
            if x[0] == '-' and x[1] not in '1234567890':
                clKeys.append(x)
                keyDex.append(i)
            i += 1
        temp = ""
        for i in xrange(len(args)):
            if i == 0:
                pass
            elif i == len(args) - 1:
                temp += args[i]
                clValues.append(temp.strip())
                temp = ''
            elif i in keyDex:
                clValues.append(temp.strip())
                temp = ''
            else:
                temp += args[i] + " "
        print clKeys
        print clValues
        if len(clKeys) != len(clValues):
            string = "Incorrect number of terms read in from command line"
            raise ValueError(string)
        for i in xrange(len(clKeys)):
            value = clValues[i]
            key = clKeys[i].lower() 
            if key[0] == '-':
                key = key[1:]
            key, value = CompatabilityCheck(key, value, params, mode)
            if key not in mandDict and key not in optDict: 
                string = "In command line, %s is not a valid keyword" % key
                raise KeyError(string)
            #overwrite file
            for j in xrange(len(keys)):
                if key == keys[j]:
                    values[j] = value
                else:
                    keys.append(key)
                    values.append(value)

    #check that all mandatory keys are present
    for k in mandDict.keys():
        if k not in keys:
            string = "Mandatory keyword %s is not present." % k
            raise KeyError(string)

    #change the value from string to appropriate type
    for i in xrange(len(keys)):
        if keys[i] in optDict:
            dic = optDict
        else:
            dic = mandDict
        if (not hasattr(dic[keys[i]], "strip") and 
            hasattr(dic[keys[i]], "__getitem__") or
            hasattr(dic[keys[i]], "__iter__")):
            values[i] = values[i].split()
            for j in xrange(len(values[i])):
                if j < len(dic[keys[i]]):
                    defaultV = dic[keys[i]][j]
                else:
                    defaultV = dic[keys[i]][-1]
                if type(defaultV).__name__ == 'int':
                    try:
                        values[i][j] = int(values[i][j])
                    except ValueError:
                        string = ("Something went wrong with key %s with "
                                  "value %s. Unable to turn value to int." % 
                                  (keys[i], values[i]))
                        raise ValueError(string)
                elif type(defaultV).__name__ == 'float':
                    try:
                        values[i][j] = float(values[i][j])
                    except ValueError:
                        string = ("Something went wrong with key %s with "
                                  "value %s. Unable to turn value to float." % 
                                  (keys[i], values[i]))
                        raise ValueError(string)
        else:
            if type(dic[keys[i]]).__name__ == 'int':
                try:
                    values[i] = int(values[i])
                except ValueError:
                    string = ("Something went wrong with key %s with value %s"
                              ". Unable to turn value to int." % 
                              (keys[i], values[i]))
                    raise ValueError(string)
            elif type(dic[keys[i]]).__name__ == 'float':
                try:
                    values[i] = float(values[i])
                except ValueError:
                    string = ("Something went wrong with key %s with value %s"
                              ". Unable to turn value to float." % 
                              (keys[i], values[i]))
                    raise ValueError(string)
        params[keys[i]] = values[i]
         
        #FIXME: is this the best way?
    if "model_param" in params:
        del params["model_param"]

    sanityCheck(params, mode)

    return params

def CompatabilityCheck(key, value, params, mode):
    """Do anything special to individual keywords or values.
    i.e. allow for backwards compatability, enforce formatting etc.

    Inputs:
        params: dictionary of parameters and values
        mode: Right now only sedfit. Controls stuff

    Outputs:
        Dictionary of parameters
    """
    if mode == "fitsed":
        optDict = deepcopy(fitsedOptParams)
        mandDict = deepcopy(fitsedMandParams)

    #backwards compatability
    if key == "model_bbsed_file":
        key = "model_file"
    if key == "model_mags":
        key = "model_flux_columns"
    if key == "data_mags":
        key = "data_flux_columns"
    if key == "data_uncertainties":
        key = "data_error_columns"
 
    if key == "model_param":
        if len(value.split()) > 4:
            string = ("Too many values in model_param. Expected at most 4. "
                      "Recieved %i: %s" % (len(value.split()), value))
            raise ValueError(string)
    if key == "data_param":
        if len(value.split()) > 3:
            string = ("Too many values in model_param. Expected at most 3. "
                      "Recieved %i: %s" % (len(value.split()), value))
            raise ValueError(string)

    

    #multiple param key functionality
    if key == "model_param":
        values = value.split()
        temp = deepcopy(optDict["model_param"])
        for i in xrange(len(values)):
            temp[i] = values[i]
        if temp[-1] == "fluxscale":
            temp[-1] = True
        else:
            temp[-1] = False
        try:
            params["model_param_column"].append(int(temp[0]))
        except ValueError:
            string = ("Something is wrong with model_param %s. Unable to turn "
                      "column into integer." % value)
            raise ValueError(string)
        params["model_param_name"].append(temp[1])
        params["model_param_output_format"].append(temp[2])
        params["model_param_fluxscale"].append(temp[3])
    if key == "data_param":
        values = value.split()
        temp = params["data_param"]
        for i in xrange(len(values)):
            temp[i] = values[i]
        try:
            params["data_param_column"].append(int(temp[0]))
        except ValueError:
            string = ("Something is wrong with data_param %s. Unable to turn "
                      "column into integer." % value)
            raise ValueError(string)
        params["data_param_name"].append(temp[1])
        params["data_param_output_format"].append(temp[2])

    return key, value

def sanityCheck(params, mode):
    #individual value checks
    for key, value in params.iteritems():

        if key in ["model_file", "data_file", "output_file"]:
            if " " in value:
                string = "Please remove space in %s" % key
                raise ValueError(string)

        if (key == "errorbar_method" and 
            value.lower() not in ["montecarlo", "dchisq","none"]):
            string = ("Value for key errorbar_method must be one of "
                      "'montecarlo' or 'dchisq'. Current value is %s." % value)
            raise ValueError(string)

        if key in ["model_flux_unit", "data_flux_unit"]:
            if value.lower() not in ["mag", "jansky"]:
                string = "%s must be one of mag or jansky." % key
                raise ValueError(string)
        
    if len(params["model_flux_columns"]) != len(params["data_flux_columns"]):
        string = ("model_flux_columns and data_flux_columns must have the same "
                  "number of columns.")
        raise ValueError(string)

    if len(params["data_error_columns"]) != len(params["data_flux_columns"]):
        string = ("data_error_columns and data_flux_columns must have the same "
                  "number of columns.")
        raise ValueError(string)
                        
                
                
                
                

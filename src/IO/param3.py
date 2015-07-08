import os
from copy import deepcopy
import re

#FIXME: add in other optional parameters

def SetParams(pfile, args=None):
    """Read in parameters from a file and overwrite default params. 
    Parameter file must have lines of the form 'keyword value'.
    Anything after Hashtags are ignored. Whitespace is ignored.
    Parameters can also be input directly from command line with command line
    over-riding text file over-riding default values.

    Inputs:
        pfile: String of file path to parameter text file
        paramList: List of Param objects. Says what keys to us and allows 
            for setting default values and error checking.

    Outputs:
        Dictionary containing parameter values.
    """
    paramList = deepcopy(fitsedParams)
    #paramList = deepcopy(paramList)
    allowedKeys = []
    usedKeys = []
    for x in paramList:
        allowedKeys.append(x.key)
 
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
                key, value = CompatabilityCheck(key, value)
                if key not in allowedKeys: 
                    string = "In %s, %s is not a valid keyword" % (pfile,key)
                    raise KeyError(string)
                
                for x in paramList:
                    if x.key == key:
                        x.setValue(value)
                        break
                if key in usedKeys and not x.multipleKey:
                    print("WARNING: Keyword %s used twice. Using value %s" %
                          (key, value))
                else:
                    usedKeys.append(key)

    if args is not None and len(args) > 0:
        clKeys = []
        clValues = []
        for x in args:
            if x[0] == '-' and x[1] not in '1234567890':
                clKeys.append(x)
        clValues = re.split("-\[a-z]+", ' '.join(args))
        clValues = map(str.strip, clValues)
        print(clValues)
        
        if len(clKeys) != len(clValues):
            string = "Incorrect number of terms read in from command line"
            raise ValueError(string)
        for i in xrange(len(clKeys)):
            value = clValues[i]
            key = clKeys[i].lower() 
            if key[0] == '-':
                key = key[1:]
            key, value = CompatabilityCheck(key, value)
            if key not in allowedKeys: 
                string = "In command line, %s is not a valid keyword" % key
                raise KeyError(string)
            #overwrite param file choices
            for x in paramList:
                if x.key == key:
                    x.setValue(value)
                    break
            if key in usedKeys and not x.multipleKey:
                print("WARNING: Keyword %s used twice. Using command line "
                      "value of %s" % (key, value))
            else:
                usedKeys.append(key)
  
    for x in paramList:
        #check that all mandatory keys are present
        if x.mandatory and x.key not in usedKeys:
            string = "Mandatory keyword %s is not present." % x.key
            raise KeyError(string)
        #check that associated keys are present
        if len(x.assocKeys) > 0:
            for key in x.assocKeys:
                if key not in usedKeys:
                    string = ("ERROR: When using %s, %s key and value must "
                              "also be declared in the param file or command "
                              "line." % (x.key, key))
                    raise KeyError(string)
        #set values of optional keywords to default if not present
        if not hasattr(x, "value"):
            x.value = x.defaultValue

    paramDict = {}
    for x in paramList:
        paramDict[x.key] = x.value

    #sanityCheck(paramDict)

    return paramDict

def CompatabilityCheck(key, value):
    """Allow for backwards compatability of key words.
    """
    if key == "model_bbsed_file":
        key = "model_file"
    if key == "model_mags":
        key = "model_flux_columns"
    if key == "data_mags":
        key = "data_flux_columns"
    if key == "data_uncertainties":
        key = "data_error_columns"

    return key, value

def sanityCheck(params, mode):
    #FIXME: move this to main module
        
    if len(params["model_flux_columns"]) != len(params["data_flux_columns"]):
        string = ("model_flux_columns and data_flux_columns must have the same "
                  "number of columns.")
        raise ValueError(string)

    if len(params["data_error_columns"]) != len(params["data_flux_columns"]):
        string = ("data_error_columns and data_flux_columns must have the same "
                  "number of columns.")
        raise ValueError(string)

#############################################################################

class Param(object):
    """Contains all information about a parameter that will be used to 
    control the flow of the main program. Paramaters are initiated with just a
    keyword (and other options for format and error handling), but no value.
    The value must be explicitly set using the setValue() method.

    Inputs:
        key -- The key word name
        defaultValue -- Optional default value
        dataType -- Data type the value should be. Can be list of types if 
            value will be a list. Default is str
        isList -- Boolean flag stating if the value can be a list
        multipleKey -- Boolean flag allowing for a key to be input multiple 
            times. Each value will be saved in a list, rather than 
            overwriting the previous value.
        mandatory -- Boolean flag stating whether it is necessary for the user
            to input a value for this keyword. i.e. the program cannot run
            without it and there is no default value.
        maxSize -- If the value will be a list, this is the max size allowed.
        allowedValues -- Optional list of possible values for this key. Any 
            other values will raise an exception.
        canBeNeg -- Boolean flag stating whether a value can be less than zero
        canHaveSpace -- Boolean flag stating whether a value can contain a 
            whitespace character.
        assocKeys -- List of associated keywords that must be present if this 
            key is used.
    """
    def __init__(self, key, defaultValue=None, dataType=str,  isList=False, 
                 multipleKey=False, mandatory=False, maxSize=None, 
                 allowedValues=None, canBeNeg=True, canHaveSpace=True,
                 assocKeys=[]):
        self.key = key
        self.defaultValue = defaultValue
        self.dataType = dataType
        self.isList = isList
        self.multipleKey = multipleKey
        self.mandatory = mandatory
        self.maxSize = maxSize
        self.allowedValues = allowedValues
        self.canBeNeg = canBeNeg
        self.canHaveSpace = canHaveSpace
        self.assocKeys = assocKeys

        if self.isList:
            if not isinstance(self.dataType, list):
                self.dataType = [self.dataType,]
            if not isinstance(self.defaultValue, list):
                self.defaultValue = [self.defaultValue,]
            if not isinstance(self.canBeNeg, list):
                self.canBeNeg = [self.canBeNeg,]

    def setValue(self, value):
        """Set the value associated with a keyword. Raise an exception if 
        the value is incompatible with the options chosen during initiation.

        Inputs:
            value -- The value to be set. Must be a string. If the value should
                be a list, then the elements should be separated by whitespace.
        """
        if not self.canHaveSpace and len(value.split()) > 1:
            string = "Please remove space in value of %s" % self.key
            raise ValueError(string)
        if self.isList:
            values = value.split()
            if self.maxSize is not None:
                if len(values) > self.maxSize:
                    string = ("Too many values in %s. Expected at most %i. "
                              "Recieved %i: %s" % 
                              (self.key, self.maxSize, len(values), value))
                    raise ValueError(string)
            newValue = []
            for i in xrange(len(values)):
                if i < len(self.defaultValue):
                    j = i 
                else:
                    j = -1
                if values[i].lower() in ["y", "yes", "t", "true"]:
                    values[i] = True
                elif values[i].lower() in ["n", "no", "f", "false"]:
                    values[i] = False
                try:
                    values[i] = self.dataType[j](values[i])
                except ValueError:
                    string = ("ERROR: Something went wrong with key %s with "
                              "value %s. Unable to convert %s to format %s." %
                              (self.key, value, values[i], self.dataType[j]))
                    raise ValueError(string)
                if i < len(self.canBeNeg):
                    k = i 
                else:
                    k = -1
                if not self.canBeNeg[k] and values[i] < 0:
                    string = "Value of %s cannot be negative." % self.key
                    raise ValueError(string)
                newValue.append(values[i])
            if len(newValue) < len(self.defaultValue):
                for x in self.defaultValue[len(newValue):]:
                    newValue.append(x)
            value = newValue
        else:
            if self.allowedValues is not None:
                if value.lower() not in self.allowedValues:
                    string = ("ERROR: value for %s must be one of %s. Instead "
                              "value is %s." % 
                              (self.key, self.allowedValues, value))
                    raise ValueError(string)
                else:
                    value = value.lower()
            if value.lower() in ["y", "yes", "t", "true"]:
                    value = True
            elif value.lower() in ["n", "no", "f", "false"]:
                value = False
            try:
                value = self.dataType(value)
            except ValueError:
                string = ("ERROR: Something went wrong with key %s with "
                          "value %s. Unable to convert to format %s." %
                          (self.key, value, self.dataType))
                raise ValueError(string)

        if not self.canBeNeg and value < 0:
            string = "Value of %s cannot be negative." % self.key
            raise ValueError(string)

        if self.multipleKey:
            if not hasattr(self, "value"):
                self.value = [value,]
            else:
                self.value.append(value)
        else:
            self.value = value

def formatConversion(x):
    print x
    try:
        test = 2
        x % test
    except:
        string = "Unsupported format conversion string %s." % x
        raise ValueError(string)
    return x

#########################################################################
#FIXME: move to main module???
fitsedParams = [Param("fitting_method", defaultValue="brute", 
                      allowedValues=["brute", "tree"]),
                Param("model_flux_unit", defaultValue="mag",
                      allowedValues=["mag", "jansky"]),
                Param("data_flux_unit", defaultValue="mag",
                      allowedValues=["mag", "jansky"]),
                Param("model_param", isList=True, multipleKey=True,
                      defaultValue=[0, "model_param", "%.4f", False],
                      dataType=[int, str, formatConversion, bool],
                      canBeNeg=False),
                Param("output_overwrite", defaultValue=False, dataType=bool,
                      allowedValues=["y", "yes", "n", "no", 
                                     "t", "true", "f", "false"]),
                Param("errorbar_method", defaultValue="montecarlo",
                      allowedValues=["montecarlo", "dchisq", "none"]),
                Param("errorbar_range", isList=True, dataType=float, 
                      canBeNeg=False),
                Param("montecarlo_iters", defaultValue=300, dataType=int,
                      canBeNeg=False),
                Param("data_param", isList=True, multipleKey=True,
                      defaultValue=[0, "data_param", "%.4f"],
                      dataType=[int, str, formatConversion],
                      canBeNeg=False),
                Param("mag_softening", defaultValue=0.03, dataType=float,
                      canBeNeg=False),
                Param("model_file", mandatory=True, canHaveSpace=False),
                Param("output_file", mandatory=True, canHaveSpace=False),
                Param("data_file", mandatory=True, canHaveSpace=False),
                Param("model_flux_columns", mandatory=True, isList=True,
                      dataType=int, canBeNeg=False),
                Param("data_flux_columns", mandatory=True, isList=True,
                      dataType=int, canBeNeg=False),
                Param("data_error_columns", mandatory=True, isList=True,
                      dataType=int, canBeNeg=False),
                Param("brute_space",allowedValues=['color','flux'])
                      ]
                


               
        
                        
                
                
                
                

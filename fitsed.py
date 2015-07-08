'''
Created by Bobby

'''
from __future__ import print_function
from __future__ import division
import sys
sys.path.append('.')
import param

def main(params):
    
    print "Now running fitsed with model file %s" % params['model_file']
    os.system('say "Welcome to pie sedfit"')
    
    # Read in models and data as specified by params, do stuff
    # Data
    dataFlux = np.genfromtxt(params.get(''))
    
 
if __name__ == "__main__":
    pfile = sys.argv[1]
    args = sys.argv[2:]
    params = param.SetParams(pfile, 'fitsed', args)
    main(params)

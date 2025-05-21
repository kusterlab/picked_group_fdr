from __future__ import print_function

import sys
import os
import json

configDir = sys.argv[1]
localDir = sys.argv[2]
with open(os.path.join(configDir, 'config.json')) as f:
    config = json.load(f)
    if "options" in config:
        optionsConfig = config["options"]
    elif "fastaDigestOptions" in config:
        optionsConfig = config["fastaDigestOptions"]
    else:
        raise ValueError("Could not find \"option\" or \"fastaDigestOptions\" in config.json")
    
    if sys.argv[3] == "PERC_RESULT_FILES":
        if config["uploads"]["prosit_target.psms"] and config["uploads"]["prosit_decoy.psms"]:
            print(os.path.join(localDir, "prosit_target.psms"), os.path.join(localDir, "prosit_decoy.psms"))
        else:
            print(os.path.join(localDir, "out/percolator/andromeda.mokapot.psms.txt"), os.path.join(localDir, "out/percolator/andromeda.mokapot.decoy.psms.txt"))
    
    elif sys.argv[3] == "PROSIT_FLAG":
        if config["uploads"]["prosit_target.psms"] and config["uploads"]["prosit_decoy.psms"]:
            print("--pout_input_type prosit")
        else:
            print()
    
    elif sys.argv[3] == "DIGEST_PARAMS":
        digest_params = []
        if optionsConfig["protease"]:
            digest_params.append("--enzyme " + optionsConfig["protease"])
        
        if optionsConfig.get("digestion", False):
            digest_params.append("--digestion " + optionsConfig["digestion"])
        
        if optionsConfig["minPeptideLength"]:
            digest_params.append("--min-length " + str(optionsConfig["minPeptideLength"]))
        
        if optionsConfig["maxPeptideLength"]:
            digest_params.append("--max-length " + str(optionsConfig["maxPeptideLength"]))
        
        if optionsConfig["specialAAs"]:
            digest_params.append("--special-aas " + optionsConfig["specialAAs"])
        else:
            digest_params.append("--special-aas \\\"\\\"")
            
        if not isinstance(optionsConfig["missedCleavages"], bool): # check for bool, as value can be 0 which also evaluates to False
            digest_params.append("--cleavages " + str(optionsConfig["missedCleavages"]))
            
        print(" ".join(digest_params))
    
    elif sys.argv[3] == "PICKED_GROUP_FDR_EXTRA_PARAMS":
        if optionsConfig.get("minLFQPeptides", False):
            print("--lfq_min_peptide_ratios " + str(optionsConfig["minLFQPeptides"]))
        else:
            print()
        
    elif sys.argv[3] == "NUM_THREADS":
        if "numThreads" in config:
            print(config["numThreads"])
        else:
            print(1)
     

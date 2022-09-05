from __future__ import print_function

import sys
import os
import json

configDir = sys.argv[1]
localDir = sys.argv[2]
with open(os.path.join(configDir, 'config.json')) as f:
    config = json.load(f)
    
    if sys.argv[3] == "PERC_RESULT_FILES":
        if config["uploads"]["prosit_target.psms"] and config["uploads"]["prosit_decoy.psms"]:
            print(os.path.join(localDir, "prosit_target.psms"), os.path.join(localDir, "prosit_decoy.psms"))
        else:
            print(os.path.join(localDir, "out/percolator/andromeda.mokapot.psms.txt"), os.path.join(localDir, "out/percolator/andromeda.mokapot.decoys.psms.txt"))
    
    elif sys.argv[3] == "PROSIT_FLAG":
        if config["uploads"]["prosit_target.psms"] and config["uploads"]["prosit_decoy.psms"]:
            print("--pout_input_type prosit")
        else:
            print()
    
    elif sys.argv[3] == "DIGEST_PARAMS":
        digest_params = []
        if config["options"]["protease"]:
            digest_params.append("--enzyme " + config["options"]["protease"])
        
        if config["options"]["minPeptideLength"]:
            digest_params.append("--min-length " + str(config["options"]["minPeptideLength"]))
        
        if config["options"]["maxPeptideLength"]:
            digest_params.append("--max-length " + str(config["options"]["maxPeptideLength"]))
        
        if config["options"]["specialAAs"]:
            digest_params.append("--special-aas " + config["options"]["specialAAs"])
        else:
            digest_params.append("--special-aas \\\"\\\"")
            
        if not isinstance(config["options"]["missedCleavages"], bool):
            digest_params.append("--cleavages " + str(config["options"]["missedCleavages"]))
            
        print(" ".join(digest_params))
     
    elif sys.argv[3] == "NUM_THREADS":
        if "numThreads" in config:
            print(config["numThreads"])
        else:
            print(1)
     

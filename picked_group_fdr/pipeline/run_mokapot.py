import sys
import os
import logging
import time
import datetime

import mokapot
from mokapot import __version__


def run_mokapot(perc_folder, test_fdr, train_fdr):
    psms_result_file = os.path.join(perc_folder, 'andromeda.mokapot.psms.txt')
    if os.path.isfile(psms_result_file):
      print(f"Found mokapot output file {psms_result_file}, remove this file to rerun mokapot.")
      sys.exit(0)

    start = time.time()

    # Setup logging
    verbosity_dict = {
        0: logging.ERROR,
        1: logging.WARNING,
        2: logging.INFO,
        3: logging.DEBUG,
    }

    logging.basicConfig(
        filename=os.path.join(perc_folder, "andromeda.log"),
        format=("[{levelname}] {message}"),
        style="{",
        level=verbosity_dict[2],
    )

    logging.info("mokapot version %s", str(__version__))
    logging.info("Written by William E. Fondrie (wfondrie@uw.edu) in the")
    logging.info(
        "Department of Genome Sciences at the University of " "Washington."
    )
    logging.info("Command issued:")
    logging.info("%s", " ".join(sys.argv))
    logging.info("")
    logging.info("Starting Analysis")
    logging.info("=================")

    psms = mokapot.read_pin(os.path.join(perc_folder, "andromeda.tab"))
    results, models = mokapot.brew(psms, test_fdr = test_fdr)
    results.to_txt(dest_dir = perc_folder, file_root="andromeda", decoys = True)

    total_time = round(time.time() - start)
    total_time = str(datetime.timedelta(seconds=total_time))

    logging.info("")
    logging.info("=== DONE! ===")
    logging.info("mokapot analysis completed in %s", total_time)


if __name__ == "__main__":
    test_fdr = float(sys.argv[1])
    train_fdr = float(sys.argv[2]) # currently not supported
    perc_folder = sys.argv[3]
    
    run_mokapot(perc_folder, test_fdr, train_fdr)

import sys
if not '-m' in sys.argv:
    from .pipeline import (
        run_picked_group_fdr_all,
        run_andromeda_to_pin,
        run_merge_pout,
        run_mokapot,
        run_picked_group_fdr,
        run_picked_group_fdr_percolator_input,
        run_update_evidence,
        run_filter_fdr_maxquant
    )

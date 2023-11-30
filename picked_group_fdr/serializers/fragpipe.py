from typing import Union

import numpy as np


def fragpipe_format_extra_columns(x: Union[str, int, float]) -> str:
    if type(x) == str:
        return x
    if type(x) == int:
        return x
    if np.isnan(x) or x == 0.0:
        return 0.0
    return "%.5g" % (x)

"""
Convert 1 2D+t RSML <-> n 2D RSML.
"""

import numpy as np
from openalea.mtg import MTG
from copy import deepcopy
from typing import List, Dict


#############################################################################################################################################################
# 1 2D+t -> n 2D
#############################################################################################################################################################

def extract_mtg_at_time_t(g: MTG, t: int) -> MTG:
    g_new = deepcopy(g)

    time_prop = g_new.property("time")
    time_h_prop = g_new.property("time_hours")
    diameter_prop = g_new.property("diameter")
    geometry_prop = g_new.property("geometry")

    to_remove = []
    for v, serie in time_prop.items():
        first_t = serie[0]
        if first_t > t:
            to_remove.append(v)
        else:
            idx = max(i for i, tau in enumerate(serie) if tau <= t)

            _truncate_lists(time_prop, idx, v)
            _truncate_lists(time_h_prop, idx, v)
            _truncate_lists(diameter_prop, idx, v)
            _truncate_lists(geometry_prop, idx, v)

            # if list is empty has 1 or less elements, remove vertex
            if len(geometry_prop[v]) <= 1:
                to_remove.append(v)

    for v in to_remove:
        try:
            g_new.remove_tree(v)
        except Exception:
            g_new.remove_vertex(v, reparent_child=False)

    return g_new


#############################################################################################################################################################
# Utils
#############################################################################################################################################################


def _truncate_lists(prop: Dict[int, List], idx: int, v: int) -> None:
    """ 
    Truncate lists in a property dictionary at a given index for a specific vertex.
    """
    val = prop.get(v)
    if isinstance(val, (list, tuple)) and len(val) > idx + 1:
        prop[v] = val[: idx + 1] 
        

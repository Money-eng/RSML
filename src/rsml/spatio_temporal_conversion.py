"""
Convert 1 2D+t RSML <-> n 2D RSML.
"""

import os
from openalea.mtg import MTG
from copy import deepcopy
from typing import List, Dict
from SimpleITK import Transform
import SimpleITK as sitk
from re import findall
from rsml import rsml2mtg



#############################################################################################################################################################
# n 2D -> 1 2D+t
#############################################################################################################################################################

def combine_mtgs_at_time(g_list: List[MTG], time_step: int = 1) -> MTG:
    """Making 1 2D+t MTG from n 2D MTGs.

    Args:
        g_list (List[MTG]): List of 2D MTGs to merge.
        time_step (int, optional): Time step for the merged MTG. Defaults to 1.

    Returns:
        MTG: Merged 2D+t MTG.
    """
    
    # Transform geometry of mtg if exists (simpleITK transform)
    # Match Roots 
    # Spatio-temporal regression
    # TODO

def transform_mtg_geometry(g: MTG, transform: Transform) -> MTG:
    """Apply a SimpleITK transform to the geometry of an MTG.

    Args:
        g (MTG): The input MTG.
        transform (Transform): The SimpleITK transform to apply.

    Returns:
        MTG: The transformed MTG.
    """
    g_transformed = deepcopy(g)

    for v in g_transformed.vertices():
        geometry = g_transformed.properties().get("geometry", {}).get(v, None)
        if geometry:
            for i, point in enumerate(geometry):
                point_2d = (point[0], point[1], 0)
                point_transformed = transform.TransformPoint(point_2d)
                geometry[i] = (point_transformed[0], point_transformed[1])  # no z-coordinate

    return g_transformed

def match_plants_list(g_list: List[MTG], max_distance: float = None) -> List[Dict]:
    """Match plants across a list of MTGs.

    Args:
        g_list (List[MTG]): List of MTGs to match.
        max_distance (float, optional): Maximum distance for matching. Defaults to None.

    Returns:
        List[Dict]: List of match results between consecutive MTGs.
    """
    from rsml.matching import match_plants

    match_results = []
    for i in range(len(g_list) - 1):
        match_result = match_plants(g_list[i], g_list[i + 1], max_distance=max_distance)
        match_results.append(match_result)

    return match_results

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
        
def load_transform_from_folder(folder: str) -> List[Transform]:
    """Load SimpleITK transforms from a folder.

    Args:
        folder (str): Path to the folder containing the transforms.
    Returns:

        List[Transform]: List of SimpleITK transforms.
    """
    transforms = []
    # sorted must be regarding number in filename
    print(sorted(os.listdir(folder), key=lambda x: int(findall(r'\d+', x)[0])))
    for file in sorted(os.listdir(folder), key=lambda x: int(findall(r'\d+', x)[0])):
        if file.endswith(".txt") or file.endswith(".tif") or file.endswith(".tfm"):
            transform_path = os.path.join(folder, file)
            transform = sitk.ReadTransform(transform_path)
            transforms.append(transform)
    return transforms

graph = rsml2mtg("/home/loai/Images/DataTest/UC1_data/Val/230629PN019/61_graph.rsml")
# define small affine transform (x = x +10, y = y + 0)
transforms = load_transform_from_folder("/home/loai/Images/DataTest/UC3/Output_Data/Process/B73_R05_01/Transforms_1")
# plot transforms[0] and see effect
print(transforms[0])
print(transforms[0].GetParameters())
print(transforms[0].GetFixedParameters())
print(transforms[0].GetName())
print(transforms[0].GetInverse())
# apply transform to geometry of mtg
transformed_graph = transform_mtg_geometry(graph, transforms[0])
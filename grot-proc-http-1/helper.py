import numpy as np

def denumpyfylist(colection):
    if len(colection) == 0:
        return colection    
    newList = []
    for item in colection:
        if type(item) is list:
            item = denumpyfylist(item)
        if type(item) is np.ndarray:
            newList.append(item.tolist())
        else:
            newList.append(item)
    return newList
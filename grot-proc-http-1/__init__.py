import logging
import base64
from .grot.solver import Build
from .grot.bmp import load_im, create_geom
from .grot.prep import materials, thicks

import azure.functions as func


def main(req: func.HttpRequest) -> func.HttpResponse:
    logging.info('Python HTTP trigger function processed a request.')

    try:
        req_body = req.get_json()
    except ValueError:
        return func.HttpResponse(f"Wrong request body {ValueError.__cause__}", status_code=400)
        pass
    else:
        problem = req_body.get('problem')
        material = req_body.get('material')
        dimensionUnit = req_body.get('dimensionUnit')
        scale = req_body.get('scale')
        thickness = req_body.get('thickness')
        load = req_body.get('load')
        solverType = req_body.get('solver')
        disp = req_body.get('disp')
        stress = req_body.get('stress')
        deformed = req_body.get('deformed')
        plast = req_body.get('plast')
        image_as_string = req_body.get('image')
        img_bytes = base64.b64decode(image_as_string)

        IMAGE = load_im(img_bytes)
        GEOM = create_geom(IMAGE)
        NODES = GEOM[0]
        ELES = GEOM[1].store()
        CONS = GEOM[2]
        BC_DICT = GEOM[3]
        PROB_DICT = GEOM[4]
        MAT = materials(ELES)
        MAT.add(material)
        MAT.assignall(1)
        MAT.set_unit(dimensionUnit)
        MAT.set_scale(float(scale))
        SCALE = float(scale)
        THICKS = thicks(ELES, MAT)
        THICKS.add(float(thickness))
        THICKS.assignall(1)
        CONS.load(BC_DICT['magenta'], load.get('X'), load.get('Y'))
        constraints = CONS.store()
        STATE = problem

        SOL = Build(NODES, ELES, constraints, STATE, load_inc=1.0, scale=SCALE)
        disp = SOL.direct()
        return func.HttpResponse(str(disp), status_code=200)



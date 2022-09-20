import logging
import base64
import copy
import sys
import gc
import json
import pathlib
import numpy as np
from .grot.solver import Build
from .grot.bmp import load_im, create_geom
from .grot.prep import materials, thicks
from .grot.plast import Prepare, search
from .grot.prob import write
from .helper import denumpyfylist

import azure.functions as func


def main(req: func.HttpRequest) -> func.HttpResponse:
    logging.info('Python HTTP trigger function processed a request.')

    try:
        req_body = req.get_json()
    except ValueError:        
        f = open(pathlib.Path(__file__).parent /  "sample.dat", "r")
        sample = f.read()
        f.close()
        return func.HttpResponse(f"Wrong request body {ValueError.__cause__} \n Example: \n {sample}", status_code=400)
        pass
    else:
        name = req_body.get('name')
        problem = req_body.get('problem')
        material = req_body.get('material')
        dimensionUnit = req_body.get('dimensionUnit')
        scale = req_body.get('scale')
        thickness = req_body.get('thickness')
        load = req_body.get('load')
        solverType = req_body.get('solver')
        res_d = req_body.get('disp')
        stress = req_body.get('stress')
        deformed = req_body.get('deformed')
        plast = req_body.get('plast')
        plastIterations = req_body.get('plastIterations')
        probe = req_body.get('probe')
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
        result = {}

        if plast == 'no':
            disp = SOL.direct()
            result['disp'] = disp
            result['strains'] = SOL.strains_calc(disp)
        if plast == 'yes':
            disp = SOL.direct_plast()
            result['disp'] = disp
            disp_el = copy.copy(disp)
            strains = SOL.strains_calc(disp, msg=0)
            result['strains'] = strains
            strains_el = copy.deepcopy(strains)
            iter_res = Prepare(disp, strains, ELES)
            step_factor = iter_res.first_step(MAT)
        if (plast =='yes') and (step_factor < 1):
            load_step = step_factor
            steps_num = int(plastIterations)
            load_inc = (1 - step_factor) / (steps_num)
            flags_list = []
            eles_list = []
            sys.stdout.write("\r" + "Nonlinear plasticity solver iteration [" + str(1) + \
                            " of " + str(steps_num) + "]")
            sys.stdout.flush()
            for i in range(steps_num):
                load_step += load_inc
                check_res = iter_res.out()
                # Runge Kutta 2nd order procedure
                SOL.plast_update([], load_inc / 2.0)
                disp = SOL.direct_plast()
                strains = SOL.strains_calc(disp, msg=0)
                halfstep_strains = iter_res.halfstep(strains)
                plast_res = search(ELES, halfstep_strains, flags_list)

                eles_list = plast_res[0]
                flags_list = plast_res[1]
                stress2plast_list = plast_res[2]  # for residuals check

                sys.stdout.write("\r" + "Nonlinear plasticity solver iteration [" + \
                                str(i + 1) + " of " + str(steps_num) + "]")
                sys.stdout.flush()

                STATE = problem
                SOL.plast_update(eles_list, load_inc)
                MAT.assignplast(eles_list)

                disp = SOL.direct_plast()
                strains = SOL.strains_calc(disp, msg=0)

                final_results = iter_res.store(MAT, disp, strains, flags_list)
                disp = final_results[0]
                strains = final_results[1]
                eff_pl_strains = final_results[2]
                eff_pl_strains_rate = final_results[3]
                result['disp'] = disp
                result['strains'] = strains
                result['eff_pl_strains'] = eff_pl_strains
                result['eff_pl_strains_rate'] = eff_pl_strains_rate

                s2plast_corrected = []  # to calculate actual ratio, not ratio in hafstep
                for val in stress2plast_list:
                    val -= val * ((load_inc / 2.0) / (load_step - (load_inc / 2.0)))
                    s2plast_corrected.append(val)
                if len(s2plast_corrected) == 0:
                    s2plast_corrected.append(0)
                min_val = str(round(min(s2plast_corrected), 3))
                max_val = str(round(max(s2plast_corrected), 3))
                mean_val = str(round(sum(s2plast_corrected) / len(s2plast_corrected), 3))
                new_eles = str(len(eles_list))
                all_eles = str(len(flags_list))
                result['mean_val'] = mean_val
                result['min_val'] = min_val
                result['max_val'] = max_val
                result['new_eles'] = new_eles
                result['all_eles'] = all_eles
                result['eff_pl_strains'] = str(round(max(eff_pl_strains), 3))
                result['eff_pl_strains_rate'] = str(round(max(eff_pl_strains_rate), 3))          
            print("\nPlasticity analysis details created stored in results")
            print("")

            check_res = iter_res.out()
            res_disp = iter_res.residual_disp(disp_el)
            res_strains = iter_res.residual_strains(strains_el)
            strains = iter_res.store_plstrain(strains)
            result['check_res'] = check_res
            result['res_disp'] = res_disp
            result['res_strains'] = res_strains
            result['strains'] = strains

            disp_el, strains_el, iter_res, plast = None, None, None, None
            halfstep_strains, plast_res, final_results = None, None, None
        gc.collect()
        if probe is not None:
            write(probe, PROB_DICT, strains, name, MAT)

        print("")
        print("Task finished")

        denumpyfied_result = {}
        for key in result:
            if type(result[key]) is np.ndarray:
                denumpyfied_result[key] = result[key].tolist()
            else:
                if type(result[key]) is list:
                    denumpyfied_result[key] = denumpyfylist(result[key])
                else:
                    denumpyfied_result[key] = result[key]

        return func.HttpResponse(json.dumps(denumpyfied_result), status_code=200)




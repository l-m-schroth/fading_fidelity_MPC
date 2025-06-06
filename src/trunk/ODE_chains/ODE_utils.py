
import numpy as np

"""
This file con be used to load the ODE parameters for the trunk that are consistent with the mujoco model.
"""

def params_dict_to_list(params_dict):
        return [params_dict["g"], params_dict["m"], params_dict["l"], params_dict["Iz"], params_dict["c"], params_dict["d"], params_dict["gear"]] 

def get_ode_params_dict():
    # params dict 
    g = 9.81
    r = 0.0005
    l_tot_2 = 0.25
    l_2 = l_tot_2/2 # l in the ode is the half length
    m_2 = 0.020158552860534508 # extracted masses and inertias directly from mujoco
    #m_2 = np.pi*l_tot_2*r**2*1000 # mass for 2 link chain default desity is 1000, mass is calculate by cylinder density times volume pi*r^2*total_length
    Iz_2_cog = 0.00011082000810461473
    Iz_2 = Iz_2_cog #+ m_2*l_2**2 apparently it is already shifted, as it behaves exactly like mojuco when removing the steiner part
    gear = 2
    params_2_links = {
        "g": g,
        "m": m_2, # default desity is 1000, mass is calculate by cylinder density times volume pi*r^2*total_length
        "l": l_2, 
        "Iz": Iz_2,#1/4*m_2*r**2 + 7/6*m_2*l_2**2,
        "c": 0.65,
        "d": 0.032,
        "gear": gear
    }
    m_4 = 0.010341075818066403 # extracted from mujoco
    Iz_4_cog = 1.5017794631863459e-05
    l_4 = l_2/2
    Iz_4 = Iz_4_cog #+ m_4*l_4**2
    params_4_links = {
        "g": g,
        "m": m_4, # default desity is 1000, mass is calculate by cylinder density times volume pi*r^2*total_length
        "l": l_4, 
        "Iz": Iz_4,#1/4*m_4*r**2 + 7/6*m_4*l_4**2,
        "c": 0.9,
        "d": 0.044,
        "gear": gear
    }
    m_8 = 0.005432337296832351 # extracted from mujoco
    l_8 = l_4/2
    Iz_8_cog = 2.2064984183718783e-06
    Iz_8 = Iz_8_cog #+ m_8*l_8**2
    params_8_links = {
        "g": g,
        "m": m_8, # default desity is 1000, mass is calculate by cylinder density times volume pi*r^2*total_length
        "l": l_8, 
        "Iz": Iz_8, #1/4*m_8*r**2 + 7/6*m_8*l_8**2,
        "c": 1.55,
        "d": 0.0765,
        "gear": gear
    }
    m_16 = 0.0029779680362153247 # extracted from mujoco
    l_16 = l_8/2
    Iz_16_cog = 3.7882422547229845e-07
    Iz_16 = Iz_16_cog #+ m_16*l_16**2
    params_16_links = {
        "g": g,
        "m": m_16, # default desity is 1000, mass is calculate by cylinder density times volume pi*r^2*total_length
        "l": l_16, 
        "Iz": Iz_16, #1/4*m_8*r**2 + 7/6*m_8*l_8**2,
        "c": 3,
        "d": 0.15,
        "gear": gear
    }
    params_dict = {
        "2": params_2_links,
        "4": params_4_links,
        "8": params_8_links,
        "16": params_16_links,
    }

    return params_dict
def get_parameters():
    exp_da = 2*3.9*10**(-6)   # diameter of active particle in m
    exp_width = 2*10**(-6) # diameter of confinement in m
    exp_rcon = 55*10**(-6)  # radius of ring in m
    exp_Dr = 0.0012         # rotational diffusion coefficient in 1/s
    exp_Dt = 0.8*10**(-14) # translational diffusion coefficient in m^2/s
    exp_visc = 0.004        # viscosity in Pa s

    l0 = exp_da             # length scale in simulations 
    t0 = 1.0/exp_Dr         # time scale in simulations
    
    sim_da = exp_da/l0 # diameter of active particle in simulations
    sim_width = exp_width/l0 # diameter of probe in simulations
    sim_rcon = exp_rcon/l0 # radius of ring in simulations
    sim_Dr = exp_Dr*t0      # rotational diffusion coefficient in simulations
    sim_Dt = exp_Dt*t0/l0**2 # translational diffusion coefficient in simulations

    parameters = {}
    parameters['da'] = sim_da
    parameters['width'] = sim_width
    parameters['rcon'] = sim_rcon
    parameters['Dr'] = sim_Dr
    parameters['Dt'] = sim_Dt
    parameters['l0'] = l0
    parameters['t0'] = t0

    return parameters
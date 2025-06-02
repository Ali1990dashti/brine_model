# Longitudinal and transverse dispersivities
disp_l = 15
disp_t = 5
start_production = 3600

[Mesh]
  type = FileMesh
  file = hex_mesh.e
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 -9.81'
[]

[Variables]
  [pp]
  []
  [C]
  []
[]

[Functions]
  [pres_func]
    type = ParsedFunction
    expression = '-(z-500) * rho * g' # z in mesh is negative
    symbol_names =  'rho g'
    symbol_values = '1000 9.81'
  []
[]

[ICs]
  [pressure_ic]
    type = FunctionIC
    variable = pp
    function = pres_func
  []
  [C_ic]
    type = ConstantIC
    variable = C
    value = 0.01
  []
[]

[AuxVariables]
  [Darcy_vel_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [Darcy_vel_y]
    order = CONSTANT
    family = MONOMIAL
  []
  [Darcy_vel_z]
    order = CONSTANT
    family = MONOMIAL
  []
  [density_water]
    order = CONSTANT
    family = MONOMIAL
  []
  [viscosity_water]
    order = CONSTANT
    family = MONOMIAL
  []
  [xnacl]
    initial_condition = 0.001 # 0.1047 # equivalent to 2.0 molality or 116.88 gr NaCl/1000 gr water
  []
[]

[AuxKernels]
  [Darcy_vel_x]
    type = PorousFlowDarcyVelocityComponent
    variable = Darcy_vel_x
    component = x
    fluid_phase = 0
  []
  [Darcy_vel_y]
    type = PorousFlowDarcyVelocityComponent
    variable = Darcy_vel_y
    component = y
    fluid_phase = 0
  []
  [Darcy_vel_z]
    type = PorousFlowDarcyVelocityComponent
    variable = Darcy_vel_z
    component = z
    fluid_phase = 0
  []
  [density_water]
    type = PorousFlowPropertyAux
    variable = density_water
    property = density
    phase = 0
    execute_on = timestep_end
  []
  [viscosity_water]
    type = PorousFlowPropertyAux
    variable = viscosity_water
    property = viscosity
    phase = 0
    execute_on = timestep_end
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp C'
    number_fluid_phases = 1
    number_fluid_components = 2
  []
[]

[Kernels]
  [mass_der_water]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = pp
  []
  [adv_pp]
    type = PorousFlowFullySaturatedDarcyFlow
    variable = pp
    fluid_component = 1
  []
  [diff_pp]
    type = PorousFlowDispersiveFlux
    fluid_component = 1
    variable = pp
    disp_trans = ${disp_t}
    disp_long = ${disp_l}
  []
  [mass_der_C]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = C
  []
  [adv_C]
    type = PorousFlowFullySaturatedDarcyFlow
    fluid_component = 0
    variable = C
  []
  [diff_C]
    type = PorousFlowDispersiveFlux
    fluid_component = 0
    variable = C
    disp_trans = ${disp_t}
    disp_long = ${disp_l}
  []
[]

[FluidProperties]
  [water]
    type = Water97FluidProperties
  []
  [watertab]
    type = TabulatedBicubicFluidProperties
    fp = water
    save_file = false
  []
[]

[Materials]
  [brine]
    type = PorousFlowBrine
    compute_enthalpy = false
    compute_internal_energy = false
    xnacl = xnacl
    phase = 0
    water_fp = watertab
  []
  [ps]
    type = PorousFlow1PhaseFullySaturated
    porepressure = pp
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.25
  []
  [permeability_bot]
    type = PorousFlowPermeabilityConst
    permeability = '1E-13 0 0   0 1E-13 0   0 0 1E-13'
    block = 'bot_vol'
  []
  [permeability_top]
    type = PorousFlowPermeabilityConst
    permeability = '1E-14 0 0   0 1E-14 0   0 0 1E-14'
    block = 'top_vol'
  []
  # [water]
  #   type = PorousFlowSingleComponentFluid
  #   fp = water
  #   phase = 0
  # []
  [massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = C
  []
  [temperature]
    type = PorousFlowTemperature
    temperature = 303.15 # 30 C
  []
  [diff]
    type = PorousFlowDiffusivityConst
    diffusion_coeff = '0 0'
    tortuosity = 0.1
  []
  [relperm]
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 0
  []
[]

[DiracKernels]
  [c2]
    type = PorousFlowSquarePulsePointSource
    point = '0 0 -10'
    mass_flux = -0.001
    variable = C
    start_time = ${start_production}
  []
  [c3]
    type = PorousFlowSquarePulsePointSource
    point = '0 0 -11'
    mass_flux = -0.001
    variable = C
    start_time = ${start_production}
  []
  [c4]
    type = PorousFlowSquarePulsePointSource
    point = '0 0 -12'
    mass_flux = -0.001
    variable = C
    start_time = ${start_production}
  []
  [c5]
    type = PorousFlowSquarePulsePointSource
    point = '0 0 -13'
    mass_flux = -0.001
    variable = C
    start_time = ${start_production}
  []
  [s2]
    type = PorousFlowSquarePulsePointSource
    point = '0 0 -10'
    mass_flux = -0.1
    variable = pp
    start_time = ${start_production}
  []
  [s3]
    type = PorousFlowSquarePulsePointSource
    point = '0 0 -11'
    mass_flux = -0.1
    variable = pp
    start_time = ${start_production}
  []
  [s4]
    type = PorousFlowSquarePulsePointSource
    point = '0 0 -12'
    mass_flux = -0.1
    variable = pp
    start_time = ${start_production}
  []
  [s5]
    type = PorousFlowSquarePulsePointSource
    point = '0 0 -13'
    mass_flux = -0.1
    variable = pp
    start_time = ${start_production}
  []
[]

[BCs]
  [pressure]
    type = FunctionDirichletBC
    variable = pp
    function = pres_func
    boundary = 'east west'
  []
  [inlet_tracer]
    type = PorousFlowSink
    variable = C
    boundary = 'east'
    flux_function = -0.001
  []
[]

[Preconditioning]
  active = 'basic'
  [basic]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = ' asm      lu           NONZERO                   2'
  []
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = 'gmres bjacobi lu NONZERO'
  []
  [mumps]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  []
[]

[Executioner]
  type = Transient
  end_time = 8640000
  dtmax = 43200
  automatic_scaling = true
  nl_abs_tol = 2e-8
  nl_rel_tol = 1e-5
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 10
  []
[]

[Outputs]
  [exodus]
    type = Exodus
  []
[]

# [Postprocessors]
#   [Li_prod]
#     type = PointValue
#     point = '0 0 -11'
#     variable = C
#   []
#   [pp_prod]
#     type = PointValue
#     point = '0 0 -11'
#     variable = pp
#   []
# []
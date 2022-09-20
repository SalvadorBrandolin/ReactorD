class ReactorBase:
    """ReactorBase class wich the other reactor classes must inherit

    Parameters
    ----------
    mix : Mix object
        Object builded from Mix class
    kinetic : Kinetic object
        Object builded from Kinetic class
        t_operation : str
            Thermal operation of the reactor. Options avaliable: 
            'isothermal', 'not-isothermal'
        p_operation : str
            Thermal operation of the reactor. Options avaliable: 
            'isobaric', 'not-isobaric'
        r_v : float, optional
            Reactor's volume [m^3], by default None
        r_dims_minmax : array-like object[float], optional
            List or narray containing the upper and lower bounds of the 
            spatial variables dominium. E.g. a two dimensions reactor: 
            r_dims_minmax = [[0,1],[0,1]], by default None
        time_minmax : array-like object[flaot], optional
            List or narray containing the upper and lower bounds of
            time dominium, by default None
        f_in : array-like object, optional
            List or array containing the constant or time dependent 
            [func(time)] inlet molar flux of each substance in mix 
            [mol/s]. The unknown inlet molar fluxes must be specified 
            as : 'var', by default None
        f_out : array-like object, optional
            List or array containing the constant or time dependent 
            [func(time)] outlet molar flux of each substance in mix 
            [mol/s]. The unknown outlet molar fluxes must be specified 
            as : 'var', by default None
        t_in : optional
            Temperature of the inlet molar flux [K]. Specify as: 'var'
            if is unknown, by default None
        t_out : optional
            Temperature of the outlet molar flux [K]. Specify as; 'var' 
            if is unknown, by default None
        n_initial : array-like object, optional
            List or array containing the mole numbers of each Substance 
            in mix at the lower bound of time_minmax [mol], The unknown 
            mole amounts must be specified as : 'var', by default None
        t_initial : optional
            Reactor's temperature at the lower bound of time_minmax [K], 
            Specify as; 'var' if is unknown, by default None
        p_initial : optional
            Reactor's pressure at the lower bound of time_minmax [Pa], 
            Specify as; 'var' if is unknown, by default None
        n_final : array-like object, optional
            List or array containing the mole numbers of each Substance 
            in mix at the upper bound of time_minmax [mol], The unknown 
            mole amounts must be specified as : 'var', by default None
        t_final : optional
            Reactor's temperature at the upper bound of time_minmax [K], 
            Specify as; 'var' if is unknown, by default None
        p_final : _type_, optional
            Reactor's pressure at the upper bound of time_minmax [Pa], 
            Specify as; 'var' if is unknown, by default Non
        catalyst_particle : Particle object, optional
            Object builded from Particle class, by default None
        q : float, optional
            Heat flow of the reactor [J/s], by default 0.0
        """
    
    def __init__(
        self, mix, kinetic, t_operation, p_operation,
        r_v=None, 
        r_dims_minmax=None, time_minmax=None,
        f_in=None,
        f_out=None,
        t_in=None, t_out=None, 
        n_initial=None, t_initial=None, p_initial=None,
        n_final=None, t_final=None, p_final=None,
        catalyst_particle=None,
        q=0.0
    ):
        
        self._mix = mix
        self._kinetic = kinetic
        
        # Checking supported thermal, pressure 
        if t_operation in ['isothermal', 'non-isothermal']:
            self._thermal_operation = t_operation
        else:
            raise Exception(f"Not supported thermal operation: {t_operation}")

        if p_operation in ['isobaric', 'non-isobaric']:
            self._pressure_operation = p_operation
        else:
            raise Exception(f"Not supported pressure operation: {p_operation}")

        self._catalyst_particle = catalyst_particle
        self._r_volume = r_v
        self._r_dims_minmax = r_dims_minmax
        self._time_minmax = time_minmax
        self._f_in = f_in
        self._f_out = f_out
        self._t_in = t_in
        self._t_out = t_out
        self._n_initial = n_initial
        self._t_initial = t_initial
        self._p_initial = p_initial
        self._n_final = n_final
        self._t_final = t_final
        self._p_final = p_final
        self._q = q

        #self._degree_of_fredoom_check()
                
    def _degree_of_fredoom_check(self):
        """Method to guarantee that the minimum information is given
        to the Reactor instantiation to solve mass, pressure and energy 
        balances
        """

        if (self._r_dims_minmax or self._r_volume) is not None:
            for i,(fin, fout) in enumerate(zip(self._f_in, self._f_out)):
                if fin == 'var' and fout == 'var':
                    raise Exception("Not border conditions given for "
                                    f"{self._mix.substances[i].name}"
                                    "'s molar flow")
        
        if self._thermal_operation == 'non-isothermal':
            if self._t_in == 'var' and self._t_out == 'var':
                raise Exception("Not border condition given for in or" 
                                " out temperature")
        else:
            pass

        if self._time_minmax is not None:
            for i,(n0, nf) in enumerate(zip(self._n_initial, self._n_final)):
                if n0 == 'var' and nf == 'var':
                        raise Exception("Not border conditions given for "
                                        "initial or final reactor's molar "
                                        "load of:" 
                                        f"{self._mix.substances[i].name}")
        
        if self._thermal_operation == 'non-isothermal':
            if self._t_initial == 'var' and self._t_final == 'var':
                raise Exception("Not border condition given for initial or " 
                                "final reactor's temperature")
        
    def _border_condition_builder(self):
        pass

    def _mass_balance(self):
        pass

    def _energy_balance(self):
        pass

    def _pressure_balance(self):
        pass

    def _particle_mass_balance(self):
        pass
    
    def _particle_energy_balance(self):
        pass
    
    def _result_report(self):
        pass

    def solve(self):
        pass
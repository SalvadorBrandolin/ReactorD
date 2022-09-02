class ReactorBase:
    
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
        self.mix = mix
        self.kinetic = kinetic
        
        # Checking supported thermal, pressure and catalysis operations
        if t_operation in ['isothermal', 'non-isothermal']:
            self._thermal_operation = t_operation
        else:
            raise Exception(f"Not supported thermal operation: {t_operation}")

        if p_operation in ['isobaric', 'non-isobaric']:
            self._pressure_operation = p_operation
        else:
            raise Exception(f"Not supported pressure operation: {p_operation}")

        self.catalyst_particle = catalyst_particle
        self.r_volume = r_v
        self.r_dims_minmax = r_dims_minmax
        self.time_minmax = time_minmax
        self.f_in = f_in
        self.f_out = f_out
        self.t_in = t_in
        self.t_out = t_out
        self.n_initial = n_initial
        self.t_initial = t_initial
        self.p_initial = p_initial
        self.n_final = n_final
        self.t_final = t_final
        self.p_final = p_final
        self.q = q

        self._degree_of_fredoom_check()
                
    def _degree_of_fredoom_check(self):
        if (self.r_dims_minmax or self.r_volume) is not None:
            for i,(fin, fout) in enumerate(zip(self.f_in, self.f_out)):
                if fin == 'var' and fout == 'var':
                    raise Exception("Not border conditions given for "
                                    f"{self.mix[i].name}'s molar flow")
        
        if self._thermal_operation == 'non-isothermal':
            if self.t_in == 'var' and self.t_out == 'var':
                raise Exception("Not border condition given for in or out " 
                                 "temperature")
        else:
            pass

        if self.time_minmax is not None:
            for i,(n0, nf) in enumerate(zip(self.n_initial, self.n_final)):
                if n0 == 'var' and nf == 'var':
                        raise Exception("Not border conditions given for "
                                        "initial or final reactor's molar "
                                        f"load of: {self.mix[i].name}")
        
        if self._thermal_operation == 'non-isothermal':
            if self.t_initial == 'var' and self.t_final == 'var':
                raise Exception("Not border condition given for initial or " 
                                "final reactor's temperature")
        else:
            pass
        pass

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
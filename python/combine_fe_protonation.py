import numpy as np


class BindingFE(object):
    @property
    def temperature(self):
        """float: Temperature of the simulation."""
        return self._temperature

    @temperature.setter
    def temperature(self, value):
        self._temperature = value

    @property
    def mc_steps(self):
        """int: Number of Monte Carlo steps for bootstrapping."""
        return self._mc_steps

    @mc_steps.setter
    def mc_steps(self, value):
        self._mc_steps = value

    @property
    def protonated(self):
        """list: The binding free-energy of the protonated ligand, [fe, sem]"""
        return self._protonated

    @protonated.setter
    def protonated(self, value):
        self._protonated = value

    @property
    def unprotonated(self):
        """list: The binding free-energy of the unprotonated ligand, [fe, sem]"""
        return self._unprotonated

    @unprotonated.setter
    def unprotonated(self, value):
        self._unprotonated = value

    @property
    def results(self):
        """dict: Results of the combined binding free-energy in dictionary"""
        return self._results

    @results.setter
    def results(self, value):
        self._results = value

    @property
    def pKa(self):
        """float: pKa value of the ligand."""
        return self._pKa

    @pKa.setter
    def pKa(self, value):
        self._pKa = value

    @property
    def pH(self):
        """float: pH value of the ligand."""
        return self._pH

    @pH.setter
    def pH(self, value):
        self._pH = value

    def __init__(self, temperature=298.15, mc_steps=1000000):
        self._temperature = temperature
        self._mc_steps = mc_steps
        self.kT = 8.314/4184 * self.temperature
        self.beta = 1.0/self.kT
        self.factor = 0

        self._pKa = None
        self._pH = None

        self._protonated = []
        self._unprotonated = []
        self._results = {}

    def _initialize(self):
        self.kT = 8.314 / 4184 * self.temperature
        self.beta = 1 / self.kT
        assert self.pKa is not None
        assert self.pH is not None
        assert len(self._protonated) > 0
        assert len(self._unprotonated) > 0

        self.factor = 10**(self.pH - self.pKa)

    def bootstrap(self):
        self._initialize()

        print("Running bootstrap on free-energy results")
        cycle_results = np.empty(self.mc_steps)

        for n in range(self._mc_steps):

            # Dissociation constant for protonated
            Kd1 = np.exp(
                -self.beta * np.random.normal(
                    loc=self.protonated[0],
                    scale=self.protonated[1],
                )
            )

            # Dissociation constant for unprotonated
            Kd2 = np.exp(
                -self.beta * np.random.normal(
                    loc=self.unprotonated[0],
                    scale=self.unprotonated[1],
                )
            )

            # Free-energy
            cycle_results[n] = -self.kT * np.log(
                (Kd1 + Kd2*self.factor) / (1 + self.factor)
            )

        dGave = np.mean(cycle_results)
        dGsem = np.std(cycle_results)

        self.results['fe'] = dGave
        self.results['sem'] = dGsem


if __name__ == '__main__':
    calc = BindingFE()
    calc.temperature = 298.15
    calc.mc_steps = 1000
    calc.pKa = 7.5
    calc.pH = 3.5
    calc.protonated = [-13.51, 0.46]
    calc.unprotonated = [-12.34, 0.83]
    calc.bootstrap()
    results = calc.results
    print(f"Combined free-energy: {results['fe']:.2f} +- {results['sem']:.2f} kcal/mol")

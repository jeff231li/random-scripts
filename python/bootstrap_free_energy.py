import numpy as np


class bootstrap(object):
    @property
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self, value):
        self._temperature = value

    @property
    def mc_steps(self):
        return self._mc_steps

    @mc_steps.setter
    def mc_steps(self, value):
        self._mc_steps = value

    @property
    def free_energy(self):
        return self._free_energy

    @free_energy.setter
    def free_energy(self, value):
        self._free_energy = value

    @property
    def enthalpy(self):
        return self._enthalpy

    @enthalpy.setter
    def enthalpy(self, value):
        self._enthalpy = value

    @property
    def results(self):
        return self._results

    @results.setter
    def results(self, value):
        self._results = value

    def __init__(self, temperature=298.15, mc_steps=1000000):
        self._temperature = temperature
        self._mc_steps = mc_steps
        self.kT = 8.314/4184 * self.temperature
        self.beta = 1.0/self.kT
        self.n_orientations = 0

        self._free_energy = []
        self._enthalpy = []
        self._results = {'dG': {}, 'dH': {}}

    def initialize(self):
        self.kT = 8.314 / 4184 * self.temperature
        self.beta = 1 / self.kT
        assert len(self._free_energy) > 0
        self.n_orientations = len(self._free_energy)

        print('Binding free energy')
        for n in range(self.n_orientations):
            self.results['dG'][f'orient{n+1}'] = self._free_energy[n]
            print(f"orient{n+1}: {self._free_energy[n][0]:.2f} +- {self._free_energy[n][1]:.2f}")

        if len(self._enthalpy) > 0:
            assert len(self._free_energy) == len(self._enthalpy)
            print('Binding enthalpy')
            for n in range(self.n_orientations):
                self.results['dH'][f'orient{n + 1}'] = self._enthalpy[n]
                print(f"orient{n + 1}: {self._enthalpy[n][0]:.2f} +- {self._enthalpy[n][1]:.2f}")

    def sampleDG(self):
        self.initialize()

        print("Running bootstrap on free-energy results")
        cycle_results = np.empty(self.mc_steps)

        for n in range(self._mc_steps):

            probability = np.empty(self.n_orientations)

            for orient in range(self.n_orientations):
                probability[orient] = np.exp(
                    -self.beta * np.random.normal(
                        loc=self.free_energy[orient][0],
                        scale=self.free_energy[orient][1],
                    )
                )
            cycle_results[n] = -self.kT*np.log(np.sum(probability))

        dGave = np.mean(cycle_results)
        dGsem = np.std(cycle_results)
        self.results['dG']['value'] = [dGave, dGsem]
        print(f"Combined dG: {dGave:.2f} +- {dGsem:.2f}")

    def sampleDH(self):

        print("Running bootstrap on enthalpy results")
        cycle_results = np.empty(self.mc_steps)

        for n in range(self._mc_steps):

            enthalpy = np.empty(self.n_orientations)
            probability = np.empty(self.n_orientations)

            for orient in range(self.n_orientations):
                enthalpy[orient] = np.random.normal(
                    loc=self.enthalpy[orient][0],
                    scale=self.enthalpy[orient][1],
                )
                probability[orient] = np.exp(
                    -self.beta * np.random.normal(
                        loc=self._free_energy[orient][0],
                        scale=self._free_energy[orient][1],
                    )
                )
            cycle_results[n] = np.sum(enthalpy*probability) / np.sum(probability)

        dHave = np.mean(cycle_results)
        dHsem = np.std(cycle_results)
        self.results['dH']['value'] = [dHave, dHsem]
        print(f"Combined dH: {dHave:.2f} +- {dHsem:.2f}")


if __name__ == '__main__':
    calc = bootstrap()
    calc.free_energy = [[-3.71, 0.27], [-4.22, 0.26]]
    calc.enthalpy = [[-2.21, 0.95], [-1.07, 0.36]]
    calc.sampleDG()
    calc.sampleDH()
    # results = calc.results

import os
import numpy as np
import matplotlib.pyplot as plt
import simtk.unit as unit
from glob import glob


class Work:
    @property
    def path(self) -> str:
        """Path to the trajectory files."""
        return self._path

    @path.setter
    def path(self, value: str):
        self._path = value

    @property
    def colvars_trajectory(self) -> str:
        """File name of the colvars trajectory files."""
        return self._colvars_trajectory

    @colvars_trajectory.setter
    def colvars_trajectory(self, value: str):
        self._colvars_trajectory = value

    @property
    def velocity(self) -> unit:
        """Velocity used in the steered-MD simulation."""
        return self._velocity

    @velocity.setter
    def velocity(self, value: unit):
        self._velocity = value

    @property
    def time_step(self) -> unit:
        """Time step between frames."""
        return self._time_step

    @time_step.setter
    def time_step(self, value: unit):
        self.time_step = value

    @property
    def temperature(self) -> unit:
        """Temperature of the simulation."""
        return self._temperature

    @temperature.setter
    def temperature(self, value: unit):
        self._temperature = value

    @property
    def reaction_coordinate(self) -> list:
        """Reaction coordinate values."""
        return self._reaction_coordinate

    @reaction_coordinate.setter
    def reaction_coordinate(self, value: list):
        self._reaction_coordinate = value

    @property
    def force_column(self) -> int:
        """Column for the force values (output applied force)."""
        return self._force_column

    @force_column.setter
    def force_column(self, value: int):
        self._force_column = value

    @property
    def position_column(self) -> int:
        """Column for the force values (output applied force)."""
        return self._position_column

    @position_column.setter
    def position_column(self, value: int):
        self._position_column = value

    @property
    def results(self) -> dict:
        """Final results in a python dictionary."""
        return self._results

    @results.setter
    def results(self, value: dict):
        self._results = value

    def __init__(self):
        self._velocity = None
        self._time_step = None
        self._temperature = None
        self._reaction_coordinate = None
        self._path = "./"
        self._colvars_trajectory = None
        self._position_column = 1
        self._force_column = 2
        self.data = []
        self._results = {}
        self.n_ensemble = 0
        self.n_data = []
        self.work = []

    def _read_file(self, file):
        data = []
        with open(os.path.join(self.path, file), "r") as f:
            for line in f.readlines():
                if not line.startswith("#"):
                    data.append(
                        [
                            float(line.split()[self.position_column]),
                            float(line.split()[self.force_column]),
                        ]
                    )

        return data

    def collect_data(self):
        file_names = glob(self.colvars_trajectory)
        self.n_ensemble = len(file_names)
        for file in file_names:
            data = self._read_file(file)
            self.n_data.append(len(data))
            self.data.append(data)
            print(f"{file}: {self.n_data} data points")

    def calculate_smd_work(self):

        vel_times_dt = self.velocity.value_in_unit(unit.angstrom / unit.nanosecond) *\
                       self.time_step.value_in_unit(unit.nanosecond)

        self.results['work'] = np.zeros((self.n_ensemble, self.n_data[0]))

        for n in range(self.n_ensemble):
            force = self.data[n]
            work = np.zeros(self.n_data[n])

            for i in range(self.n_data[n]):
                if i == 0:
                    self.results['work'][n, i] = force[i] * vel_times_dt
                else:
                    self.results['work'][n, i] = work[i-1] + force[i] * vel_times_dt

    def calculate_pmf_ensemble(self):
        kT = (self.temperature * unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA).value_in_unit(unit.kilocalorie_per_mole)
        beta = 1.0 / kT

        probability = np.exp(-beta * self.results['work'])

        self.results['pmf'] = kT*np.log(np.mean(probability, axis=0))

    def plot_ensemble(self):
        pass


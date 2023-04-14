from pyiron_contrib.tinybase.container import (
            AbstractOutput,
            StructureInput,
)
from pyiron_contrib.tinybase.task import (
            AbstractTask,
            ListTask,
            ListInput,
            ReturnStatus
)

from copy import deepcopy

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si
import scipy.optimize as so

from pyiron_atomistics.atomistics.structure.has_structure import HasStructure

MurnaghanInputBase = StructureInput.from_attributes(
        "MurnaghanInputBase",
        "strains",
        "task"
)

class MurnaghanInput(MurnaghanInputBase, ListInput):
    def check_ready(self):
        structure_ready = self.structure is not None
        strain_ready = len(self.strains) > 0
        task = self.task
        task.input.structure = self.structure
        return structure_ready and strain_ready and task.input.check_ready()

    def set_strain_range(self, range, steps):
        self.strains = (1 + np.linspace(-range, range, steps))**(1/3)

    def _create_tasks(self):
        cell = self.structure.get_cell()
        tasks = []
        for s in self.strains:
            n = deepcopy(self.task)
            n.input.structure = self.structure.copy()
            n.input.structure.set_cell(cell * s, scale_atoms=True)
            tasks.append(n)
        return tasks

MurnaghanOutputBase = AbstractOutput.from_attributes(
        "MurnaghanOutputBase",
        'base_structure',
        volumes=list,
        energies=list
)

class MurnaghanOutput(MurnaghanOutputBase, HasStructure):
    def plot(self):
        plt.plot(self.volumes, self.energies)

    @property
    def equilibrium_volume(self):
        inter = si.interp1d(self.volumes, self.energies)
        return so.minimize_scalar(inter, bounds=(np.min(self.volumes), np.max(self.volumes))).x

    def _number_of_structures(self):
        return 1

    def _get_structure(self, frame, wrap_atoms=True):
        s = self.base_structure
        s.set_cell(s.get_cell() * (self.equilibrium_volume/s.get_volume())**(1/3))
        return s

class MurnaghanTask(ListTask):

    def _get_input(self):
        return MurnaghanInput()

    def _get_output(self):
        return MurnaghanOutput()

    def _execute(self, output):
        output.base_structure = self.input.structure
        return super()._execute(output)

    def _extract_output(self, output, step, node, ret, node_output):
        if len(output.energies) == 0:
            output.energies = np.zeros(len(self.input.strains))
        if len(output.volumes) == 0:
            output.volumes = np.zeros(len(self.input.strains))
        if ret.is_done():
            output.energies[step] = task_output.energy_pot
            output.volumes[step] = task.input.structure.get_volume()

# Capsul import
from capsul.api import Process

# Trait import
from traits.api import Float


class AvailableProcesses(list):
    # List of all the processes in this file
    def __init__(self):
        super(AvailableProcesses, self).__init__()
        self.append(Addition)
        self.append(Substraction)


class Addition(Process):

    def __init__(self, node_name):
        super(Addition, self).__init__()

        self.node_name = node_name

        self.add_trait(node_name + "_in_1", Float(output=False))
        self.add_trait(node_name + "_in_2", Float(output=False))
        self.add_trait(node_name + "_out", Float(output=True))

        self.in_1 = getattr(self, node_name + "_in_1")
        self.in_2 = getattr(self, node_name + "_in_2")
        self.out = getattr(self, node_name + "_out")

    def _run_process(self):
        self.get_plugs_values()

        self.out = self.in_1 + self.in_2
        print('Addition - ', self.node_name, '\n...\nInputs: {', self.in_1, ', ',
              self.in_2, '}\nOutput: ', self.out, '\n...\n')

        self.set_plugs_values()

    def get_plugs_values(self):
        self.in_1 = getattr(self, self.node_name + "_in_1")
        self.in_2 = getattr(self, self.node_name + "_in_2")
        self.out = getattr(self, self.node_name + "_out")

    def set_plugs_values(self):
        setattr(self, self.node_name + "_in_1", self.in_1)
        setattr(self, self.node_name + "_in_2", self.in_2)
        setattr(self, self.node_name + "_out", self.out)


class Substraction(Process):

    def __init__(self, node_name):
        super(Substraction, self).__init__()

        self.node_name = node_name

        self.add_trait(node_name + "_in_1", Float(output=False))
        self.add_trait(node_name + "_in_2", Float(output=False))
        self.add_trait(node_name + "_out", Float(output=True))

        self.in_1 = getattr(self, node_name + "_in_1")
        self.in_2 = getattr(self, node_name + "_in_2")
        self.out = getattr(self, node_name + "_out")

    def _run_process(self):
        self.get_plugs_values()

        self.out = self.in_1 - self.in_2
        print('Substraction - ', self.node_name, '\n...\nInputs: {', self.in_1, ', ',
              self.in_2, '}\nOutput: ', self.out, '\n...\n')

        self.set_plugs_values()

    def get_plugs_values(self):
        self.in_1 = getattr(self, self.node_name + "_in_1")
        self.in_2 = getattr(self, self.node_name + "_in_2")
        self.out = getattr(self, self.node_name + "_out")

    def set_plugs_values(self):
        setattr(self, self.node_name + "_in_1", self.in_1)
        setattr(self, self.node_name + "_in_2", self.in_2)
        setattr(self, self.node_name + "_out", self.out)

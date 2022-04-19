import numpy as np

from names import Actions as act
from names import ParamNames as pn
from names import GlobNames as gn
from names import PropNames as prn
from model import Model, ModelFactory

class SampleModel(Model):
    def take_action(self, action, params, globdat):
        print('BarModel taking action', action)

    def configure(self, props, globdat):
        pass

def declare(factory):
    factory.declare_model('Sample', SampleModel)

import numpy as np

from model import *

MODELS = 'models'
TYPE = 'type'


class MultiModel(Model):
    def take_action(self, action, params, globdat):
        for model in self._models:
            model.take_action(action, params, globdat)

    def configure(self, props, globdat):
        models = props[MODELS].strip('[').strip(']').split(',')

        mfac = globdat['modelFactory']

        self._models = []

        for m in models:
            modelprops = props[m]

            print('multimodel ', modelprops[TYPE])
            model = mfac.get_model(modelprops[TYPE], m)
            model.configure(modelprops, globdat)

            self._models.append(model)


def declare(factory):
    factory.declare_model('Multi', MultiModel)

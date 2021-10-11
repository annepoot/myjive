import numpy as np

from model import *

class BarModel (Model):
    def take_action (self, action, params, globdat):
        print('BarModel takes action', action)

        if 'stiffness' in action:
            self.__stiffness(params, globdat)

    def configure (self, props, globdat):
        print('BarModel::configure')

    def __stiffness (self, params, globdat):
        print('Someday this will compute K... if they let me work in peace')


def declare (factory):
    factory.declare_model ('Bar', BarModel)


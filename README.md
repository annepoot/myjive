# MyJive
A personal adaptation of [PyJive](https://gitlab.tudelft.nl/fmeer/pyjive) that is set up in a more "Pythonic" way.

## Getting started
Run the following anaconda commands to take care of all dependencies:

```
conda env create -f ENVIRONMENT.yml
conda activate myjive
```

## Differences with PyJive
PyJive is a Python adaptation of the C++ Jive library that stays very close to C++ setup.
In MyJive, more liberty is taken, since Python offers a lot more flexibility that C++ does not have.

Key differences are:

- `MultiModel()` has been removed. In (Py)Jive, each module can only have a single model associated with it. The `MultiModel` class was created to connect multiple models to a single module. In MyJive, each module has a list of models associated with it, rather than a single model, so the `MultiModel` class is redundant.
- `take_action(action)` has been removed. In (Py)Jive, modules use the `take_action` function to determine which model functions are called. The `take_action` function is called on all child models, and based on the value of the `action` parameter, each child model decides whether it should do something or not. To make it less opaque what models come into play at which actions, MyJive has been set up so that modules collect all relevant models that are able to perform a certain action, and then call the action only on these models. All methods with an all-caps name are considered actions.
- `params` has been fully removed. All input and output and output has now been made explicit for all actions that are executed by the modules in the models. Optional arguments are handled with `**kwargs` arguments.
- Individual `declare_model` and `declare_module` functions have been replaced with the `declare` class method in the `Module` and `Model` classes.
- Jive and non-Jive functionality has been separated into different folders (`/jive/` and `/core/`, respectively).

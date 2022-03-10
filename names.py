class GlobNames:
    MODEL = 'model'
    DOFSPACE = 'dofSpace'
    MESHSHAPE = 'meshShape'
    MESHRANK = 'meshRank'
    ESET = 'elemSet'
    NSET = 'nodeSet'
    NGROUPS = 'nodeGroups'
    EGROUPS = 'elemGroups'
    MODELFACTORY = 'modelFactory'
    MODULEFACTORY = 'moduleFactory'
    SHAPEFACTORY = 'shapeFactory'
    TIMESTEP = 'timeStep'
    STATE0 = 'state0'
    MATRIX0 = 'matrix0'
    TABLES = 'tables'


class PropNames:
    TYPE = 'type'


class Actions:
    GETMATRIX0 = 'getMatrix0'
    GETMATRIXLB = 'getMatrixLB'
    GETINTFORCE = 'getIntForce'
    GETEXTFORCE = 'getExtForce'
    GETCONSTRAINTS = 'getConstraints'
    GETTABLE = 'getTable'


class ParamNames:
    MATRIX0 = 'matrix0'
    MATRIX1 = 'matrix1'
    INTFORCE = 'intForce'
    EXTFORCE = 'extForce'
    CONSTRAINTS = 'constraints'
    TABLE = 'table'
    TABLENAME = 'tableName'
    TABLEWEIGHTS = 'tableWeights'

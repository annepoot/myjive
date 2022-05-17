class GlobNames:
    MODEL = 'model'
    DOFSPACE = 'dofSpace'
    MESHSHAPE = 'meshShape'
    MESHRANK = 'meshRank'
    COARSEMESH = 'coarseMesh'
    ESET = 'elemSet'
    NSET = 'nodeSet'
    NGROUPS = 'nodeGroups'
    EGROUPS = 'elemGroups'
    MODELFACTORY = 'modelFactory'
    MODULEFACTORY = 'moduleFactory'
    SHAPEFACTORY = 'shapeFactory'
    TIMESTEP = 'timeStep'
    STATE0 = 'state0'
    EXTFORCE = 'extForce'
    CONSTRAINTS = 'constraints'
    OLDSTATE0 = 'oldstate0'
    BACKUPSTATE0 = 'backupstate0'
    HISTORY = 'history'
    MATRIX0 = 'matrix0'
    MATRIX2 = 'matrix2'
    TABLES = 'tables'
    LAMBDA = 'lambda'
    MODULEDATA = 'module'
    SLIDERS = 'sliders'
    ACCEPTED = 'accepted'
    LBFACTORS = 'lbFactors'


class PropNames:
    TYPE = 'type'


class Actions:
    GETMATRIX0 = 'getMatrix0'
    GETMATRIX2 = 'getMatrix2'
    GETMATRIXLB = 'getMatrixLB'
    GETINTFORCE = 'getIntForce'
    GETEXTFORCE = 'getExtForce'
    GETUNITFORCE = 'getUnitForce'
    CHECKCOMMIT = 'checkCommit'
    GETCONSTRAINTS = 'getConstraints'
    GETTABLE = 'getTable'
    ADVANCE = 'advance'


class ParamNames:
    MATRIX0 = 'matrix0'
    MATRIX1 = 'matrix1'
    MATRIX2 = 'matrix2'
    INTFORCE = 'intForce'
    EXTFORCE = 'extForce'
    UNITFORCE = 'unitForce'
    CONSTRAINTS = 'constraints'
    TABLE = 'table'
    TABLENAME = 'tableName'
    TABLEWEIGHTS = 'tableWeights'


class GPActions:
    GETPRIORMEAN = 'getPriorMean'
    GETPRIORCOVARIANCE = 'getPriorCovariance'
    GETPOSTERIORMEAN = 'getPosteriorMean'
    GETPOSTERIORCOVARIANCE = 'getPosteriorCovariance'
    GETPRIORSAMPLES = 'getPriorSamples'
    GETPOSTERIORSAMPLES = 'getPosteriorSamples'


class GPParamNames:
    PRIORMEAN = 'priorMean'
    PRIORCOVARIANCE = 'priorCovariance'
    POSTERIORMEAN = 'posteriorMean'
    POSTERIORCOVARIANCE = 'posteriorCovariance'
    PRIORSAMPLES = 'priorSamples'
    POSTERIORSAMPLES = 'posteriorSamples'
    FIELD = 'field'
    FULLCOVARIANCE = 'fullCovariance'
    NSAMPLE = 'nsamples'

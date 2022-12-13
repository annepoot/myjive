gpinit =
{
  type = GPInit;

  mesh =
  {
    type = gmsh;
    file = tri6mesh.msh;
  };

  coarseMesh =
  {
    type = gmsh;
    file = tri3mesh.msh;
  };

  nodeGroups = [ left, right, bottom ];

  left =
  {
    xtype = min;
  };

  right =
  {
    xtype = max;
  };

  bottom =
  {
    ytype = min;
  };
};

gpsolver =
{
  type = GPSolver;
  solver =
  {
    type = CG;
    maxIter = 100;
    allowMaxIter = False;
  };
  preconditioner =
  {
    type = ichol;
  };
  nsample = 3;
  seed = 0;
  storeMatrix = True;
  storeConstraints = True;
  getUnitMassMatrix = True;
  getForceResults = True;
};

model =
{
  type = Multi;

  models = [ poisson, gp, diri ];

  poisson =
  {
    type = XPoisson;

    elements = all;

    kappa = 1.0;
    rho = 1;
    q = (x<1)*1+(x>1)*-0.8;

    shape =
    {
      type = Triangle6;
      intScheme = Gauss3;
    };
  };


  gp =
  {
    type = GPf;

    prior =
    {
      type = SPDE;
      func = alpha**2 * M;
      hyperparams =
      {
        alpha = opt;
      };
    };

    obsNoise = 1e-5;

    shape =
    {
      type = Triangle3;
      intScheme = Gauss1;
    };
  };

  diri =
  {
    type = Dirichlet;

    groups = [ bottom ];
    dofs   = [ u ];
    values = [ 0.0 ];
  };

  neum =
  {
    type = Neumann;

    groups = [ right ];
    dofs   = [ u ];
    values = [ 0.4 ];
  };
};

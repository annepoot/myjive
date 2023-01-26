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
  type = GPSampler;
  solver =
  {
    type = cholmod;
  };
  nsample = 100;
  seed = 0;
};

model =
{
  type = Multi;

  models = [ poisson, gp, load, diri ];

  poisson =
  {
    type = XPoisson;

    elements = all;

    kappa = 1.0;

    shape =
    {
      type = Triangle6;
      intScheme = Gauss3;
    };
  };

  gp =
  {
    type = GPf;

    explicitInverse = False;
    coarseInit = False;

    prior =
    {
      type = SPDE;
      func = alpha**2 * M;
      hyperparams =
      {
        alpha = 0.3;
      };
    };

    solver =
    {
      type = CG;
      maxIter = 200;
      allowMaxIter = False;
    };

    preconditioner =
    {
      type = diag;
    };

    obsNoise = 1e-5;

    shape =
    {
      type = Triangle3;
      intScheme = Gauss1;
    };
  };

  load =
  {
    type = Load;

    elements = all;

    dofs   = [ u ];
    values = [ (x<1)*1+(x>1)*-0.8 ];

    shape =
    {
      type = Triangle6;
      intScheme = Gauss3;
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

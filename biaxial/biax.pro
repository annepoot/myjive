init =
{
  type = Init;
  model = femodel;

  nodeGroups = [ top, left, bottom, right ];

  mesh =
  {
    type = gmsh;
    file = biax_fine.msh;
  };

  top =
  {
    ytype = max;
  };

  bottom =
  {
    ytype = min;
  };

  left =
  {
    xtype = min;
  };

  right =
  {
    xtype = max;
  };
};

solver =
{
  type = Solver;
  model = femodel;

  nsteps = 1;
  storeMatrix = True;
  storeConstraints = True;
  getUnitMassMatrix = True;
};

femodel =
{
  type = Multi;

  models = [ elastic, diri ];

  elastic =
  {
    type = Elastic;

    elements = all;

    young = 10000.;
    thickness = 0.2;
    poisson = 0.2;
    rho = 0.0;
    state = plane_stress;

    shape =
    {
      type = Triangle3;
      intScheme = Gauss1;
    };
  };

  diri =
  {
    type = Dirichlet;

    groups = [ top, bottom, left, right ];
    dofs   = [ dy, dy, dx, dx ];
    values = [ 1., 0., 0., 1. ];
  };
};

gpinit =
{
  type = GPInit;
  model = gpmodel;

  nodeGroups = [ top, left, bottom, right ];

  mesh =
  {
    type = gmsh;
    file = biax_coarse.msh;
  };

  top =
  {
    ytype = max;
  };

  bottom =
  {
    ytype = min;
  };

  left =
  {
    xtype = min;
  };

  right =
  {
    xtype = max;
  };
};

gaussian =
{
  type = Gaussian;
  model = gpmodel;

  storeMatrix = True;
  getUnitMassMatrix = True;
};

sampler =
{
  type = Sampler;
  model = gpmodel;

  nsample = 3;
  seed = 0;
};

gpmodel =
{
  type = GP;

  obsNoise = 1e-10;
  alpha = 100;

  shape =
  {
    type = Triangle3;
    intScheme = Gauss1;
  };
};

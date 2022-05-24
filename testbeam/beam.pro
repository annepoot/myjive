init =
{
  type = Init;
  model = femodel;

  nodeGroups = [ lb, rb, tm ];

  mesh =
  {
    type = gmsh;
    file = beam.msh;
  };

  lb =
  {
    xtype = min;
    ytype = min;
  };

  rb =
  {
    xtype = max;
    ytype = min;
  };

  tm =
  {
    ytype = max;
    xtype = mid;
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

  models = [ elastic, diri, neum ];

  elastic =
  {
    type = Elastic;

    elements = all;

    young = 10000.;
    thickness = 0.2;
    poisson = 0.2;
    rho = 1;
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

    groups = [ lb, lb, rb ];
    dofs   = [ dx, dy, dy ];
    values = [ 0., 0., 0. ];
  };

  neum =
  {
    type = Neumann;

    groups = [ tm ];
    dofs   = [ dy ];
    values = [ -1. ];
  };
};

gpinit =
{
  type = GPInit;
  model = gpmodel;

  nodeGroups = [ lb, rb, tm ];

  mesh =
  {
    type = gmsh;
    file = beam_coarse.msh;
  };

  lb =
  {
    xtype = min;
    ytype = min;
  };

  rb =
  {
    xtype = max;
    ytype = min;
  };

  tm =
  {
    ytype = max;
    xtype = mid;
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

  nsample = 10;
};

gpmodel =
{
  type = GP;

  obsNoise = 1e-5;
  alpha = opt;

  shape =
  {
    type = Triangle3;
    intScheme = Gauss1;
  };
};

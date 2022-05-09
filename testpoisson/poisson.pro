init =
{
  type = Init;
  model = femodel;

  mesh =
  {
    type = gmsh;
    file = mesh.msh;
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

solver =
{
  type = Solver;
  model = femodel;

  nsteps = 1;
  storeMatrix = True;
  storeConstraints = True;
  getMassMatrix = True;
};

femodel =
{
  type = Multi;

  models = [ poisson, diri, neum ];

  poisson =
  {
    type = Poisson;

    elements = all;

    kappa = 1.0;
    rho = 1;

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

gpinit =
{
  type = GPInit;
  model = gpmodel;

  mesh =
  {
    type = gmsh;
    file = mesh.msh;
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

gaussian =
{
  type = Gaussian;
  model = gpmodel;

  storeMatrix = True;
  getMassMatrix = True;
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
  alpha = 1;

  shape =
  {
    type = Triangle3;
    intScheme = Gauss1;
  };
};

view =
{
  type = View;
  model = gpmodel;

  plot = std_u_post[u];
  ncolors = 100;
};

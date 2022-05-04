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

view =
{
  type = View;
  model = femodel;

  plot = solution[u];
  ncolors = 100;
};

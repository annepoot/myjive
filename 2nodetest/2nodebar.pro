init =
{
  type = Init;
  model = femodel;

  nodeGroups = [ left, right, mid ];

  mesh =
  {
    type = manual;
    file = 2nodebar.mesh;
  };

  left =
  {
    xtype = min;
  };

  right =
  {
    xtype = max;
  };

  mid =
  {
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
  getMassMatrix = True;
};

femodel =
{
  type = Multi;
  models = [ bar, diri, neum ];

  bar =
  {
    type = Bar;

    elements = all;

    EA = 1.0;
    k = 1.0;
    q = 1.0;

    shape =
    {
      type = Line2;
      intScheme = Gauss2;
    };
  };

  diri =
  {
    type = Dirichlet;

    groups = [ left, right ];
    dofs   = [ dx, dx ];
    values = [ 0.0, 0.0 ];
  };

  neum =
  {
    type = Neumann;

    groups = [ mid ];
    dofs = [ dx ];
    values = [ 1.0 ];
  };
};

gaussian =
{
  type = Gaussian;
  model = gpmodel;

  storeMatrix = True;
  getMassMatrix = True;

  nobs = 9;
  obsNoise = 1e-5;
  alpha = 1;
};

gpmodel =
{
  type = GP;

  obsNoise = 1e-5;
  alpha = opt;
  nobs = 9;
}
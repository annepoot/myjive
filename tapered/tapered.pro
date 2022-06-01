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
  getUnitMassMatrix = True;
};

femodel =
{
  type = Multi;
  models = [ bar, diri ];

  bar =
  {
    type = Bar;

    elements = all;

    A = 1.0;
    E = 3.0 - 0.29 * x;
    k = 0.0;
    q = 0.1;

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
    values = [ 0.0, 1.0 ];
  };

  neum =
  {
    type = Neumann;

    groups = [ right ];
    dofs = [ dx ];
    values = [ 1.0 ];
  };
};

gpinit =
{
  type = GPInit;
  model = gpmodel;

  nodeGroups = [ left, right, mid ];

  mesh =
  {
    type = manual;
    file = 2nodebar_coarse.mesh;
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

  nsample = 30;
  seed = 0;
};

gpmodel =
{
  type = GP;

  obsNoise = 1e-10;
  alpha = 0.1;

  shape =
  {
    type = Line2;
    intScheme = Gauss2;
  };
}

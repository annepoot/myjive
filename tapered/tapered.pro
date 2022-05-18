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

    A = 1.0;
    E = 3.0 - 0.58 * abs(x-5);
    k = 0.0;
    q = 0.0;

    shape =
    {
      type = Line2;
      intScheme = Gauss2;
    };
  };

  diri =
  {
    type = Dirichlet;

    groups = [ left ];
    dofs   = [ dx ];
    values = [ 0.0 ];
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
};

gpmodel =
{
  type = GP;

  obsNoise = 1e-5;
  alpha = opt;

  shape =
  {
    type = Line2;
    intScheme = Gauss2;
  };
}

init =
{
  model = femodel;

  nodeGroups = [ left, right, mid ];

  mesh =
  {
    type = manual;
    file = 3nodebar.mesh;
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
  model = femodel;
  nsteps = 1;
  storeMatrix = True;
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
    k = 0.0;
    q = 1.0;

    shape =
    {
      type = Line3;
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
  alpha = opt;
};

gpmodel =
{
  type = Multi;

  models = [ sampler ];

  sampler =
  {
    type = Sample;
  };
}
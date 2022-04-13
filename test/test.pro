init =
{
  nodeGroups = [ left, right, mid ];

  mesh =
  {
    type = manual;
    file = bar.mesh;
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

model =
{
  type = Multi;

  models = [ bar, diri, neum ];

  bar =
  {
    type = Bar;

    elements = all;

    EA = 1.0;
    k = 1.0;

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
  nsteps = 1;
  storeMatrix = True;

  nobs = 21;
  obsNoise = 1e-5;
  alpha = 1;
};
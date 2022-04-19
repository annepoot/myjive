init =
{
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
  storeMatrix = True;

  nobs = 9;
  obsNoise = 1e-5;
  alpha = opt;
};
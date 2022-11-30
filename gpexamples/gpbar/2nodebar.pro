gpinit =
{
  type = GPInit;

  mesh =
  {
    type = manual;
    file = 2nodebar.mesh;
  };

  coarseMesh =
  {
    type = manual;
    file = 2nodebar_coarse.mesh;
  };

  nodeGroups = [ left, right, mid ];

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

gpsolver =
{
  type = GPSolver;
  nsample = 30;
  seed = 0;
  storeMatrix = True;
  storeConstraints = True;
  getUnitMassMatrix = True;
  getForceResults = True;
};

model =
{
  type = Multi;
  models = [ bar, gp, diri, neum ];

  bar =
  {
    type = XBar;

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

  gp =
  {
    type = GPf;

    prior =
    {
      type = SPDE;
      func = alpha**2 * M;
      hyperparams =
      {
        alpha = opt;
      };
    };

    obsNoise = 1e-5;

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
    values = [ 5.0 ];
  };
};

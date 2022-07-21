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

  nsteps = 1;
  storeMatrix = True;
  storeConstraints = True;
  getUnitMassMatrix = True;
  getForceResults = True;
};

gpsampler =
{
  type = GPSampler;

  nsample = 30;
  seed = 0;
};

model =
{
  type = Multi;
  models = [ bar, gp, diri ];

  bar =
  {
    type = XBar;

    elements = all;

    EA = 3.0 - 0.29 * x;
    k = 0.0;
    q = 0.1;

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
        alpha = 0.1;
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
    values = [ 0.0, 1.0 ];
  };
};

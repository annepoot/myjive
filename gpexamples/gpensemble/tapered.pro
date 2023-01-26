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
};

model =
{
  type = Multi;
  models = [ bar, gp, load, diri ];

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
    type = GPEnKF;

    prior =
    {
      type = SPDE;
      func = alpha**2 * M;
      hyperparams =
      {
        alpha = 0.1;
      };
      premultiplier = K;
    };

    ensemble = 1000;
    seed = 0;
    obsNoise = 1e-5;

    shape =
    {
      type = Line2;
      intScheme = Gauss2;
    };
  };

  load =
  {
    type = Load;

    elements = all;

    dofs = [ dx ];
    values = [ 0.1 ];

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

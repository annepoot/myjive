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
  getFullCovariance = True;
};

model =
{
  type = Multi;
  models = [ solid, gp, diri ];

  solid =
  {
    type = Solid;

    elements = all;

    material =
    {
      type = Heterogeneous;
      rank = 1;
      anmodel = bar;

      E = 3.0 - 0.29 * x;
      rho = 0.1;
    };

    gravity = True;

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
    values = [ 0.0, 1.0 ];
  };
};

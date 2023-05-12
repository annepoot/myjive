gpinit =
{
  type = GPInit;

  mesh =
  {
    type = manual;
    file = bar_fine.mesh;
  };

  coarseMesh =
  {
    type = manual;
    file = bar_coarse.mesh;
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
  priorMean = dirichlet;
  nsample = 30;
  seed = 0;
};

model =
{
  type = Multi;
  models = [ solid, gp, load, diri ];

  solid =
  {
    type = Solid;

    elements = all;

    material =
    {
      type = Heterogeneous;
      rank = 1;
      anmodel = bar;

      E = 1.0 - 0.99 * x;
    };

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
        alpha = 1.0;
      };
      mean = dirichlet;
    };

    obsNoise = 1e-8;
    pdNoise = 1e-8;

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

    dofs   = [ dx ];
    values = [ 3. ];

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

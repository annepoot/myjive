gpinit =
{
  type = GPInit;

  mesh =
  {
    type = gmsh;
    file = meshes/bridge-q4-r2.msh;
  };

  coarseMesh =
  {
    type = gmsh;
    file = meshes/bridge-q4-r0.msh;
  };

  nodeGroups = [ lb, rb, m1, m2 ];

  lb =
  {
    xtype = min;
    ytype = min;
  };

  rb =
  {
    xtype = max;
    ytype = min;
  };

  m1 =
  {
    xtype = 3;
    ytype = min;
  };

  m2 =
  {
    xtype = 7;
    ytype = min;
  };
};

gpsolver =
{
  type = GPSolver;
  nsample = 100;
  seed = 0;
  tables = [ stress, strain ];
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
      type = Isotropic;
      rank = 2;
      anmodel = plane_stress;

      E = 10000.;
      nu = 0.2;
    };

    shape =
    {
      type = Quad4;
      intScheme = Gauss4;
    };
  };

  load =
  {
    type = Load;

    elements = all;

    dofs   = [ dy ];
    values = [ -10. ];

    shape =
    {
      type = Quad4;
      intScheme = Gauss4;
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

    obsNoise = 1e-10;

    shape =
    {
      type = Quad4;
      intScheme = Gauss4;
    };
  };

  diri =
  {
    type = Dirichlet;

    groups = [ lb, lb, m1, m2, rb ];
    dofs   = [ dx, dy, dy, dy, dy ];
    values = [ 0., 0., 0., 0., 0. ];
  };
};

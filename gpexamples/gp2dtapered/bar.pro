gpinit =
{
  type = GPInit;

  mesh =
  {
    type = gmsh;
    file = meshes/4-node/bar_fine2.msh;
  };

  coarseMesh =
  {
    type = gmsh;
    file = meshes/4-node/bar_coarse.msh;
  };

  nodeGroups = [ lb, rb ];

  lb =
  {
    xtype = min;
  };

  rb =
  {
    xtype = max;
  };
};

gpsolver =
{
  type = GPSolver;
  nsample = 3;
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
      type = Heterogeneous;
      rank = 2;
      anmodel = plane_stress;

      E = 1.0 - 0.99 * x;
      nu = 0.2;
    };

    thickness = 4.;

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

  load =
  {
    type = Load;

    elements = all;

    dofs   = [ dx ];
    values = [ 12. ];

    shape =
    {
      type = Quad4;
      intScheme = Gauss4;
    };
  };

  diri =
  {
    type = Dirichlet;

    groups = [ lb, lb, rb ];
    dofs   = [ dx, dy, dx ];
    values = [ 0., 0., 1. ];
  };
};

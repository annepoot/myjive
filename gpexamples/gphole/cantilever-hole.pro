gpinit =
{
  type = GPInit;

  mesh =
  {
    type = gmsh;
    file = meshes/q9/hole-9.msh;
  };

  coarseMesh =
  {
    type = gmsh;
    file = meshes/q4/hole-9.msh;
  };

  nodeGroups = [ l, lt, lb, r, rt, rb ];

  l =
  {
    xtype = min;
  };

  lt =
  {
    xtype = min;
    ytype = max;
  };

  lb =
  {
    xtype = min;
    ytype = min;
  };

  r =
  {
    xtype = max;
  };

  rt =
  {
    xtype = max;
    ytype = max;
  };

  rb =
  {
    xtype = max;
    ytype = min;
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

  models = [ solid, gp, load, diri, neum ];

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

    thickness = 0.2;

    shape =
    {
      type = Quad9;
      intScheme = Gauss9;
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
      intScheme = Gauss9;
    };
  };

  load =
  {
    type = Load;

    elements = all;

    dofs   = [ dy ];
    values = [ 0. ];

    shape =
    {
      type = Quad9;
      intScheme = Gauss9;
    };
  };

  diri =
  {
    type = Dirichlet;

    groups = [ l,  lb ];
    dofs   = [ dx, dy ];
    values = [ 0., 0. ];
  };

  neum =
  {
    type = Neumann;

    groups = [ rt ];
    dofs   = [ dy ];
    values = [ 1. ];
  };
};

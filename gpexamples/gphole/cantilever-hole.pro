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

  nodeGroups = [ l, rm ];

  l =
  {
    xtype = min;
  };

  rm =
  {
    xtype = max;
    ytype = mid;
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

  models = [ solid, gp, diri, neum];

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
        alpha = 1.;
      };
    };

    obsNoise = 1e-10;

    shape =
    {
      type = Quad4;
      intScheme = Gauss9;
    };
  };

  diri =
  {
    type = Dirichlet;

    groups = [ l,  l , rm ];
    dofs   = [ dx, dy, dx ];
    values = [ 0., 0., 1. ];
  };

  neum =
  {
    type = Neumann;

    groups = [];
    dofs   = [];
    values = [];
  };
};

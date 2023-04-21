gpinit =
{
  type = GPInit;

  mesh =
  {
    type = gmsh;
    file = meshes-void/nfib-16_r1.msh;
  };

  coarseMesh =
  {
    type = gmsh;
    file = meshes-void/nfib-16_r0.msh;
  };

  nodeGroups = [ l, b, r, t ];

  l =
  {
    xtype = min;
  };

  b =
  {
    ytype = min;
  };

  r =
  {
    xtype = max;
  };

  t =
  {
    ytype = max;
  };
};

gpsolver =
{
  type = GPSolver;
  nsample = 50;
  seed = 0;
  tables = [ strain ];
  priorMean = dirichlet;
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
      type = Isotropic;
      rank = 2;
      anmodel = plane_stress;

      E = 3000.;
      nu = 0.2;
      thickness = 1.;
    };

    shape =
    {
      type = Triangle3;
      intScheme = Gauss1;
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

    obsNoise = 1e-2;
    pdNoise = 0;
    bcNoise = 1e-8;

    boundary =
    {
      type = dirichlet;
      groups = [ r , t ];
      dofs   = [ dx, dy ];
      covs   = [ 0.1, 0.1 ];
    };

    shape =
    {
      type = Triangle3;
      intScheme = Gauss1;
    };
  };

  diri =
  {
    type = Dirichlet;

    groups = [ l , b , r , t  ];
    dofs   = [ dx, dy, dx, dy ];
    values = [ 0., 0., 1., 1. ];
  };
};

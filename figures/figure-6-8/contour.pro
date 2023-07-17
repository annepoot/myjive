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

      E = 3.;
      nu = 0.2;
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

    obsNoise = 1e-8;
    pdNoise = 1e-6;
    bcNoise = 1e-6;

    boundary =
    {
      type = dirichlet;
      groups = [ r , t ];
      dofs   = [ dx, dy ];
      covs   = [ 1.0, 1.0 ];
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

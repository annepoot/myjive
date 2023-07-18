gpinit =
{
  type = GPInit;

  mesh =
  {
    type = gmsh;
    file = meshes-mixed/plate_r2.msh;
  };

  coarseMesh =
  {
    type = gmsh;
    file = meshes-mixed/plate_r01.msh;
  };

  nodeGroups = [ l, lb, r ];

  l =
  {
    xtype = min;
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
};

gpsolver =
{
  type = GPSolver;
  nsample = 100;
  seed = 0;
  tables = [ strain ];
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

    obsNoise = 1e-2;
    pdNoise = 1e-6;
    bcNoise = 1e-4;

    boundary =
    {
      type = dirichlet;
      groups = [ l , r ];
      dofs   = [ dx, dx ];
      covs   = [ 100., 100. ];
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

    groups = [ l , lb , r ];
    dofs   = [ dx, dy, dx ];
    values = [ 0., 0., 1. ];
  };
};

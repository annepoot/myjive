gpinit =
{
  type = GPInit;

  mesh =
  {
    type = gmsh;
    file = meshes/plate_r1.msh;
  };

  coarseMesh =
  {
    type = gmsh;
    file = meshes/plate_r0.msh;
  };

  nodeGroups = [ l, lb ];

  l =
  {
    xtype = min;
  };

  lb =
  {
  	xtype = min;
    ytype = min;
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

  models = [ solid, load, gp, diri ];

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

  load =
  {
    type = Load;

    elements = all;

    dofs   = [ dx  ];
    values = [ 1.0 ];

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

    obsNoise = 1e-6;
    pdNoise = 1e-6;
    bcNoise = 1e-6;

    shape =
    {
      type = Triangle3;
      intScheme = Gauss1;
    };
  };

  diri =
  {
    type = Dirichlet;

    groups = [ l , lb ];
    dofs   = [ dx, dy ];
    values = [ 0., 0. ];
  };
};

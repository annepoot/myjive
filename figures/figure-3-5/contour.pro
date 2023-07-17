gpinit =
{
  type = GPInit;

  mesh =
  {
    type = gmsh;
    file = meshes/singularity-r3.msh;
  };

  coarseMesh =
  {
    type = gmsh;
    file = meshes/singularity-r1.msh;
  };

  nodeGroups = [ l, lb, tm ];

  l =
  {
    xtype = min;
  };

  lb =
  {
    xtype = max;
    ytype = min;
  };

  tm =
  {
    ytype = max;
    xtype = mid;
  };
};

gpsolver =
{
  type = GPSolver;
  priorMean = dirichlet;
};

model =
{
  type = Multi;

  models = [ solid, gp, diri, neum ];

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
      func = alpha**2 * M + beta**2 * F;
      hyperparams =
      {
        alpha = 1.;
        beta = 1.;
      };
    };

    obsNoise = 1e-7;
    pdNoise = 1e-5;
    bcNoise = 1e-8;

    shape =
    {
      type = Quad4;
      intScheme = Gauss4;
    };
  };

  diri =
  {
    type = Dirichlet;

    groups = [ l , l  ];
    dofs   = [ dx, dy ];
    values = [ 0., 0. ];
  };

  neum =
  {
    type = Neumann;

	groups = [ tm ];
	dofs   = [ dx ];
	values = [ 1. ];
  };
};

gpinit =
{
  type = GPInit;

  mesh =
  {
    type = gmsh;
    file = beam.msh;
  };

  coarseMesh =
  {
    type = gmsh;
    file = beam_coarse.msh;
  };

  nodeGroups = [ lb, rb, tm ];

  lb =
  {
    xtype = min;
  };

  rb =
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

  nsteps = 1;
  storeMatrix = True;
  storeConstraints = True;
  getUnitMassMatrix = True;
};

gpsampler =
{
  type = GPSampler;

  nsample = 50;
  seed = 5; // 3, 5
};

model =
{
  type = Multi;

  models = [ elastic, gp, diri ];

  elastic =
  {
    type = XElastic;

    elements = all;

    young = 10000.;
    thickness = 0.2;
    poisson = 0.2;
    rho = 1;
    state = plane_stress;

    shape =
    {
      type = Triangle3;
      intScheme = Gauss1;
    };
  };

  gp =
  {
     type = GP;

     obsNoise = 1e-5;
     alpha = 0.2;

     shape =
     {
       type = Triangle3;
       intScheme = Gauss1;
     };
  };

  diri =
  {
    type = Dirichlet;

    groups = [ lb, lb ];
    dofs   = [ dx, dy ];
    values = [ 0., 0. ];
  };

  neum =
  {
    type = Neumann;

    groups = [ tm ];
    dofs   = [ dy ];
    values = [ -1. ];
  };
};

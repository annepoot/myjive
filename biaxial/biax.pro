gpinit =
{
  type = GPInit;

  mesh =
  {
    type = gmsh;
    file = biax_fine.msh;
  };

  coarseMesh =
  {
    type = gmsh;
    file = biax_coarse.msh;
  };

  nodeGroups = [ top, left, bottom, right ];

  top =
  {
    ytype = max;
  };

  bottom =
  {
    ytype = min;
  };

  left =
  {
    xtype = min;
  };

  right =
  {
    xtype = max;
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

  nsample = 3;
  seed = 0;
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
    rho = 1000;
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
     alpha = opt;

     shape =
     {
       type = Triangle3;
       intScheme = Gauss1;
     };
  };

  diri =
  {
    type = Dirichlet;

    groups = [ top, bottom, left, right ];
    dofs   = [ dy, dy, dx, dx ];
    values = [ 0., 0., 0., 0. ];
  };
};

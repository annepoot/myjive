gpinit =
{
  type = GPInit;

  mesh =
  {
    type = gmsh;
    file = beam_fine2.msh;
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
    ytype = min;
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
  solver = CG;
  preconditioner = ichol;
  storeMatrix = True;
  storeConstraints = True;
  getUnitMassMatrix = True;
  getForceResults = True;
  tables = [ stress, strain ];
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

      E = 10000.;
      nu = 0.2;
      rho = 1.0;
      thickness = 0.2;
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
        alpha = opt;
      };
    };

    obsNoise = 1e-10;

    shape =
    {
      type = Triangle3;
      intScheme = Gauss1;
    };
  };

  diri =
  {
    type = Dirichlet;

    groups = [ lb, lb, rb ];
    dofs   = [ dx, dy, dy ];
    values = [ 0., 0., 0. ];
  };
};

vtkout =
{
  file = results;
  type = VTKOut;
  tables = [ stress, strain ];
};

gpinit =
{
  type = GPInit;

  mesh =
  {
    type = gmsh;
    file = tri6mesh.msh;
  };

  coarseMesh =
  {
    type = gmsh;
    file = mesh.msh;
  };

  nodeGroups = [ left, right, bottom ];

  left =
  {
    xtype = min;
  };

  right =
  {
    xtype = max;
  };

  bottom =
  {
    ytype = min;
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

  nsample = 10;
  seed = None;
};

model =
{
  type = Multi;

  models = [ poisson, gp, diri, neum ];

  poisson =
  {
    type = Poisson;

    elements = all;

    kappa = 1.0;
    rho = 1;

    shape =
    {
      type = Triangle6;
      intScheme = Gauss3;
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

    groups = [ bottom ];
    dofs   = [ u ];
    values = [ 0.0 ];
  };

  neum =
  {
    type = Neumann;

    groups = [ right ];
    dofs   = [ u ];
    values = [ 0.4 ];
  };
};

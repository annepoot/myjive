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
    file = tri3mesh.msh;
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

  nsample = 3;
  seed = 0;
};

model =
{
  type = Multi;

  models = [ poisson, gp, diri ];

  poisson =
  {
    type = XPoisson;

    elements = all;

    kappa = 1.0;
    rho = 1;
    q = (x<1)*1+(x>1)*-0.8;

    shape =
    {
      type = Triangle6;
      intScheme = Gauss3;
    };
  };

  gp =
  {
     type = GP;

     obsNoise = 1e-4;
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

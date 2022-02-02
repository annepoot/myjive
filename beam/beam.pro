init =
{
  nodeGroups = [ lb, rb, tm ];

  mesh =
  {
    type = gmsh;
    file = beam.msh;
  };

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

model =
{
  type = Multi;

  models = [ elastic, diri, neum ];

  elastic =
  {
    type = Elastic;

    elements = all;

    young = 10000.;
    thickness = 0.2;
    poisson = 0.2;
    state = plane_stress;

    shape =
    {
      type = Triangle3;
      intScheme = Gauss1;
    };
  };

  diri =
  {
    type = Dirichlet; 

    groups = [ lb, lb, rb, rb ];
    dofs   = [ dx, dy, dx, dy ];
    values = [ 0., 0., 0., 0. ];
  };
  
  neum =
  {
    type = Neumann; 

    groups = [ tm ];
    dofs   = [ dy ];
    values = [ -1. ];
  };
};

solver =
{
  nsteps = 1;
};

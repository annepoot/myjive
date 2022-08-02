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

  models = [ solid, diri, neum ];

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
//      rho = 1.0;
      thickness = 0.2;
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

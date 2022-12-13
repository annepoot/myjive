init =
{
  nodeGroups = [ lb, rb, tm ];

  mesh =
  {
    type = gmsh;
    file = ../beam.msh;
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

  models = [ solid, diri ];

  solid =
  {
    type = Solid;

    elements = all;

    material =
    {
      type = Deteriorated;
      rank = 2;
      anmodel = plane_stress;

      E = 10000.;
      nu = 0.2;
      rho = 1.0;
      thickness = 0.2;

      deteriorations = 5;
      scale = 0.1;
      seed = 0;
      locX = x;
      locY = y;
      stdX = np.random.uniform(0.5, 1.0);
      stdY = np.random.uniform(0.5, 1.0);
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

    groups = [ lb, lb, rb ];
    dofs   = [ dx, dy, dy ];
    values = [ 0., 0., 0. ];
  };
};

solver =
{
  type = Linsolve;
  tables = [stiffness];
};

vtkout =
{
  type = VTKOut;
  tables = [stiffness];
  file = stiffness;
};

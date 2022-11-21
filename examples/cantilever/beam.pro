init =
{
  nodeGroups = [ lb, rb, bottom, top, left, right, rt ];

  mesh =
  {
    type = gmsh;
    file = beam.msh;
  };

  bottom =
  {
    ytype = min;
  };

  top =
  {
    ytype = max;
  };

  left =
  {
    xtype = min;
  };

  right =
  {
    xtype = max;
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

  rt =
  {
    ytype = max;
    xtype = max;
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

    young = 1000.;
    poisson = 0.3;
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

    groups = [ left, left ];
    dofs   = [ dx, dy ];
    values = [ 0., 0. ];
  };

  neum =
  {
    type = Neumann;

    groups = [ rt ];
    dofs = [ dy ];
    values = [ -1.0 ];
  };
};

solver =
{
  type = Linsolve;
  nsteps = 1;
};

vtkout =
{
  type = VTKOut;
  file = results;
  tables = [ stress ];
};

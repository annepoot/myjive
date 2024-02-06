init =
{
  type = Init;

  mesh =
  {
    type = gmsh;
    file = 2partbeam.msh;
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

model =
{
  type = Multi;

  models = [ left, right, diri, neum ];

  left =
  {
    type = Solid;

    elements = left;

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
      type = Triangle3;
      intScheme = Gauss1;
    };
  };

  right =
  {
    type = Solid;

    elements = right;

    material =
    {
      type = Isotropic;
      rank = 2;
      anmodel = plane_stress;

      E = 100000000.;
      nu = 0.2;
    };

    thickness = 0.2;

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
  type = Linsolve;
};

init =
{
  nodeGroups = [ lb, rb ];

  mesh =
  {
    type = gmsh;
    file = elem4.msh;
  };

  lb =
  {
    xtype = min;
  };

  rb =
  {
    xtype = max;
  };
};

model =
{
  type = Multi;

  models = [ solid, load, diri, neum ];

  solid =
  {
    type = Solid;

    elements = all;

    material =
    {
      type = Isotropic;
      rank = 2;
      anmodel = plane_stress;

      E = 1.;
      nu = 0.;
    };

    thickness = 1.;

    shape =
    {
      type = Quad4;
      intScheme = Gauss4;
    };
  };

  load =
  {
    type = Load;

    elements = all;

    dofs   = [ dy ];
    values = [ 0. ];

    shape =
    {
      type = Quad4;
      intScheme = Gauss4;
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

    groups = [ rb ];
    dofs   = [ dx ];
    values = [ 1. ];
  };
};

solver =
{
  type = Linsolve;
};

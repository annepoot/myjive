modules = [ init, solver ];

init =
{
  nodeGroups = [ lb, rb ];

  mesh =
  {
    type = gmsh;
    file = elem9.msh;
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
  models = [ solid, load, diri ];

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
      type = Quad9;
      intScheme = Gauss9;
    };
  };

  load =
  {
    type = Load;

    elements = all;

    dofs   = [ dx ];
    values = [ 1. ];

    shape =
    {
      type = Quad9;
      intScheme = Gauss9;
    };
  };

  diri =
  {
    type = Dirichlet;

    groups = [ lb, lb ];
    dofs   = [ dx, dy ];
    values = [ 0., 0. ];
  };
};

solver =
{
  type = Linsolve;
};

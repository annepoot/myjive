init =
{
  type = Init;

  mesh =
  {
    type = gmsh;
    file = meshes/beam_coarse.msh;
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

solver =
{
  type = Linsolve;
  solver = CG;
};

model =
{
  type = Multi;

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

  load =
  {
    type = Load;

    elements = all;

    dofs   = [ dy ];
    values = [ -0.2 ];

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

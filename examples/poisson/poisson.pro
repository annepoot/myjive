init =
{
  nodeGroups = [ left, right ];

  mesh =
  {
    type = gmsh;
    file = mesh.msh;
  };

  left =
  {
    xtype = min;
  };

  right =
  {
    xtype = max;
  };
};

model =
{
  type = Multi;

  models = [ bar, diri ];

  bar =
  {
    type = Poisson;

    elements = all;

    kappa = 10.0;

    shape =
    {
      type = Triangle3;
      intScheme = Gauss1;
    };
  };

  diri =
  {
    type = Dirichlet;

    groups = [ left, right ];
    dofs   = [ u, u ];
    values = [ 0.0, 1.0 ];
  };
};

solver =
{
  type = Linsolve;
  solver = CG;
  preconditioner = ichol;
};

init =
{
  nodeGroups = [ left, right ];

  mesh =
  {
    type = geo;
    file = bar.geom;
  };

  left = 0;
  right = [1];
};

model =
{
  type = Multi;

  models = [ bar, diri, neum ];

  bar =
  {
    type = Bar;

    elements = all;

    EA = 1.0;

    k = 1.0;

    shape =
    {
      type = Line2;
      intScheme = Gauss1;
    };
  };

  diri =
  {
    type = Dirichlet;

    groups = [ right ];
    dofs   = [ dx ];
    values = [ 0.0 ];
  };

  neum =
  {
    type = Neumann;

    groups = [ left ];
    dofs = [ dx ];
    values = [ -1.0 ];
  };
};

solver =
{
  type = Linsolve;
};

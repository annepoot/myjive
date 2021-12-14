init =
{
  mesh = bar.msh;

  nodeGroups = [ left, right ];

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

  models = [ bar, diri, neum ];

  bar =
  {
    type = Bar;

    elements = all;

    young = 1.0;

    shape =
    {
      type = Line2;
      intScheme = Gauss1;
    };
  };

  diri =
  {
    type = Dirichlet; 

    groups = [ left ];
    dofs   = [ dx ];
    values = [ 0.0 ];
  };

  neum = 
  {
    type = Neumann;

    groups = [ right ];
    dofs = [ dx ];
    values = [ 1.0 ];
  };
};

solver =
{
  nsteps = 1;
};

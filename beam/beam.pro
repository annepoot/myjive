init =
{
  mesh = beam.msh;

  nodeGroups = [ lb, rb, bottom, top, left, right, tm ];

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

  tm =
  {
    ytype = max;
    xtype = mid;
  };
};

model =
{
  type = Multi;

  models = [ elastic, diri ];

  elastic =
  {
    type = Elastic;

    elements = all;

    young = 1.;
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

    groups = [ lb, lb, rb, tm ];
    dofs   = [ dx, dy, dy, dy ];
    values = [ 0., 0., 0., -1. ];
  };
};

solver =
{
  nsteps = 1;
};

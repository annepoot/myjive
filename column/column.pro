init =
{
  nodeGroups = [ bot, top ];

  mesh = 
  {
    type = geo;
    file = column.geom;
  };

  bot = [0];
  top = [1];
};

model =
{
  type = Multi;

  models = [ frame, diri, neum ];

  frame =
  {
    type = Frame;

    elements = all;

    subtype = linear;

    EA = 1.e6;

    GAs = 20.e3;

    EI = 3.e3;

    shape =
    {
      type = Line2;
      intScheme = Gauss1;
    };
  };

  diri =
  {
    type = Dirichlet; 

    groups = [ bot, bot, top ];
    dofs   = [ dx, dy, dx ];
    values = [ 0.0, 0.0, 0.0 ];
  };

  neum = 
  {
    type = Neumann;

    groups = [ top ];
    dofs = [ dy ];
    values = [ -1. ];
  };
};

solver =
{
  nsteps = 1;
  storeMatrix = False;
};

linbuck = 
{
};

frameview = 
{
  deform = 10.;
};

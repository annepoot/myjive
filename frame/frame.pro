init =
{
  nodeGroups = [ bl, tl, tr, br ];

  mesh = 
  {
    type = geo;
    file = frame.geom;
  };

  bl = [0];
  tl = [1];
  tr = [2];
  br = [3];
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

    EA = 20000;

    GAs = 10000;

    EI = 10;

    shape =
    {
      type = Line2;
      intScheme = Gauss1;
    };
  };

  diri =
  {
    type = Dirichlet; 

    groups = [ bl, bl, br, br ];
    dofs   = [ dx, dy, dx, dy ];
    values = [ 0.0, 0.0, 0.0, 0.0 ];
  };

  neum = 
  {
    type = Neumann;

    groups = [ tl, tl, tr ];
    dofs = [ dx, dy, dy ];
    values = [ 1000, -5000, -5000 ];
  };
};

solver =
{
  nsteps = 1;
  storeMatrix = False;
};

init =
{
  nodeGroups = [ left, right, mid ];

  mesh = 
  {
    type = geo;
    file = euler.geom;
  };

  left = [0];
  mid = [1];
  right = [2];
};

model =
{
  type = Multi;

  models = [ frame, diri, neum ];

  frame =
  {
    type = Frame;

    elements = all;

    subtype = nonlin;

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

    groups = [ left, left, right ];
    dofs   = [ dy, dx, dy ];
    values = [ 0.0, 0.0, 0.0 ];
  };

  neum = 
  {
    type = Neumann;

    groups = [ right, mid ];
    dofs = [ dx, dy ];
    values = [ -.1, 0.001 ];
    loadIncr = [ -.1, 0. ];
  };
};

nonlin =
{
  nsteps = 10;
  itermax = 10;
  tolerance = 1e-6;
};

loaddisp = 
{
  groups = [ left, mid, right ];
};

frameview =
{
  deform = 10.;
  interactive = True;
};

graph = 
{
  xData = [loaddisp.right.disp.dx,loaddisp.mid.disp.dy];
  yData = [loaddisp.right.load.dx,loaddisp.right.load.dx];
};

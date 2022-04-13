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

    groups = [ bl, bl, br, br, tr ];
    dofs   = [ dx, dy, dx, dy, dy ];
    values = [ 0.0, 0.0, 0.0, 0.0, 0.0 ];
    dispIncr = [ 0.0, 0.0, 0.0, 0.0, -0.0002 ];
  };

  neum = 
  {
    type = Neumann;

    groups = [ tl ];
    dofs = [ dx ];
    values = [ -0.2 ];
    loadIncr = [ 0. ];
  };
};

nonlin =
{
  nsteps = 10;
  itermax = 10;
  tolerance = 1e-6;
};

frameview =
{
  deform = 10.;
};

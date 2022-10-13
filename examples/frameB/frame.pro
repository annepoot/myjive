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
    plastic = True;

    EA = 20000;
    GAs = 10000;
    EI = 10;
    Mp = 0.4;

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
    dispIncr = [ 0.0, 0.0, 0.0, 0.0 ];
  };

  neum =
  {
    type = Neumann;

    groups = [ tl, tl, tr ];
    dofs = [ dx, dy, dy ];
    values = [ 0, 0, 0 ];
    loadIncr = [ 1000, -5000, -5000 ];
  };
};

arclen =
{
  nsteps = 100;
  itermax = 10;
  tolerance = 1e-6;
  beta = 0.1;
  dl = 0.02;
};

loaddisp =
{
  type = LoadDisp;
  groups = [ tl ];
};

frameview =
{
  type = FrameView;
  deform = 1.;
};

graph =
{
  xData = [loaddisp.tl.disp.dx,loaddisp.tl.disp.dy];
  yData = [loaddisp.tl.load.dx,loaddisp.tl.load.dy];
};

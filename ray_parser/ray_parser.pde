///////////////////////////////////////////////////////////////////////
  //
  //Command Line Interface (CLI) Parser  
  //
  ///////////////////////////////////////////////////////////////////////
  String gCurrentFile = new String("data/t0.cli"); // A global variable for holding current active file name.
  
  ///////////////////////////////////////////////////////////////////////
  //
  //Press key 1 to 9 and 0 to run different test cases.
  //
  ///////////////////////////////////////////////////////////////////////
  public void keyPressed() {
  switch(key) {
  case '1':  gCurrentFile = new String("data/t0.cli"); interpreter(); break;
  case '2':  gCurrentFile = new String("data/t1.cli"); interpreter(); break;
  case '3':  gCurrentFile = new String("data/t2.cli"); interpreter(); break;
  case '4':  gCurrentFile = new String("data/t3.cli"); interpreter(); break;
  case '5':  gCurrentFile = new String("data/c0.cli"); interpreter(); break;
  case '6':  gCurrentFile = new String("data/c1.cli"); interpreter(); break;
  case '7':  gCurrentFile = new String("data/c2.cli"); interpreter(); break;
  case '8':  gCurrentFile = new String("data/c3.cli"); interpreter(); break;
  case '9':  gCurrentFile = new String("data/c4.cli"); interpreter(); break;
  case '0':  gCurrentFile = new String("data/c5.cli"); interpreter(); break;
  }
  }
  
  ///////////////////////////////////////////////////////////////////////
  //
  //Get a float value from a string.
  //
  ///////////////////////////////////////////////////////////////////////
  float get_float(String str) { return Float.parseFloat(str); }
  
  ///////////////////////////////////////////////////////////////////////
  //
  //Parser core. It parses the CLI file and processes it based on each 
  //token. Only "color", "rect", and "write" tokens are implemented. 
  //You should start from here and add more functionalities for your
  //ray tracer.
  //
  //Note: Function "splitToken()" is only available in processing 1.25 
  //or higher.
  //
  ///////////////////////////////////////////////////////////////////////
  void interpreter() {
  
  reset(); // clear all shapes
  String str[] = loadStrings(gCurrentFile);
  if (str == null) println("Error! Failed to read the file.");
  for (int i=0; i<str.length; i++) {
  
  String[] token = splitTokens(str[i], " "); // Get a line and parse tokens.
  if (token.length == 0) continue; // Skip blank line.
  
  if (token[0].equals("fov")) {
    fov = get_float(token[1]);
  }
  else if (token[0].equals("background")) {
    bgr = get_float(token[1]);
    bgg = get_float(token[2]);
    bgb = get_float(token[3]);
  }
  else if (token[0].equals("light")) {
    lights.add(new Light(get_float(token[1]), get_float(token[2]), get_float(token[3]), get_float(token[4]), get_float(token[5]), get_float(token[6])));
  }
  else if (token[0].equals("surface")) {
    CDr = get_float(token[1]);
    CDg = get_float(token[2]);
    CDb = get_float(token[3]);
    CAr = get_float(token[4]);
    CAg = get_float(token[5]);
    CAb = get_float(token[6]);
    CSr = get_float(token[7]);
    CSg = get_float(token[8]);
    CSb = get_float(token[9]);
    P = get_float(token[10]);
    KREFL = get_float(token[11]);
  }    
  else if (token[0].equals("begin")) {
    A = B = C = null;
  }
  else if (token[0].equals("vertex")){
    if (A == null){
      A = new Vertex(get_float(token[1]), get_float(token[2]), get_float(token[3]));
    } else if (B == null) {
      B = new Vertex(get_float(token[1]), get_float(token[2]), get_float(token[3]));
    } else {
      C = new Vertex(get_float(token[1]), get_float(token[2]), get_float(token[3]));
      shapes.add(new Polygon());
      A = B = C = null;
    }
  }
  else if (token[0].equals("end")) {
    A = B = C = null;
  }
  else if (token[0].equals("sphere")) {
    Sphere newSphere = new Sphere(get_float(token[1]), get_float(token[2]), get_float(token[3]), get_float(token[4]));
    shapes.add(newSphere);
  }
  else if (token[0].equals("color")) {
  float r =get_float(token[1]);
  float g =get_float(token[2]);
  float b =get_float(token[3]);
  fill(r, g, b);
  }
  else if (token[0].equals("rect")) {
  float x0 = get_float(token[1]);
  float y0 = get_float(token[2]);
  float x1 = get_float(token[3]);
  float y1 = get_float(token[4]);
  rect(x0, height-y1, x1-x0, y1-y0);
  }
  else if (token[0].equals("write")) {
  if ( !gCurrentFile.equals("D:\\Dropbox\\EclipseWorkspace\\RayTracer\\src\\data\\rect_test.cli")) {
    render();
  }
  save(token[1]);
  }
  }
  }
  
  
  //RENDERING ETC
  void render(){
    
    loadPixels();
    
    double k = Math.tan(Math.toRadians(fov/2.0));
    
    for (int r = 0; r < 300; r++){
      for (int c = 0; c < 300; c++){
        int pix = ((299-c)*300) + r;
        Ray ray = new Ray(new Vertex(0,0,0), new Vertex( (r-150) * (k/150.0), (c-150) * (k/150.0), -1));
        
        Intersectable closest = null;
        Vertex point = null;
        
        for (Intersectable shape : shapes){
          Vertex tempPoint = shape.intersect(ray);
          if (tempPoint != null && (point == null || tempPoint.distance(ray.origin) < point.distance(ray.origin))){
            closest = shape;
            point = tempPoint;
          }
        }
        
        if (closest == null){
          pixels[pix] = color(bgr, bgg, bgb);
        } else {
          float red = 0;
          float green = 0;
          float blue = 0;
          float[] colors = determineColor(ray, closest, point, red, green, blue, 0);
          pixels[pix] = color(colors[0], colors[1], colors[2]);
        }
        
      }
    }
    
    updatePixels();
  }
  
  float[] determineColor(Ray ray, Intersectable surface, Vertex point, float red, float green, float blue, int depth){
    Vertex normal = surface.getNormal(point);
    if (depth < 12 && surface.krefl > 0) {
      Vertex unitraydirec = ray.slope.unit();
      Ray reflDirec = new Ray(point, normal.scale(2 * normal.dot(unitraydirec)).minus(unitraydirec).scale(-1));
      Intersectable reflClosest = null;
      Vertex reflPoint = null;
      for (Intersectable shape : shapes){
        if (shape != surface){
          Vertex tempPoint = shape.intersect(reflDirec);
          if (tempPoint != null && ( reflPoint == null || tempPoint.distance(reflDirec.origin) < reflPoint.distance(reflDirec.origin))){
            reflClosest = shape;
            reflPoint = tempPoint;
          }
        }
      }
      if (reflClosest != null){
        float[] reflcolors = determineColor(reflDirec, reflClosest, reflPoint, red, green, blue, depth + 1);
        red += surface.krefl * reflcolors[0];
        green += surface.krefl * reflcolors[1];
        blue += surface.krefl * reflcolors[2];
      } else {
        red += surface.krefl * bgr;
        green += surface.krefl * bgg;
        blue += surface.krefl * bgb;
      }
    }
    
    
    
    red += surface.car;
    green += surface.cag;
    blue += surface.cab;
    
    for (Light l : lights){
      Ray lightray = new Ray(point,l.loc.minus(point));
      Vertex lPoint = null;
      for (Intersectable shape : shapes){
        Vertex tempPoint = shape.intersect(lightray);
        if (tempPoint != null && (lPoint == null || tempPoint.distance(lightray.origin) < point.distance(ray.origin))) {
          lPoint = tempPoint;
        }
      }
      
      
      
      if (lPoint == null || lPoint.distance(lightray.origin) > l.loc.distance(lightray.origin)){ //IF VISIBLE
        
        double costheta = Math.abs((lightray.slope.dot(normal) / (lightray.slope.length() * normal.length())));
        red += costheta * surface.cdr * l.r;
        green += costheta * surface.cdg * l.g;
        blue += costheta * surface.cdb * l.b;
        
        Vertex L = lightray.slope.unit();
        Vertex R = normal.scale(2 * normal.dot(L)).minus(L).unit();
        Vertex E = ray.slope.unit();
        if (R.dot(E) < 0){
          double phong = Math.pow(R.dot(E), surface.p);
          red += phong * surface.csr * l.r;
          green += phong * surface.csg * l.g;
          blue += phong * surface.csb * l.b;
        }
      }
    }
    float[] result = { red, green, blue };
    return result;
  }
  
  
  //NEEDED VARIABLES
  float P, KREFL;
  float CDr, CDg, CDb, CAr, CAg, CAb, CSr, CSg, CSb;
  float bgr, bgg, bgb;
  float fov;
  ArrayList<Light> lights;
  ArrayList<Intersectable> shapes;
  Vertex A, B, C;
  
  //EXTRA CLASSES
  abstract class Intersectable {
    public abstract Vertex intersect(Ray r);
    float cdr, cdg, cdb, car, cag, cab, csr, csg, csb, p, krefl;
    abstract Vertex getNormal(Vertex v);
  }
  
  class Ray {
    Vertex origin, slope;
    
    Ray(Vertex origin, Vertex slope){
      this.origin = origin;
      this.slope = slope;
    }
    
    Vertex getPoint(double t){
      return new Vertex(origin.x + t*slope.x, origin.y + t*slope.y, origin.z + t*slope.z);
    }
  }
  
  class Light {
    Vertex loc;
    float r, g, b;
    Light(float x, float y, float z, float r, float g, float b){
      loc = new Vertex(x, y, z);
      this.r = r; this.g = g; this.b = b;
    }
  }
  
  class Sphere extends Intersectable{
    float r;
    Vertex loc;
    Sphere(float r, float x, float y, float z){
      loc = new Vertex(x, y, z);
      this.r = r;
      cdr = CDr; cdg = CDg; cdb = CDb; car = CAr; cag = CAg; cab = CAb; csr = CSr; csg = CSg; csb = CSb;
      p = P; krefl = KREFL;
    }
    
    public Vertex intersect(Ray ray){
      
      double A = ray.slope.dot(ray.slope);
      double B = 2 * ((ray.origin.minus(loc)).dot(ray.slope));
      double C = (ray.origin.minus(loc)).dot(ray.origin.minus(loc)) - r*r;
      double disc = (B*B) - (4*A*C);
      
      
      if (disc < 0){
        return null;
      }
      double root1 = (-B + Math.sqrt(disc)) / (2*A);
      double root2 = (-B - Math.sqrt(disc)) / (2*A);
      double root = Math.min(root1, root2);
      if (root > -.000001){
        return ray.getPoint(root);
      } else {
        return null;
      }
    }

    Vertex getNormal(Vertex v) {
      return v.minus(loc).scale(-1).unit();
    }
  }
  
  class Vertex {
    double x, y, z;
    
    Vertex(double x, double y, double z){
      this.x = x; this.y = y; this.z = z;
    }
    
    double dot(Vertex v){
      return (x*v.x + y*v.y + z*v.z);
    }
    
    Vertex minus(Vertex v){
      return new Vertex(x-v.x, y-v.y, z-v.z);
    }
    
    double distance(Vertex a){
      return Math.sqrt( (a.x - x)*(a.x - x) + (a.y - y)*(a.y - y) + (a.z - z)*(a.z - z) );
    }
    
    Vertex crossProduct(Vertex v){
      return new Vertex(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
    }
    
    Vertex unit(){
      double length = Math.sqrt(dot(this));
      return new Vertex(x / length, y / length, z / length);
    }
    
    double length() {
      return Math.sqrt(dot(this));
    }
    
    Vertex scale(double s){
      return new Vertex(s*x, s*y, s*z);
    }
  }
  
  class Polygon extends Intersectable{
    Vertex a, b, c, normal;
    double d;
    Polygon(){
      a = A; b = B; c = C;
      cdr = CDr; cdg = CDg; cdb = CDb; car = CAr; cag = CAg; cab = CAb; csr = CSr; csg = CSg; csb = CSb;
      p = P; krefl = KREFL;
      normal = b.minus(a).crossProduct(c.minus(a)).unit();
      d = normal.dot(a);
    }
    
    public Vertex intersect(Ray r) {
      double denom = normal.dot(r.slope);
      if (denom == 0) {
        return null;
      }
      double t = (d - normal.dot(r.origin)) / denom;
      if (t <= .000001) {
        return null;
      }
      Vertex point = r.getPoint(t);
      if ( b.minus(a).crossProduct(point.minus(a)).dot(normal) >= -.0000000000001 &&
          c.minus(b).crossProduct(point.minus(b)).dot(normal) >= -.0000000000001 &&
          a.minus(c).crossProduct(point.minus(c)).dot(normal) >= -.0000000000001) {
        //JANKY FLOATING POINT MATH GAWD
        return point;
      }
      return null;
    }

    @Override
    Vertex getNormal(Vertex v) {
      return normal;
    }
  }
  
  void reset(){
    lights = new ArrayList<Light>();
    shapes = new ArrayList<Intersectable>();
    A = B = C = null;
    bgr = bgg = bgb = 0;
    int bg = color(bgr, bgg, bgb);
    loadPixels();
    for (int i = 0; i < pixels.length; i++){
      pixels[i] = bg;
    }
    updatePixels();
  }
  
  ///////////////////////////////////////////////////////////////////////
  //
  //Some initializations for the scene.
  //
  ///////////////////////////////////////////////////////////////////////
  public void setup() {
  size(300, 300);
  noStroke();
  colorMode(RGB, 1.0f);
  background(0, 0, 0);
  interpreter();
//  fill(255,255,255);
//  rect(100,100,100,100);
  }
  
  ///////////////////////////////////////////////////////////////////////
  //
  //Draw frames.
  //
  ///////////////////////////////////////////////////////////////////////
  public void draw() {
  
  }

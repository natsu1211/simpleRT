/*------------------------------------------------------------------------------
 * a simple C# version raytracer for study purpose
 * the setting of cornell box is from the famous 99 lines pathtracer smallpt
 * 2014 Created by Long
-------------------------------------------------------------------------------*/

using System;
using System.IO;

namespace simpleRT
{
    using Color = simpleRT.Vec3;
    class simpleRT
    {
        const double PI = 3.14159265358979323846;
        const double INF = 1e20;
        const double EPS = 1e-6;
        const double MaxDepth = 5;
        public struct Vec3
        {
            public double x,y,z;
            public Vec3(double x_ = 0.0, double y_ = 0.0, double z_ = 0.0) { x = x_; y = y_; z = z_; }
            public static Vec3 operator +(Vec3 lhs, Vec3 rhs) { return new Vec3(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z); }
            public static Vec3 operator -(Vec3 lhs, Vec3 rhs) { return new Vec3(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z); }
            public static Vec3 operator *(Vec3 lhs, double n) { return new Vec3(lhs.x * n, lhs.y * n, lhs.z * n); }
            public static Vec3 operator *(double n, Vec3 rhs) { return new Vec3(rhs.x * n, rhs.y * n, rhs.z * n); }
            public static Vec3 operator /(Vec3 lhs, double n) { return new Vec3(lhs.x / n, lhs.y / n, lhs.z / n); }
            public double Length() { return Math.Sqrt(x * x + y * y + z * z); }
            //public double SquaredLength() { return x * x + y * y + z * z; }
            public double Dot(Vec3 v) { return x * v.x + y * v.y + z * v.z; }
            public Vec3 Cross(Vec3 v) { return new Vec3((y * v.z) - (z * v.y), (z * v.x) - (x * v.z), (x * v.y) - (y * v.x)); }
            /// <summary>
            /// for color blend
            /// </summary>
            public Vec3 Multi(Vec3 v) { return new Vec3(x * v.x,  y * v.y, z * v.z); }
            public Vec3 Norm() { return this = this * (1 / Math.Sqrt(x * x + y * y + z * z)); }
        }
        
        public struct Ray
        {
            public Vec3 org, dir;
            public Ray(Vec3  org_, Vec3 dir_) { org = org_; dir = dir_; }
        }

        /// <summary>
        /// The Reflection type of material, using in Radiance()
        /// </summary>
        public enum ReflectionType
        {
            DIFFUSE,
            SPECULAR,
            REFRACTION,
        };
        public class Sphere
        {
            public double rad;
            public Vec3 pos;
            public ReflectionType reflType;
            public Vec3 emission, color;
            public Sphere(double rad_, Vec3 pos_, Vec3 emission_, Vec3 color_, ReflectionType refltype_) { rad = rad_; pos = pos_; reflType = refltype_; emission = emission_; color = color_; }
            public double Intersect(Ray ray)
            {
                Vec3 op = pos - ray.org;
                var b= ray.dir.Dot(op);
                var d4 = b * b - op.Dot(op) + rad * rad;
                /// <summary>
                /// see if the ray is intersectd with sphere
                /// </summary>

                if (d4 >= 0.0)
                {
                    double sqrtD4 = Math.Sqrt(d4);
                    double s1 = b - sqrtD4;
                    double s2 = b + sqrtD4;
                    if (s1 > EPS)
                        return s1;
                    else if (s2 > EPS)
                        return s2;
                }
                //else no solution
                 return 0.0;
            }
        }

        /// <summary>
        /// setting of cornell box
        /// </summary>
        public static Sphere[] Spheres =
    {
        new Sphere(1e5, new Vec3( 1e5+1,40.8,81.6), new Color(), new Color(0.75, 0.25, 0.25), ReflectionType.DIFFUSE),// left wall
        new Sphere(1e5, new Vec3(-1e5+99,40.8,81.6),new Color(), new Color(0.25, 0.25, 0.75), ReflectionType.DIFFUSE),// right wall
        new Sphere(1e5, new Vec3(50,40.8, 1e5), new Color(), new Color(0.75, 0.75, 0.75), ReflectionType.DIFFUSE),// back wall
        new Sphere(1e5, new Vec3(50,40.8,-1e5+170), new Color(), new Color(), ReflectionType.DIFFUSE),// front wall
        new Sphere(1e5, new Vec3(50, 1e5, 81.6), new Color(), new Color(0.75, 0.75, 0.75), ReflectionType.DIFFUSE),// bed
        new Sphere(1e5, new Vec3(50,-1e5+81.6,81.6), new Color(), new Color(0.75, 0.75, 0.75), ReflectionType.DIFFUSE),// ceiling
        new Sphere(16.5,new Vec3(65,20,20), new Color(), new Color(0.25,0.75,0.25), ReflectionType.DIFFUSE),// ball
        new Sphere(16.5,new Vec3(27,16.5,47), new Color(), new Color(1,1,1)*.99, ReflectionType.SPECULAR),// mirror 
        new Sphere(16.5,new Vec3(73,16.5,78), new Color(), new Color(1,1,1)*.99, ReflectionType.REFRACTION),//glass
    };
        
        static bool IntersectScene(Ray ray, ref double distance, ref int id)
        {
            double d;
            distance = INF;
            for (int i = Spheres.Length - 1; i >= 0; i--)
                if ((d = Spheres[i].Intersect(ray)) != 0 && d < distance) 
                { 
                    distance = d; 
                    id = i; 
                }
            return distance < INF;
        }

        /// <summary>
        /// point light source
        /// </summary>
        public static Vec3 lightPos = new Vec3(50.0, 75.0, 81.6);
        public static Color lightColor = new Color(256, 256, 256);
        public static Color Radiance(Ray ray, int depth)
        {
            int id = 0;
            double distance = 0;
            if (! IntersectScene(ray, ref distance, ref id))//no intersect
                return new Color(0.0,0.0,0.0);//return background color(black)
            var intersectedObj = Spheres[id];
            Vec3 hitPoint = ray.org + ray.dir * distance;
            Vec3 normal = (hitPoint - intersectedObj.pos).Norm(); 
            Vec3 orientedNormal = normal.Dot(ray.dir) < 0.0 ? normal : (normal * -1.0);

            //terimation of recursion
            if (depth > MaxDepth)
                return new Color();

            switch (intersectedObj.reflType)
            {
                //no recursive raytrace by diffuse reflection here, since it's too time consuming
                case ReflectionType.DIFFUSE:
                    {
                        double distance1 = 0;
                        int id1 = 0;
                        Vec3 lightDir = lightPos - hitPoint;
                        double length = lightDir.Length();
                        IntersectScene(new Ray(hitPoint, lightDir/length), ref distance1, ref id1);
                        //visibility judge
                        if(distance1 >= length)//visible
                        {
                            return intersectedObj.emission + Math.Max(0.0, orientedNormal.Dot(lightDir / length)) * lightColor.Multi(intersectedObj.color) / (length * length);
                        }
                        else//invisible
                        {
                            return new Color();
                        }
                        //break;
                    }
                case ReflectionType.SPECULAR:
                    {
                        var a = intersectedObj.emission + intersectedObj.color.Multi( Radiance(new Ray(hitPoint, ray.dir - 2.0 * normal * normal.Dot(ray.dir)), depth+1));
                        return a;
                    }
                case ReflectionType.REFRACTION:
                    {
                        var reflRay = new Ray(hitPoint, ray.dir -  2.0 * normal * normal.Dot(ray.dir));
                        //see whether the ray is from inside(out object) or outside(into object)
                        bool intoObject = normal.Dot(orientedNormal) > 0.0;
                        // refractive index of vacuum
                        const double n1 = 1.0;
                        // refractive index of glass ball
                        const double n2 = 1.5; 
                        // Snell's Law
                        double nn = intoObject ? (n1 / n2) : (n2 / n1);
                        double dn = ray.dir.Dot(orientedNormal);
                        double cos2Theta2 = 1.0 - nn * nn * (1.0 - dn * dn);
                        //theta2 > 90 degree, totle reflection(no refraction)
                        if (cos2Theta2 < 0.0)
                        {
                            return intersectedObj.emission + Radiance(reflRay, depth + 1);
                        }
                        // else refraction
                        Vec3 refrDir = (ray.dir * nn - normal * (intoObject ? 1.0 : -1.0) * (dn * nn + Math.Sqrt(cos2Theta2))).Norm();

                        // Schlick approximation of Fresnel Equation
                        double a = n2 - n1, b = n2 + n1;
                        double f0 = (a * a) / (b * b);
                        double c = 1.0 - (intoObject ? -dn : refrDir.Dot(normal));
                        double fRefl = f0 + (1.0 - f0) * Math.Pow(c, 5.0);
                        double fRefr = 1.0 - fRefl;

                        return intersectedObj.emission +
                            intersectedObj.color.Multi(Radiance(reflRay, depth + 1) * fRefl + Radiance(new Ray(hitPoint, refrDir), depth + 1) * fRefr);
                    }
                default:
                    return new Color();
                    
            }

        }
        static double Clamp(double x) { return x < 0 ? 0 : (x > 1 ? 1 : x); }
        static int ToInt(double x) { return (int)(Math.Pow(Clamp(x), 1 / 2.2) * 255 + .5); }
        static void Main(string[] args)
        {
            DateTime startTime = DateTime.Now;
            const int width = 640;
	        const int height = 480;
	        //camera setting
	        var camera = new Ray(new Vec3(50.0, 52.0, 295.6), new Vec3(0.0, -0.042612, -1.0).Norm());
	        //x,y vector of screen
	        Vec3 cx = new Vec3(width * .5135 / height);
            Vec3 cy = cx.Cross(camera.dir).Norm() * .5135 ; 
            //image buffer
            Vec3[] image = new Vec3[width * height];
            for (int y = 0; y < height;++y)
            {
                Console.Write("\rRendering {0:F2}%", 100.0 * y / (height - 1));
                for (int x = 0 ; x < width; ++x)
                {
                    int index = (height - y - 1) * width + x;
                    // 2x2 subpixel sampling
                    for (int sy = 0; sy < 2; ++sy)
                    {
                        for (int sx = 0; sx < 2; ++sx)
                        {
                            double dx = sx / 2.0;
                            double dy = sy / 2.0;
                            Vec3 dir = cx * (((sx + 0.5 + dx) / 2.0 + x) / width - 0.5) +
                                        cy * (((sy + 0.5 + dy) / 2.0 + y) / height - 0.5) + camera.dir;
                            image[index] = image[index] + Radiance(new Ray(camera.org + dir * 130.0, dir.Norm()), 0);
                            
                        }
                        
                    }
                }
            }
                Console.WriteLine("\n{0} sec", (DateTime.Now - startTime).TotalSeconds);
                using (StreamWriter sw = new StreamWriter("image.ppm"))
                {
                    sw.Write("P3\r\n{0} {1}\r\n{2}\r\n", width, height, 255);
                    for (int i = 0; i < width * height; i++)
                        sw.Write("{0} {1} {2}\r\n", ToInt(image[i].x), ToInt(image[i].y), ToInt(image[i].z));
                    sw.Close();
                }
        }
    }
}

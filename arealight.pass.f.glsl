precision highp float;
precision highp int;

in vec3 vpos;
in vec3 vnormal;

in vec3 vcolor;
in vec3 wcoord;
in vec3 wvnormal;
in vec3 vdepth;

uniform vec3 eone;
uniform float rone;
uniform float exposure;
uniform vec4 pointlight;

out vec4 color;

void main() {
    vec3 bright = vec3(0,0,0);
    float pi = 3.1415926535897932384626433832795;
    float incdegree = pi / 12;
    vec3 z = vec3(0,0,1);

    //center point of area light based on pointlight vector
    vec3 l = vec3(pointlight[0] - vpos[0], pointlight[1] - vpos[1], pointlight[2] - vpos[2]);



    //vcolor copied from pass.f.glsl
    vec3 colorrgb = vec3(115.0/255.0, 161.0/255.0, 227.0/255.0);

    for (float r = 0; r < rone; r+= 0.1){
        for(float theta = 0; theta < 2 * pi; theta += incdegree){

            //pp
            vec3 pp = vec3(r * cos(theta),r * sin(theta),0);
            
            //rotation plane for rodrigues rotation
            //crossing z x eone
            vec3 plane = cross( z / length(z), eone / length(eone));
            //normalizing rotation plane
            plane = normalize(plane);

            //rotation angle calculation
            float rot = acos( dot(eone/length(eone),z / length(z)) );

            //rotation vector of pp + transform (l)
            vec3 Wpp =
                pp * cos(rot)  
                + normalize(cross(plane, 
                pp))
                + plane * (dot(plane,
                pp)) * (1-cos(rot)) 
                + l;
 
            //scale by area divided by number of samples
            //Light intensity calculation//pointlight[3] = inherent light intensity  
            float intense = (pi * (pow(rone,2))) / ( (rone / 0.1)  *   ((2*pi) / incdegree) );
            intense *= pointlight[3] * dot(-l / length(l)), eone / length(eone)); 

            //intensity is multiplied into equation
            //Wpp *= intense;
            
            //equation is multiplied into rgb values
            //colorrgb *= Wpp;
            //colorrgb = clamp(colorrgb, 0.0, 1.0);
            

            //diffuse highlight
            vec3 tmp = (normalize(vnormal));
            colorrgb *= dot(tmp , l / length(l));
            colorrgb /= length(l) * length(l);
            colorrgb = clamp(colorrgb, 0.0, 1.0);


            //phong specular highlight
            float rscale = 2.0 * dot( l / r, tmp);
            vec3 reflect = vec3(0.0, 0.0, 0.0);
            if (rscale > 0.0) {	reflect = rscale * tmp - (l / r);}
            float rdotv = clamp(dot(reflect, - vpos / length(vpos)), 0.0, 1.0);
            float specular = 0.1 * clamp(pow(rdotv , 5.0) * pointlight[3], 0.0, 1.0);
            specular /= r*r;
            colorrgb += specular;  

            //ambient component
            colorrgb += vec3(0.01, 0.01, 0.01);
    
            colorrgb = clamp(colorrgb, 0.0, 1.0);
            //add shaded region to bright
            bright += colorrgb;
        }
    }

    //outputted color vector is rgb vector + 1.0 (100%) alpha value
    color = vec4(bright, 1);
}